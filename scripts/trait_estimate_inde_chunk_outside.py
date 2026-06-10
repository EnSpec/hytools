import json
import os
import warnings
import sys
import ray
from ray.util.multiprocessing import Pool
import numpy as np
import time
import psutil

import hytools as ht
from hytools.io.envi import *
from hytools.io.netcdf import *
from hytools.masks import mask_dict
from hytools.glint import set_glint_parameters_single

warnings.filterwarnings("ignore")

def main():
    start_ram = read_current_ram_usage()
    time_start = time.perf_counter()

    config_file = sys.argv[1]
    image_order = int(sys.argv[2])
    trait_order = int(sys.argv[3])
    n_cores_outside = int(sys.argv[4])
    chunk_row_outside = sys.argv[5]
    chunk_col_outside = sys.argv[6]

    with open(config_file, 'r') as outfile:
        config_dict = json.load(outfile)

    if len(sys.argv)>7:
        meta_file = sys.argv[7]
        with open(meta_file, 'r') as outfile:
            meta_dict = json.load(outfile)
            config_dict["outside_metadata"] = meta_dict
    else:
        if "outside_metadata" in config_dict:
            if not isinstance(config_dict["outside_metadata"],dict):
                with open(config_dict["outside_metadata"], 'r') as outfile:  
                # load json and replace it by a dict
                    meta_dict = json.load(outfile)
                    config_dict["outside_metadata"] = meta_dict            
        else:
            config_dict["outside_metadata"] = None

    image= config_dict["input_files"][image_order]

    config_dict['num_cpus'] = int(n_cores_outside) # overwrite
    config_dict['chunk_size'] = [chunk_row_outside,chunk_col_outside]

    worker_count = min(os.cpu_count()-1,config_dict['num_cpus'])

    if worker_count>1:
        if ray.is_initialized():
            ray.shutdown()
        print("Using %s CPUs." % worker_count)
        context=ray.init(num_cpus = worker_count)
        print(context.dashboard_url)
    else:
        print("Using %s CPUs." % worker_count)

    HyTools = ht.HyTools()
    #actors = [HyTools.remote() for image in images]

    # Load data
    if config_dict['file_type'] in ('envi','emit','ncav'):
        if config_dict["anc_files"]:
            anc_file = config_dict["anc_files"][image]
        else:
            anc_file = {}

        if "glt_files" in config_dict:
            if bool(config_dict["glt_files"]):
            #print(config_dict["glt_files"])
            #return
                HyTools.read_file(image,config_dict['file_type'],anc_path=anc_file,glt_path=config_dict["glt_files"][image]) # chunk_glt writing is not supported
            else:
                HyTools.read_file(image,config_dict['file_type'],anc_path=anc_file)
        else:
            HyTools.read_file(image,config_dict['file_type'],anc_path=anc_file)
        #_ = ray.get([a.read_file.remote(image,config_dict['file_type'],
        #                                anc_files[image]) for a,image in zip(actors,images)])
    elif config_dict['file_type'] == 'neon':
        HyTools.read_file(image,config_dict['file_type'])
        #_ = ray.get([a.read_file.remote(image,config_dict['file_type']) for a,image in zip(actors,images)])

    ###############
    default_export_type = "envi"
    if "export_type" in config_dict:
        if not config_dict["export_type"] in ["envi","netcdf"]:
            print("Image export file type is not recognized.")
            return
    else:
        config_dict["export_type"]=default_export_type

    ###############

    trait = config_dict['trait_models'][trait_order]

    #for trait in config_dict['trait_models']:
    with open(trait, 'r') as json_file:
        trait_model = json.load(json_file)
        print("\t %s" % trait_model["name"])

    meta_result = apply_single_trait_models(HyTools,config_dict,trait_order,worker_count,start_ram_gb=start_ram)

    #_ = ray.get([a.do.remote(apply_trait_models,config_dict) for a in actors])
    if worker_count>1:
        ray.shutdown()

    time_end = time.perf_counter()
    print("Total Time: {} sec.\n=_=".format(time_end - time_start))  
    print(f"<<>@Total Time,{meta_result[0]},{meta_result[1]},{meta_result[2]},{meta_result[3]},{(time_end - time_start):.3f}\n^_^")

@ray.remote
class FileWriter:
    def __init__(self, outside_writer):
        self.writer = outside_writer

    def write_data(self, data,x,y):
        self.writer.write_chunk(data,
                    y,
                    x)
        return f"Wrote data chunk"

def read_current_ram_usage(print_flag=True):
    ram = psutil.virtual_memory()
    if print_flag:
        print("RAM usage (%):", ram.percent)
        print("RAM used (GB):", round(ram.used / 1e9, 2))
    return round(ram.used / 1e9, 2)

def calculate_einsum(inputs):
    ram_trajectory = []
    ram_trajectory+=[read_current_ram_usage(print_flag=False)]
    hy_obj, ind_list, coeffs,intercept, out_band_count, trait_transform, model_diagnostics,masks, iterator_corrections,iterator_resample, wave_mask = inputs
    chunk = hy_obj.get_chunk(ind_list[0],ind_list[1], ind_list[2],ind_list[3],
                                            corrections = iterator_corrections,
                                            resample = iterator_resample)
    # inputs can be a tuple of (subscripts, array1, array2) or similar
    #chunk,coeffs = inputs

    if not wave_mask is None:
        chunk = chunk[:,:,wave_mask]

    trait_est = np.zeros((chunk.shape[0],
                            chunk.shape[1],
                            out_band_count))

    # Apply spectrum transforms
    for transform in trait_transform:
        if  transform== "vector":    #vnorm
            norm = np.linalg.norm(chunk,axis=2)
            chunk = chunk/norm[:,:,np.newaxis]
        if transform == "absorb":
            chunk = np.log(1/chunk)
        if transform == "mean":
            mean = chunk.mean(axis=2)
            chunk = chunk/mean[:,:,np.newaxis]

    trait_pred = np.einsum('jkl,ml->jkm',chunk,coeffs, optimize='optimal')

    trait_pred = trait_pred + intercept
    trait_est[:,:,0] = trait_pred.mean(axis=2)
    trait_est[:,:,1] = trait_pred.std(ddof=1,axis=2)

    range_mask = (trait_est[:,:,0] > model_diagnostics['min']) & \
                (trait_est[:,:,0] < model_diagnostics['max'])
    trait_est[:,:,2] = range_mask.astype(int)

    ram_trajectory+=[read_current_ram_usage(print_flag=False)]
    # Subset and assign custom masks
    for i,(mask,args) in enumerate(masks):
        mask = hy_obj.mask[mask][ind_list[2]:ind_list[3],
                                    ind_list[0]:ind_list[1]]
        trait_est[:,:,3+i] = mask.astype(int)


    nd_mask = hy_obj.mask['no_data'][ind_list[2]:ind_list[3],
                                    ind_list[0]:ind_list[1]]
    trait_est[~nd_mask,:2] = -9999
    trait_est[~nd_mask,2:] = 255 #-9999

    ram_trajectory+=[read_current_ram_usage(print_flag=False)]
    return trait_est,ind_list[2],ind_list[0],ram_trajectory


    #return None #np.einsum('jkl,ml->jkm',chunk,coeffs, optimize='optimal')

def generate_chunk_ind(chunk_size,lines,columns):
    ind_list = []
    current_column = -1
    current_line = -1
    complete = False

    while not complete:
        if current_column == -1:
            current_column +=1
            current_line +=1
        else:
            current_column += chunk_size[1]
        if current_column >= columns:
            current_column = 0
            current_line += chunk_size[0]

        y_start = current_line
        y_end = current_line + chunk_size[0]
        if y_end >= lines:
            y_end = lines
        x_start = current_column
        x_end = current_column + chunk_size[1]
        if x_end >= columns:
            x_end = columns

        if (y_end == lines) and (x_end == columns):
            complete = True   

        ind_list+=[[x_start,x_end, y_start,y_end]]

    return ind_list

def apply_single_trait_models(hy_obj,config_dict,trait_order,worker_count,start_ram_gb=None):
    '''Apply trait model(s) to image and export to file.
    '''

    hy_obj.create_bad_bands(config_dict['bad_bands'])
    hy_obj.corrections  = config_dict['corrections']

    # Load correction coefficients
    if 'topo' in  hy_obj.corrections:
        hy_obj.load_coeffs(config_dict['topo'][hy_obj.file_name],'topo')

    if 'brdf' in hy_obj.corrections:
        hy_obj.load_coeffs(config_dict['brdf'][hy_obj.file_name],'brdf')

    if 'glint' in hy_obj.corrections:
        set_glint_parameters_single(hy_obj, config_dict)

    hy_obj.resampler['type'] = config_dict["resampling"]['type']

    for trait in [config_dict['trait_models'][trait_order]]:
        with open(trait, 'r') as json_file:
            trait_model = json.load(json_file)
            coeffs = np.array(trait_model['model']['coefficients'])
            intercept = np.array(trait_model['model']['intercepts'])
            model_waves = np.array(trait_model['wavelengths'])

        #Check if wavelengths match
        resample = not all(x in hy_obj.wavelengths for x in model_waves)

        if resample:
            hy_obj.resampler['out_waves'] = model_waves
            hy_obj.resampler['out_fwhm'] = trait_model['fwhm']
        else:
            wave_mask = [np.argwhere(x==hy_obj.wavelengths)[0][0] for x in model_waves]

        #############
        use_glt_output_bool=False
        if 'use_glt' in config_dict:
            use_glt_output_bool = config_dict['use_glt']
            if use_glt_output_bool is True:
                header_dict = hy_obj.get_header(warp_glt=True)
            else:
                header_dict = hy_obj.get_header()
        else:
            header_dict = hy_obj.get_header()
        ##############

        # Build trait image file
        header_dict['wavelength'] = []
        header_dict['data ignore value'] = -9999
        header_dict['data type'] = 4
        header_dict['trait unit'] = trait_model['units']
        header_dict['band names'] = ["%s_mean" % trait_model["name"],
                                     "%s_std" % trait_model["name"],
                                     'range_mask'] + [mask[0] for mask in config_dict['masks']]
        header_dict['bands'] = len(header_dict['band names'] )

        header_dict['file_type'] = config_dict['file_type']
        header_dict['transform'] = hy_obj.transform
        header_dict['projection'] = hy_obj.projection

        #Generate masks
        for mask,args in config_dict['masks']:
            mask_function = mask_dict[mask]
            hy_obj.gen_mask(mask_function,mask,args)

        output_name = config_dict['output_dir']

        ###
        if config_dict["export_type"]=="envi":
            output_name += os.path.splitext(os.path.basename(hy_obj.file_name))[0] + "_%s" % trait_model["name"]
            writer = WriteENVI(output_name,header_dict)
        else:
            output_name += os.path.splitext(os.path.basename(hy_obj.file_name))[0] + "_%s_chunk_test.nc" % trait_model["name"]
            header_dict['lines_glt'] = hy_obj.lines_glt
            header_dict['samples_glt'] = hy_obj.columns_glt
            writer = WriteNetCDF(output_name,header_dict, attr_dict=None, glt_bool=use_glt_output_bool, type_tag="trait", band_name=trait_model["name"])
            
            if use_glt_output_bool is False and config_dict['file_type'] == 'emit':
                writer.write_glt_dataset(hy_obj.glt_x,hy_obj.glt_y,dim_x_name="ortho_x",dim_y_name="ortho_y")

        ###

        time_before_loop = time.perf_counter()

        if "chunk_size" in config_dict:
            chunk_row, chunk_col = config_dict["chunk_size"]
            if not isinstance(chunk_row, int):
                if chunk_row =="row":
                    chunk_row = hy_obj.lines
            if not isinstance(chunk_col, int):
                if chunk_col =="col":
                    chunk_col = hy_obj.columns                
            iterator = hy_obj.iterate(by = 'chunk',
                      #chunk_size = (hy_obj.lines,16), #hy_obj.columns  32,32
                      chunk_size = (int(chunk_row),int(chunk_col)),
                      corrections =  hy_obj.corrections,
                      resample=resample)

        elif config_dict['file_type'] == 'envi' or config_dict['file_type'] == 'emit':

            if hy_obj.columns>hy_obj.lines:
                iterator = hy_obj.iterate(by = 'chunk',
                      #chunk_size = (hy_obj.lines,16), #hy_obj.columns  32,32
                      chunk_size = (8,hy_obj.columns),
                      corrections =  hy_obj.corrections,
                      resample=resample)
            else:
                iterator = hy_obj.iterate(by = 'chunk',
                      #chunk_size = (16,hy_obj.columns), #hy_obj.columns  32,32
                      chunk_size = (8,hy_obj.columns),
                      corrections =  hy_obj.corrections,
                      resample=resample)
        elif config_dict['file_type'] == 'neon':
            #print("Chunk size: ",math.ceil(hy_obj.lines/32),math.ceil(hy_obj.columns/32))
            if hy_obj.columns>hy_obj.lines:
                iterator = hy_obj.iterate(by = 'chunk',
                      #chunk_size = (int(np.ceil(hy_obj.lines/32)),int(np.ceil(hy_obj.columns/32))), #hy_obj.columns  32,32, math.ceil(hy_obj.lines/32),math.ceil(hy_obj.columns/32)  #math.ceil(hy_obj.lines/32)//2,math.ceil(hy_obj.columns/32)//2
                      chunk_size = (hy_obj.lines,int(np.ceil(hy_obj.columns/32))),
                      corrections =  hy_obj.corrections,
                      resample=resample)
            else:
                iterator = hy_obj.iterate(by = 'chunk',
                      chunk_size = (int(np.ceil(hy_obj.lines/32)),int(np.ceil(hy_obj.columns/32))),
                      corrections =  hy_obj.corrections,
                      resample=resample)
        elif config_dict['file_type'] == 'ncav':
            iterator = hy_obj.iterate(by = 'chunk',
                      chunk_size = (512,512),
                      corrections =  hy_obj.corrections,
                      resample=resample)

        #read_current_ram_usage()
        if worker_count==1:
            ##if config_dict["export_type"]=="netcdf":
            out_stack = np.zeros((header_dict['bands'],header_dict['lines'],header_dict['samples'])).astype(np.float32)

            ram_trajectory = []
            while not iterator.complete:
                chunk = iterator.read_next()
                if not resample:
                    chunk = chunk[:,:,wave_mask]

                trait_est = np.zeros((chunk.shape[0],
                                        chunk.shape[1],
                                        header_dict['bands']))

                # Apply spectrum transforms
                for transform in  trait_model['model']["transform"]:
                    if  transform== "vector":    #vnorm
                        norm = np.linalg.norm(chunk,axis=2)
                        chunk = chunk/norm[:,:,np.newaxis]
                    if transform == "absorb":
                        chunk = np.log(1/chunk)
                    if transform == "mean":
                        mean = chunk.mean(axis=2)
                        chunk = chunk/mean[:,:,np.newaxis]

                trait_pred = np.einsum('jkl,ml->jkm',chunk,coeffs, optimize='optimal')

                trait_pred = trait_pred + intercept
                trait_est[:,:,0] = trait_pred.mean(axis=2)
                trait_est[:,:,1] = trait_pred.std(ddof=1,axis=2)

                range_mask = (trait_est[:,:,0] > trait_model["model_diagnostics"]['min']) & \
                            (trait_est[:,:,0] < trait_model["model_diagnostics"]['max'])
                trait_est[:,:,2] = range_mask.astype(int)

                ram_trajectory+=[read_current_ram_usage(print_flag=False)]
                # Subset and assign custom masks
                for i,(mask,args) in enumerate(config_dict['masks']):
                    mask = hy_obj.mask[mask][iterator.current_line:iterator.current_line+chunk.shape[0],
                                                iterator.current_column:iterator.current_column+chunk.shape[1]]
                    #print('i:',i,header_dict['bands'])
                    trait_est[:,:,3+i] = mask.astype(int)


                nd_mask = hy_obj.mask['no_data'][iterator.current_line:iterator.current_line+chunk.shape[0],
                                                iterator.current_column:iterator.current_column+chunk.shape[1]]
                trait_est[~nd_mask,:2] = -9999
                trait_est[~nd_mask,2:] = 255 #-9999
                    
                ram_trajectory+=[read_current_ram_usage(print_flag=False)]

                x_start = iterator.current_column
                x_end = iterator.current_column + trait_est.shape[1]
                y_start = iterator.current_line
                y_end = iterator.current_line + trait_est.shape[0]    
                out_stack[:,y_start:y_end,x_start:x_end] = np.moveaxis(trait_est,-1,0)

                #if iterator.current_line>50: 
                #    break 
            max_ram = np.array(ram_trajectory).max()
            print("RAM max :",max_ram,round(max_ram-start_ram_gb,4))    
            print(f"^^^@RAM general usage,{hy_obj.interleave},{worker_count},{chunk_row},{chunk_col},{round(max_ram-start_ram_gb,4)}")      
        else:
            ##if config_dict["export_type"]=="netcdf":
            out_stack = np.zeros((header_dict['bands'],header_dict['lines'],header_dict['samples'])).astype(np.float32)
            #file_writer_actor = FileWriter.remote("output.txt")
            chunk_start_end_list = generate_chunk_ind(iterator.chunk_size,hy_obj.lines,hy_obj.columns)
 
            list_of_tasks= []
            for ind_list in chunk_start_end_list:
                if not resample:
                    list_of_tasks += [(hy_obj, ind_list, coeffs,intercept, header_dict['bands'],trait_model['model']["transform"],trait_model["model_diagnostics"],config_dict['masks'],iterator.corrections,iterator.resample,wave_mask)]
                else:
                    list_of_tasks += [(hy_obj, ind_list, coeffs,intercept, header_dict['bands'],trait_model['model']["transform"],trait_model["model_diagnostics"],config_dict['masks'],iterator.corrections,iterator.resample,None)]

            read_current_ram_usage()
            ram_trajectory_merge=[]
            with Pool() as p:
                results = p.map(calculate_einsum, list_of_tasks[:])
                read_current_ram_usage()
                for ii, each_result in enumerate(results):
                    trait_est, current_line, current_column, ram_trajectory = each_result
                    ram_trajectory_merge+=[ram_trajectory]
                    x_start = current_column
                    x_end = current_column + trait_est.shape[1]
                    y_start = current_line
                    y_end = current_line + trait_est.shape[0]    
                    out_stack[:,y_start:y_end,x_start:x_end] = np.moveaxis(trait_est,-1,0)

            max_ram = np.array(ram_trajectory_merge).max()
            print("RAM max :",max_ram,round(max_ram-start_ram_gb,4))
            print(f"^^^@RAM general usage,{hy_obj.interleave},{worker_count},{chunk_row},{chunk_col},{round(max_ram-start_ram_gb,4)}")

        if use_glt_output_bool:

            if config_dict["export_type"]=="envi":
                for iband in range(out_stack.shape[0]):
                    writer.write_band_glt(out_stack[iband,:,:],iband, (hy_obj.glt_y[hy_obj.fill_mask]-1,hy_obj.glt_x[hy_obj.fill_mask]-1),hy_obj.fill_mask)
                writer.close()
            else:
                for iband in range(2):
                    writer.write_netcdf_band_glt(out_stack[iband,:,:],iband, (hy_obj.glt_y[hy_obj.fill_mask]-1,hy_obj.glt_x[hy_obj.fill_mask]-1),hy_obj.fill_mask)
                writer.close()

                for iband in range(len(header_dict['band names'][2:])):
                    writer = WriteNetCDF(output_name,header_dict, attr_dict=config_dict["outside_metadata"], glt_bool=use_glt_output_bool, type_tag="mask", band_name=header_dict['band names'][2:][iband])
                    writer.write_mask_band_glt(out_stack[2+iband,:,:], (hy_obj.glt_y[hy_obj.fill_mask]-1,hy_obj.glt_x[hy_obj.fill_mask]-1),hy_obj.fill_mask)
                    writer.close()
        else:
            if config_dict["export_type"]=="netcdf":
                for iband in range(2):
                    writer.write_band(out_stack[iband,:,:],iband)
                writer.close()

                for iband in range(len(header_dict['band names'][2:])):
                    writer = WriteNetCDF(output_name,header_dict, attr_dict=config_dict["outside_metadata"], glt_bool=use_glt_output_bool, type_tag="mask", band_name=header_dict['band names'][2:][iband])
                    writer.write_mask_band(out_stack[2+iband,:,:])
                    writer.close()
            else:
                for iband in range(out_stack.shape[0]):
                    writer.write_band(out_stack[iband,:,:],iband)
                writer.close()

        time_after_loop = time.perf_counter()
        print("Total Time after header: {} sec.".format(time_after_loop - time_before_loop))
        print(f"@_@Total Time after header,{hy_obj.interleave},{worker_count},{chunk_row},{chunk_col},{(time_after_loop - time_before_loop):.3f}")
        return hy_obj.interleave,worker_count,chunk_row,chunk_col

if __name__== "__main__":
    main()
