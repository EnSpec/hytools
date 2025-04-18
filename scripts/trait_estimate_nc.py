import json
import os
import warnings
import sys
import ray
import numpy as np

import hytools as ht
from hytools.io.envi import *
from hytools.io.netcdf import *
from hytools.masks import mask_dict

warnings.filterwarnings("ignore")

def main():

    config_file = sys.argv[1]

    with open(config_file, 'r') as outfile:
        config_dict = json.load(outfile)

    if len(sys.argv)>2:
        meta_file = sys.argv[2]
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

    images= config_dict["input_files"]

    if ray.is_initialized():
        ray.shutdown()
    print("Using %s CPUs." % config_dict['num_cpus'])
    ray.init(num_cpus = config_dict['num_cpus'])

    HyTools = ray.remote(ht.HyTools)
    actors = [HyTools.remote() for image in images]

    # Load data
    if config_dict['file_type'] == 'envi':
        anc_files = config_dict["anc_files"]
        _ = ray.get([a.read_file.remote(image,config_dict['file_type'],
                                        anc_files[image]) for a,image in zip(actors,images)])
    elif config_dict['file_type'] == 'neon':
        _ = ray.get([a.read_file.remote(image,config_dict['file_type']) for a,image in zip(actors,images)])
    elif config_dict['file_type'] == 'emit' or config_dict['file_type'] == 'ncav':
        anc_files = config_dict["anc_files"]
        if bool(config_dict["glt_files"]):
            glt_files = config_dict["glt_files"]
            _ = ray.get([a.read_file.remote(image,config_dict['file_type'],
                                        anc_path=anc_files[image],glt_path=glt_files[image]) for a,image in zip(actors,images)])
        else:
            _ = ray.get([a.read_file.remote(image,config_dict['file_type'],
                                        anc_path=anc_files[image]) for a,image in zip(actors,images)])
    else:
        print("Image file type is not recognized.")
        return


    default_export_type = "envi"
    if "export_type" in config_dict:
        if not config_dict["export_type"] in ["envi","netcdf"]:
            print("Image export file type is not recognized.")
            return
    else:
        config_dict["export_type"]=default_export_type

    if not "use_glt" in config_dict:
        config_dict["use_glt"]=False

    print("Estimating %s traits:" % len( config_dict['trait_models']))
    for trait in config_dict['trait_models']:
        with open(trait, 'r') as json_file:
            trait_model = json.load(json_file)
            print("\t %s" % trait_model["name"])

    _ = ray.get([a.do.remote(apply_trait_models,config_dict) for a in actors])
    ray.shutdown()

def apply_trait_models(hy_obj,config_dict):
    '''Apply trait model(s) to image and export to file.
    '''

    hy_obj.create_bad_bands(config_dict['bad_bands'])
    hy_obj.corrections  = config_dict['corrections']

    # Load correction coefficients
    if 'topo' in  hy_obj.corrections:
        hy_obj.load_coeffs(config_dict['topo'][hy_obj.file_name],'topo')

    if 'brdf' in hy_obj.corrections:
        hy_obj.load_coeffs(config_dict['brdf'][hy_obj.file_name],'brdf')

    hy_obj.resampler['type'] = config_dict["resampling"]['type']

    for trait in config_dict['trait_models']:
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


        use_glt_output_bool=False
        if 'use_glt' in config_dict:
            use_glt_output_bool = config_dict['use_glt']
            if use_glt_output_bool==True:
                header_dict = hy_obj.get_header(warp_glt=True)
            else:
                header_dict = hy_obj.get_header()
        else:
            header_dict = hy_obj.get_header()

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

        if config_dict["export_type"]=="envi":
            output_name += os.path.splitext(os.path.basename(hy_obj.file_name))[0] + "_%s" % trait_model["name"]
            writer = WriteENVI(output_name,header_dict)
        else:
            output_name += os.path.splitext(os.path.basename(hy_obj.file_name))[0] + "_%s.nc" % trait_model["name"]
            header_dict['lines_glt'] = hy_obj.lines_glt
            header_dict['samples_glt'] = hy_obj.columns_glt
            writer = WriteNetCDF(output_name,header_dict,
                                 attr_dict=None,
                                 glt_bool=use_glt_output_bool,
                                 type_tag="trait",
                                 band_name=trait_model["name"])

            if (not use_glt_output_bool) and config_dict['file_type'] == 'emit':
                writer.write_glt_dataset(hy_obj.glt_x,hy_obj.glt_y,dim_x_name="ortho_x",dim_y_name="ortho_y")

        if config_dict['file_type'] == 'envi' or config_dict['file_type'] == 'emit':
            iterator = hy_obj.iterate(by = 'chunk',
                      chunk_size = (2,hy_obj.columns),
                      corrections =  hy_obj.corrections,
                      resample=resample)
        elif config_dict['file_type'] == 'neon':
            iterator = hy_obj.iterate(by = 'chunk',
                      chunk_size = (int(np.ceil(hy_obj.lines/32)),int(np.ceil(hy_obj.columns/32))),
                      corrections =  hy_obj.corrections,
                      resample=resample)

        elif config_dict['file_type'] == 'ncav':

            iterator = hy_obj.iterate(by = 'chunk',
                      chunk_size = (256,hy_obj.columns),
                      corrections =  hy_obj.corrections,
                      resample=resample)

        out_stack = np.zeros((header_dict['bands'],header_dict['lines'],header_dict['samples'])).astype(np.float32)

        while not iterator.complete:
            chunk = iterator.read_next()
            if not resample:
                chunk = chunk[:,:,wave_mask]

            trait_est = np.zeros((chunk.shape[0],
                                    chunk.shape[1],
                                    header_dict['bands']))

            # Apply spectrum transforms
            for transform in  trait_model['model']["transform"]:
                if  transform== "vector":
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


            # Subset and assign custom masks
            for i,(mask,args) in enumerate(config_dict['masks']):
                mask = hy_obj.mask[mask][iterator.current_line:iterator.current_line+chunk.shape[0],
                                              iterator.current_column:iterator.current_column+chunk.shape[1]]

                trait_est[:,:,3+i] = mask.astype(int)


            nd_mask = hy_obj.mask['no_data'][iterator.current_line:iterator.current_line+chunk.shape[0],
                                             iterator.current_column:iterator.current_column+chunk.shape[1]]

            trait_est[~nd_mask,:2] = -9999
            trait_est[~nd_mask,2:] = 255

            x_start = iterator.current_column
            x_end = iterator.current_column + trait_est.shape[1]
            y_start = iterator.current_line
            y_end = iterator.current_line + trait_est.shape[0]            
            out_stack[:,y_start:y_end,x_start:x_end] = np.moveaxis(trait_est,-1,0)

        if use_glt_output_bool:
            for iband in range(2):
                writer.write_netcdf_band_glt(out_stack[iband,:,:],iband, (hy_obj.glt_y[hy_obj.fill_mask]-1,hy_obj.glt_x[hy_obj.fill_mask]-1),hy_obj.fill_mask)
            writer.close()


            for iband in range(len(header_dict['band names'][2:])):
                writer = WriteNetCDF(output_name,header_dict,
                                     attr_dict=config_dict["outside_metadata"],
                                     glt_bool=use_glt_output_bool,
                                     type_tag="mask",
                                     band_name=header_dict['band names'][2:][iband])
                writer.write_mask_band_glt(out_stack[2+iband,:,:], (hy_obj.glt_y[hy_obj.fill_mask]-1,hy_obj.glt_x[hy_obj.fill_mask]-1),hy_obj.fill_mask)
                writer.close()
        else:
            for iband in range(2):
                writer.write_band(out_stack[iband,:,:],iband)
            writer.close()

            for iband in range(len(header_dict['band names'][2:])):
                writer = WriteNetCDF(output_name,header_dict,
                                     attr_dict=config_dict["outside_metadata"],
                                     glt_bool=use_glt_output_bool,
                                     type_tag="mask",
                                     band_name=header_dict['band names'][2:][iband])
                writer.write_mask_band(out_stack[2+iband,:,:])
                writer.close()


if __name__== "__main__":
    main()
