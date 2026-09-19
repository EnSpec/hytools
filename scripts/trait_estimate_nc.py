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
from hytools.glint import set_glint_parameters_single

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

    pass_bool, anc_required_bool = check_anc_requirement(config_dict)
    if not pass_bool:
        return

    if ray.is_initialized():
        ray.shutdown()
    print("Using %s CPUs." % config_dict['num_cpus'])
    ray.init(num_cpus = config_dict['num_cpus'])

    HyTools = ray.remote(ht.HyTools)
    actors = [HyTools.remote() for image in images]

    # Load data
    if config_dict['file_type'] == 'envi':
        if bool(config_dict["glt_files"]):
            glt_files = config_dict["glt_files"]
            if anc_required_bool:
                _ = ray.get([a.read_file.remote(image,config_dict['file_type'],
                                        anc_path=anc_files[image],
                                        glt_path=glt_files[image]) for a,image in zip(actors,images)])
            else:
                _ = ray.get([a.read_file.remote(image,config_dict['file_type'],
                                        glt_path=glt_files[image]) for a,image in zip(actors,images)])
        else:
            if anc_required_bool:
                _ = ray.get([a.read_file.remote(image,config_dict['file_type'],
                                        anc_path=anc_files[image]) for a,image in zip(actors,images)])
            else:
                _ = ray.get([a.read_file.remote(image,config_dict['file_type']) for a,image in zip(actors,images)])

    elif config_dict['file_type'] == 'neon':
        _ = ray.get([a.read_file.remote(image,config_dict['file_type']) for a,image in zip(actors,images)])
    elif config_dict['file_type'] in ['emit','ncav','tanager']:
        anc_files = config_dict["anc_files"]
        if bool(config_dict["glt_files"]):
            glt_files = config_dict["glt_files"]
            if anc_required_bool:
                _ = ray.get([a.read_file.remote(image,config_dict['file_type'],
                                        anc_path=anc_files[image],
                                        glt_path=glt_files[image]) for a,image in zip(actors,images)])
            else:
                _ = ray.get([a.read_file.remote(image,config_dict['file_type'],
                                        glt_path=glt_files[image]) for a,image in zip(actors,images)])

        else:
            if anc_required_bool:
                _ = ray.get([a.read_file.remote(image,config_dict['file_type'],
                                        anc_path=anc_files[image]) for a,image in zip(actors,images)])
            else:
                _ = ray.get([a.read_file.remote(image,config_dict['file_type']) for a,image in zip(actors,images)])

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

    print(f"Estimating {len( config_dict['trait_models'])} traits:")
    for trait in config_dict['trait_models']:
        with open(trait, 'r') as json_file:
            trait_model = json.load(json_file)
            print("\t %s" % trait_model["name"])

    _ = ray.get([a.do.remote(apply_trait_models,config_dict) for a in actors])
    ray.shutdown()

def check_anc_requirement(config_dict):
    '''Check if ANC files are required and provided in the config file, if they are required.
    '''
    anc_files = config_dict["anc_files"]
    pass_bool=False
    if ('topo' in config_dict['corrections']) or ('brdf' in config_dict['corrections']):
        anc_required_bool=True
        if config_dict['file_type'] in ['envi','emit','ncav']:
            if bool(anc_files):
                pass_bool=True
            else:
                print("'anc' files are required for correction, but they are not provided.")
                #pass_bool=False # default
        elif config_dict['file_type'] in ['neon']:
            pass_bool=True
        elif config_dict['file_type'] in ['tanager']:
            if 'topo' in config_dict['corrections']:
                if bool(anc_files):
                    pass_bool=True
                else:
                    print("External 'anc' files with slope and aspect are required for TOPO correction, but they are not provided.")
            else:
                pass_bool=True
        else:
            # not accepted image format
            pass # do not pass, pass_bool is still False
    else:
        anc_required_bool=False
        pass_bool=True

    return pass_bool, anc_required_bool

def apply_trait_models(hy_obj,config_dict):
    '''Apply trait model(s) to image and export to file.

    Models that share a wavelength set are applied in a single pass over the image
    (one read, correction and resampling of the image instead of one per model), and
    models that also share a spectrum transform chain are stacked into one coefficient
    matrix, so each chunk is transformed once and multiplied once. The output files and
    values are unchanged: one ENVI or NetCDF file per model.

    Every model of a pass keeps its own output stack (bands x lines x samples, float32)
    in memory until the pass ends; the optional config key "models_per_pass" caps how
    many models share a pass when that does not fit (default: all).
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

    #Generate masks once, they do not depend on the model
    for mask,args in config_dict['masks']:
        mask_function = mask_dict[mask]
        hy_obj.gen_mask(mask_function,mask,args)

    # Load the models and group them by wavelength set: one image pass per group
    groups = {}
    for trait in config_dict['trait_models']:
        with open(trait, 'r') as json_file:
            trait_model = json.load(json_file)
        key = (tuple(trait_model['wavelengths']),tuple(trait_model['fwhm']))
        groups.setdefault(key,[]).append(trait_model)

    models_per_pass = config_dict.get("models_per_pass",0)
    for (model_waves,model_fwhm),trait_models in groups.items():
        if models_per_pass and models_per_pass > 0:
            batches = [trait_models[i:i+models_per_pass] for i in range(0,len(trait_models),models_per_pass)]
        else:
            batches = [trait_models]
        for batch in batches:
            apply_model_group(hy_obj,config_dict,np.array(model_waves),list(model_fwhm),batch)

def apply_model_group(hy_obj,config_dict,model_waves,model_fwhm,trait_models):
    '''Apply every model in trait_models (all on model_waves) in one pass over the image.
    '''

    #Check if wavelengths match
    resample = not all(x in hy_obj.wavelengths for x in model_waves)

    if resample:
        hy_obj.resampler['out_waves'] = model_waves
        hy_obj.resampler['out_fwhm'] = model_fwhm
    else:
        wave_mask = [np.argwhere(x==hy_obj.wavelengths)[0][0] for x in model_waves]

    use_glt_output_bool=False
    if 'use_glt' in config_dict:
        use_glt_output_bool = config_dict['use_glt']
        if use_glt_output_bool==True:
            base_header = hy_obj.get_header(warp_glt=True)
        else:
            base_header = hy_obj.get_header()
    else:
        base_header = hy_obj.get_header()

    # Build trait image headers, one per model
    base_header['wavelength'] = []
    base_header['data ignore value'] = -9999
    base_header['data type'] = 4
    base_header['file_type'] = config_dict['file_type']
    base_header['transform'] = hy_obj.transform
    base_header['projection'] = hy_obj.projection
    n_bands = 3 + len(config_dict['masks'])

    headers = {}
    for trait_model in trait_models:
        header_dict = dict(base_header)
        header_dict['trait unit'] = trait_model['units']
        header_dict['band names'] = ["%s_mean" % trait_model["name"],
                                     "%s_std" % trait_model["name"],
                                     'range_mask'] + [mask[0] for mask in config_dict['masks']]
        header_dict['bands'] = n_bands
        headers[trait_model["name"]] = header_dict

    # Models with the same transform chain share the transformed chunk and one
    # matrix product: their coefficients are stacked along the iteration axis.
    stacks = {}
    for trait_model in trait_models:
        stacks.setdefault(tuple(trait_model['model']['transform']),[]).append(trait_model)
    stacked = []
    for transforms,members in stacks.items():
        coeffs = np.concatenate([np.array(m['model']['coefficients']) for m in members])
        intercept = np.concatenate([np.array(m['model']['intercepts']) for m in members])
        sizes = [len(m['model']['intercepts']) for m in members]
        offsets = np.cumsum([0] + sizes)
        stacked.append((transforms,members,coeffs,intercept,offsets))

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

    elif config_dict['file_type'] == 'tanager':
        iterator = hy_obj.iterate(by = 'chunk',
                  chunk_size= (int(np.ceil(hy_obj.lines/16)),int(np.ceil(hy_obj.columns/16))),
                  corrections =  hy_obj.corrections,
                  resample=resample)

    out_stacks = {}
    for trait_model in trait_models:
        out_stacks[trait_model["name"]] = np.zeros((n_bands,base_header['lines'],base_header['samples'])).astype(np.float32)

    while not iterator.complete:
        chunk = iterator.read_next()
        if not resample:
            chunk = chunk[:,:,wave_mask]
        lines = chunk.shape[0]
        columns = chunk.shape[1]

        # Mask bands and the no-data mask are the same for every model
        mask_bands = np.zeros((lines,columns,len(config_dict['masks'])))
        for i,(mask,args) in enumerate(config_dict['masks']):
            mask_bands[:,:,i] = hy_obj.mask[mask][iterator.current_line:iterator.current_line+lines,
                                                  iterator.current_column:iterator.current_column+columns]
        nd_mask = hy_obj.mask['no_data'][iterator.current_line:iterator.current_line+lines,
                                         iterator.current_column:iterator.current_column+columns]

        x_start = iterator.current_column
        x_end = iterator.current_column + columns
        y_start = iterator.current_line
        y_end = iterator.current_line + lines

        for transforms,members,coeffs,intercept,offsets in stacked:
            # Apply spectrum transforms
            transformed = chunk
            for transform in transforms:
                if  transform== "vector":
                    norm = np.linalg.norm(transformed,axis=2)
                    transformed = transformed/norm[:,:,np.newaxis]
                if transform == "absorb":
                    transformed = np.log(1/transformed)
                if transform == "mean":
                    mean = transformed.mean(axis=2)
                    transformed = transformed/mean[:,:,np.newaxis]

            # One product for every iteration of every model in the stack
            trait_pred = np.einsum('jkl,ml->jkm',transformed,coeffs, optimize='optimal')
            trait_pred = trait_pred + intercept

            for trait_model,start,stop in zip(members,offsets[:-1],offsets[1:]):
                pred = trait_pred[:,:,start:stop]
                trait_est = np.zeros((lines,columns,n_bands))
                trait_est[:,:,0] = pred.mean(axis=2)
                trait_est[:,:,1] = pred.std(ddof=1,axis=2)

                range_mask = (trait_est[:,:,0] > trait_model["model_diagnostics"]['min']) & \
                             (trait_est[:,:,0] < trait_model["model_diagnostics"]['max'])
                trait_est[:,:,2] = range_mask.astype(int)
                trait_est[:,:,3:] = mask_bands

                trait_est[~nd_mask,:2] = -9999
                trait_est[~nd_mask,2:] = 255

                out_stacks[trait_model["name"]][:,y_start:y_end,x_start:x_end] = np.moveaxis(trait_est,-1,0)

    for trait_model in trait_models:
        export_trait(hy_obj,config_dict,trait_model,headers[trait_model["name"]],
                     out_stacks.pop(trait_model["name"]),use_glt_output_bool)

def export_trait(hy_obj,config_dict,trait_model,header_dict,out_stack,use_glt_output_bool):
    '''Write one model's output stack as an ENVI or NetCDF file.
    '''
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
                writer = WriteNetCDF(output_name,header_dict,
                                     attr_dict=config_dict["outside_metadata"],
                                     glt_bool=use_glt_output_bool,
                                     type_tag="mask",
                                     band_name=header_dict['band names'][2:][iband])
                writer.write_mask_band_glt(out_stack[2+iband,:,:], (hy_obj.glt_y[hy_obj.fill_mask]-1,hy_obj.glt_x[hy_obj.fill_mask]-1),hy_obj.fill_mask)
                writer.close()
    else:
        if config_dict["export_type"]=="envi":
            for iband in range(out_stack.shape[0]):
                writer.write_band(out_stack[iband,:,:],iband)
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
