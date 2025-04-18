import json
import os
import warnings
import sys
import numpy as np

import hytools as ht
from hytools.io.envi import *
from hytools.masks import mask_dict

warnings.filterwarnings("ignore")

def main():

    config_file = sys.argv[1]
    image_order = int(sys.argv[2])
    trait_order = int(sys.argv[3])

    with open(config_file, 'r') as outfile:
        config_dict = json.load(outfile)

    image= config_dict["input_files"][image_order]

    actor = ht.HyTools()

    # Load data
    if config_dict['file_type'] in ('envi','emit','ncav'):
        anc_file = config_dict["anc_files"][image]
        if "glt_files" in config_dict:
            if bool(config_dict["glt_files"]):
                actor.read_file(image,config_dict['file_type'],anc_path=anc_file,glt_path=config_dict["glt_files"][image]) # chunk_glt writing is not supported
            else:
                actor.read_file(image,config_dict['file_type'],anc_path=anc_file)
        else:
            actor.read_file(image,config_dict['file_type'],anc_path=anc_file)

    elif config_dict['file_type'] == 'neon':
        actor.read_file(image,config_dict['file_type'])

    trait = config_dict['trait_models'][trait_order]

    with open(trait, 'r') as json_file:
        trait_model = json.load(json_file)
        print("\t %s" % trait_model["name"])

    apply_single_trait_models(actor,config_dict,trait_order)


def apply_single_trait_models(hy_obj,config_dict,trait_order):
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

        # Build trait image file
        header_dict = hy_obj.get_header()
        header_dict['wavelength'] = []
        header_dict['data ignore value'] = -9999
        header_dict['data type'] = 4
        header_dict['trait unit'] = trait_model['units']
        header_dict['band names'] = ["%s_mean" % trait_model["name"],
                                     "%s_std" % trait_model["name"],
                                     'range_mask'] + [mask[0] for mask in config_dict['masks']]
        header_dict['bands'] = len(header_dict['band names'] ) 

        #Generate masks
        for mask,args in config_dict['masks']:
            mask_function = mask_dict[mask]
            hy_obj.gen_mask(mask_function,mask,args)

        output_name = config_dict['output_dir']
        output_name += os.path.splitext(os.path.basename(hy_obj.file_name))[0] + "_%s" % trait_model["name"]

        writer = WriteENVI(output_name,header_dict)

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
                      chunk_size = (512,512),
                      corrections =  hy_obj.corrections,
                      resample=resample)

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


            # Subset and assign custom masks
            for i,(mask,args) in enumerate(config_dict['masks']):
                mask = hy_obj.mask[mask][iterator.current_line:iterator.current_line+chunk.shape[0],
                                              iterator.current_column:iterator.current_column+chunk.shape[1]]
                trait_est[:,:,3+i] = mask.astype(int)


            nd_mask = hy_obj.mask['no_data'][iterator.current_line:iterator.current_line+chunk.shape[0],
                                             iterator.current_column:iterator.current_column+chunk.shape[1]]
            trait_est[~nd_mask] = -9999
            writer.write_chunk(trait_est,
                               iterator.current_line,
                               iterator.current_column)

        writer.close()

if __name__== "__main__":
    main()
