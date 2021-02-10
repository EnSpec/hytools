import glob
import ray
import numpy as np
import json
import os
import warnings
import sys
import hytools as ht
from hytools.io.envi import *
from hytools.topo.scsc import *
from hytools.topo.modminn import *
from hytools.topo.cosine import *
from hytools.topo.c import *
from hytools.topo.scs import *
from hytools.brdf.standard import *
from hytools.masks import mask_dict

warnings.filterwarnings("ignore")

def main():

    config_file = sys.argv[1]

    with open(config_file, 'r') as outfile:
        config_dict = json.load(outfile)

    images= config_dict["input_files"]

    if ray.is_initialized():
        ray.shutdown()
    print("Using %s CPUs." % config_dict['num_cpus'])
    ray.init(num_cpus = config_dict['num_cpus'])

    HyTools = ray.remote(ht.HyTools)
    actors = [HyTools.remote() for image in images]

    if config_dict['file_type'] == 'envi':
        anc_files = config_dict["anc_files"]
        _ = ray.get([a.read_file.remote(image,config_dict['file_type'],anc_files[image]) for a,image in zip(actors,images)])
    elif config_dict['file_type'] == 'neon':
        _ = ray.get([a.read_file.remote(image,config_dict['file_type']) for a,image in zip(actors,images)])

    print("Estimating %s traits:" % len( config_dict['trait_models']))
    for trait in config_dict['trait_models']:
        print("\t %s" % os.path.splitext(os.path.basename(trait))[0])

    _ = ray.get([a.do.remote(apply_trait_models,config_dict) for a in actors])
    ray.shutdown()



def apply_trait_models(hy_obj,config_dict):
    '''Apply correction to image and export
        to file.

    Args:
        hy_obj (TYPE): DESCRIPTION.
        config_dict (TYPE): DESCRIPTION.
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
        else:
            wave_mask = [np.argwhere(x==hy_obj.wavelengths)[0][0] for x in model_waves]

        header_dict = hy_obj.get_header()
        header_dict['wavelength'] = []
        header_dict['data ignore value'] = 0
        header_dict['data type'] = 4
        header_dict['band names'] = ["%s_mean" % trait_model["trait_name"],
                                     "%s_std" % trait_model["trait_name"],
                                     'range_mask'] + [mask[0] for mask in config_dict['masks']]
        header_dict['bands'] = len(header_dict['band names'] )

        #Generate masks
        for mask,args in config_dict['masks']:
            mask_function = mask_dict[mask]
            hy_obj.gen_mask(mask_function,mask,args)

        output_name = config_dict['output_dir']
        output_name += os.path.splitext(os.path.basename(hy_obj.file_name))[0] + "_%s" % trait_model["trait_name"]

        writer = WriteENVI(output_name,header_dict)

        iterator = hy_obj.iterate(by = 'chunk',
                      chunk_size = (100,100),
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
                if  transform== "vnorm":
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

            range_mask = (trait_est[:,:,0] > trait_model["model_diagnostics"]['trait_min']) & \
                         (trait_est[:,:,0] < trait_model["model_diagnostics"]['trait_max'])
            trait_est[:,:,3] = range_mask.astype(int)

            # Subset and assign custom masks
            for i,(mask,args) in enumerate(config_dict['masks']):
                mask = hy_obj.mask[mask][iterator.current_line:iterator.current_line+chunk.shape[0],
                                              iterator.current_column:iterator.current_column+chunk.shape[1]]

                trait_est[:,:,3+i] = mask.astype(int)


            nd_mask = hy_obj.mask['no_data'][iterator.current_line:iterator.current_line+chunk.shape[0],
                                             iterator.current_column:iterator.current_column+chunk.shape[1]]
            trait_est[~nd_mask] = 0
            writer.write_chunk(trait_est,
                               iterator.current_line,
                               iterator.current_column)
        writer.close()



if __name__== "__main__":
    main()




















