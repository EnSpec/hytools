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

    trait_models = {}

    for trait in config_dict['trait_models']:
        trait_name =os.path.splitext(os.path.basename(trait))[0]
        with open(trait, 'r') as json_file:
            trait_model = json.load(json_file)
            trait_model['coefficients'] = np.array(trait_model['coefficients'])
            trait_model['intercept'] = np.array(trait_model['intercept'])
            trait_models[trait_name]  =trait_model
            
    hy_obj.resampler['type'] = 'cubic'        
    hy_obj.resampler['out_waves'] = trait_models[trait_name]['model_wavelengths']       
            
    iterator = hy_obj.iterate(by = 'chunk',
                  chunk_size = (100,100),
                  corrections =  hy_obj.corrections,
                  resample=True)
    
    header_dict = hy_obj.get_header()
    header_dict['bands'] = len(config_dict['trait_models'])
    header_dict['wavelength'] = []
    header_dict['data ignore value'] = 0
    header_dict['data type'] = 4
    output_name = config_dict['output_dir'] 
    output_name += os.path.splitext(os.path.basename(hy_obj.file_name))[0] + "_traits"
    
    writer = WriteENVI(output_name,header_dict)
    
    no_data = hy_obj.get_band(0) == hy_obj.no_data

    while not iterator.complete:  
        chunk = iterator.read_next() 
                
        trait_mean = np.zeros((chunk.shape[0],
                                chunk.shape[1],
                                len(config_dict['trait_models'])))
        
        mask = no_data[iterator.current_line:iterator.current_line+chunk.shape[0],
                       iterator.current_column:iterator.current_column+chunk.shape[1]]
    
        if trait_models[trait_name]['vector_norm']:
            chunk /= np.linalg.norm(chunk,axis = 2)[:,:,np.newaxis]
            chunk *= trait_models[trait_name]['vector_scaler']
        
        for t,trait in enumerate(trait_models):
            trait_model = trait_models[trait]
            trait_pred = np.einsum('jkl,ml->jkm',chunk,trait_model['coefficients'], optimize='optimal')
            trait_pred = trait_pred + trait_model['intercept']
            trait_mean[:,:,t] = trait_pred.mean(axis=2)
        
        trait_mean[mask] = 0
        writer.write_chunk(trait_mean,
                           iterator.current_line,
                           iterator.current_column)
    writer.close()



if __name__== "__main__":
    main()




















