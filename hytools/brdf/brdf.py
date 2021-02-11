# -*- coding: utf-8 -*-
"""
    BRDF Correction
"""
import json
import ray
import numpy as np
from .universal import universal_brdf,apply_universal
from .flex import flex_brdf,apply_flex
from ..masks import mask_create
from ..misc import set_brdf, update_brdf

def apply_brdf_correct(hy_obj,data,dimension,index):
    ''' Apply in memory BRDF correction.
    '''

    if hy_obj.brdf['type'] == 'universal':
        data = apply_universal(hy_obj,data,dimension,index)
    elif hy_obj.brdf['type'] == 'flex':
        data = apply_flex(hy_obj,data,dimension,index)
    elif hy_obj.brdf['type'] == 'local':
        print('Local/class BRDF correction....under development')
    return data

def load_brdf_precomputed(hy_obj,brdf_dict):
    with open(brdf_dict['coeff_files'][hy_obj.file_name], 'r') as outfile:
        hy_obj.brdf = json.load(outfile)

def set_solar_zn(hy_obj):
    solar_zn = hy_obj.get_anc('solar_zn')
    solar_zn = np.mean(solar_zn[hy_obj.mask['no_data']])
    hy_obj.brdf['solar_zn_norm_radians'] = float(solar_zn)
    return solar_zn

def calc_brdf_coeffs(actors,config_dict):

    brdf_dict = config_dict['brdf']

    if brdf_dict['type'] == 'precomputed':
        print("Using precomputed BRDF coefficients")
        _ = ray.get([a.do.remote(load_brdf_precomputed,
                                 config_dict['brdf']) for a in actors])
    else:
        # Set BRDF dict
        _ = ray.get([a.do.remote(set_brdf,brdf_dict) for a in actors])

        # Create masks used for calculating coefficients
        _ = ray.get([a.gen_mask.remote(mask_create,'calc_brdf',
                                       brdf_dict['calc_mask']) for a in actors])
        # Calculate mean solar zenith
        if isinstance(brdf_dict['solar_zn_type'],str):

            # Assign per line mean solar zenith
            solar_zn_samples = ray.get([a.do.remote(set_solar_zn) for a in actors])
            # Calculate and assign scene average solar zenith
            if brdf_dict['solar_zn_type'] == 'scene':
                scene_mean = float(np.mean(solar_zn_samples))
                _ =  ray.get([a.do.remote(update_brdf,{'key':'solar_zn_norm_radians',
                                                      'value': scene_mean }) for a in actors])
                print("Scene average solar zenith angle : %s degrees" % round(np.degrees(scene_mean),3))

        elif isinstance(brdf_dict['solar_zn_type'],float):
            _ =  ray.get([a.do.remote(update_brdf,{'key':'solar_zn_norm_radians',
                                                   'value': brdf_dict['solar_zn_type']}) for a in actors])
        else:
            print('Unrecognized solar zenith angle normalization')

        print("Calculating BRDF coefficients")
        if brdf_dict['type']== 'universal':
            universal_brdf(actors,config_dict)
        elif brdf_dict['type'] == 'flex':
            flex_brdf(actors,config_dict)
        elif brdf_dict['type'] == 'local':
            print('Local/class BRDF correction....under development')

    _ = ray.get([a.do.remote(lambda x: x.corrections.append('brdf')) for a in actors])










