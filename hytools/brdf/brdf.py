# -*- coding: utf-8 -*-
"""
    BRDF Correction
"""
import json
import ray
import numpy as np
from .universal import apply_universal,universal_brdf
from .flex import apply_flex,flex_brdf

def apply_brdf_correct(hy_obj,data,dimension,index):
    '''
    Args:
        hy_obj (TYPE): DESCRIPTION.
        band (TYPE): DESCRIPTION.
        index (TYPE): DESCRIPTION.

    Returns:
        band (TYPE): DESCRIPTION.

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

def get_solar_zn(hy_obj):
    solar_zn = hy_obj.get_anc('solar_zn')
    return np.mean(solar_zn[hy_obj.mask['no_data']])

def brdf_coeffs(actors,config_dict):
    brdf_dict = config_dict['brdf']

    if brdf_dict['type'] == 'precomputed':
        print("Using precomputed BRDF coefficients")
        _ = ray.get([a.do.remote(load_brdf_precomputed,config_dict['brdf']) for a in actors])

    else:
        #Calculate mean solar zenith
        if brdf_dict['solar_zn_norm']:
            solar_zn_samples = ray.get([a.do.remote(get_solar_zn) for a in actors])
            config_dict['brdf']['mean_solar_zenith'] = float(np.mean(solar_zn_samples))

        print("Calculating BRDF coefficients")
        if brdf_dict['type']== 'universal':
            universal_brdf(actors,config_dict)
        elif brdf_dict['type'] == 'flex':
            flex_brdf(actors,config_dict)
        elif brdf_dict['type'] == 'local':
            print('Local/class BRDF correction....under development')

    _ = ray.get([a.do.remote(lambda x: x.corrections.append('brdf')) for a in actors])
