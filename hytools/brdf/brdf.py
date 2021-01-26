# -*- coding: utf-8 -*-
"""
    BRDF Correction
"""
import ray
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
    elif hy_obj.brdf['type'] == 'local':
        data = apply_local(hy_obj,data,dimension,index)
    elif hy_obj.brdf['type'] == 'flex':
        data = apply_flex(hy_obj,data,dimension,index)
    return data

def load_brdf_precomputed(hy_obj,brdf_dict):
    with open(brdf_dict['coeff_files'][hy_obj.file_name], 'r') as outfile:
        hy_obj.brdf = json.load(outfile)

def brdf_coeffs(actors,brdf_dict):

    if brdf_dict['type'] == 'precomputed':
        print("Using precomputed BRDF coefficients")
        _ = ray.get([a.do.remote(load_brdf_precomputed,brdf_dict) for a in actors])

    else:
        print("Calculating BRDF coefficients")
        if brdf_dict['type']== 'universal':
            universal_brdf(actors,brdf_dict)
        elif brdf_dict['type'] == 'flex':
            flex_brdf(actors,brdf_dict)
        elif brdf_dict['type'] == 'local':
            local_brdf(actors,brdf_dict)

    _ = ray.get([a.do.remote(lambda x: x.corrections.append('brdf')) for a in actors])
