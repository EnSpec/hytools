# -*- coding: utf-8 -*-
"""
    BRDF Correction
"""
import ray
from .standard import apply_standard,global_brdf
from .dynamic import apply_class,class_brdf

def apply_brdf_correct(hy_obj,data,dimension,index):
    '''

    Args:
        hy_obj (TYPE): DESCRIPTION.
        band (TYPE): DESCRIPTION.
        index (TYPE): DESCRIPTION.

    Returns:
        band (TYPE): DESCRIPTION.

    '''
    if hy_obj.brdf['type'] == 'standard':
        data = apply_standard(hy_obj,data,dimension,index)
    elif hy_obj.brdf['type'] == 'class':
        data = apply_class(hy_obj,data,dimension,index)
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
        if brdf_dict['type']== 'standard':
            global_brdf(actors,brdf_dict)
        elif brdf_dict['type'] == 'class':
            class_brdf(actors,brdf_dict)

    _ = ray.get([a.do.remote(lambda x: x.corrections.append('brdf')) for a in actors])
