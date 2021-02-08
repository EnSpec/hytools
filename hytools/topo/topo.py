# -*- coding: utf-8 -*-
"""
Topographic correction

"""
import json
import numpy as np
import ray
from .modminn import apply_modminn,calc_modminn_coeffs
from .scsc import apply_scsc,calc_scsc_coeffs
from .cosine import apply_cosine,calc_cosine_coeffs
from .c import apply_c,calc_c_coeffs
from .scs import apply_scs,calc_scs_coeffs
from ..masks import mask_create

def calc_cosine_i(solar_zn, solar_az, aspect ,slope):
    """Generate cosine i image. The cosine of the incidence angle (i) is
       defined as the angle between the normal to the pixel surface
       and the solar zenith direction.
       All input geometry units must be in radians.

    Args:
        solar_az (numpy.ndarray): Solar azimuth angle.
        solar_zn (numpy.ndarray): Solar zenith angle.
        aspect (numpy.ndarray): Ground aspect.
        slope (numpy.ndarray): Ground slope.

    Returns:
        cnumpy.ndarray: Cosine i image.

    """

    relative_az = aspect - solar_az
    cosine_i = np.cos(solar_zn)*np.cos(slope) + np.sin(solar_zn)*np.sin(slope)*  np.cos(relative_az)
    return cosine_i

def apply_topo_correct(hy_obj,data,dimension,index):
    '''

    Args:
        hy_obj (TYPE): DESCRIPTION.
        band (TYPE): DESCRIPTION.
        index (TYPE): DESCRIPTION.

    Returns:
        band (TYPE): DESCRIPTION.

    '''

    if ('apply_topo' not in hy_obj.mask) & ('apply_mask' in hy_obj.topo):
        hy_obj.gen_mask(mask_create,'apply_topo',hy_obj.topo['apply_mask'])

    if hy_obj.topo['type'] == 'mod_minneart':
        data = apply_modminn(hy_obj,data,dimension,index)
    elif hy_obj.topo['type']  == 'scs+c':
        data = apply_scsc(hy_obj,data,dimension,index)
    elif hy_obj.topo['type']  == 'cosine':
        data = apply_cosine(hy_obj,data,dimension,index)
    elif hy_obj.topo['type']  == 'c':
        data = apply_c(hy_obj,data,dimension,index)
    elif hy_obj.topo['type']  == 'scs':
        data = apply_scs(hy_obj,data,dimension,index)
    return data

def load_topo_precomputed(hy_obj,topo_dict):
    with open(topo_dict['coeff_files'][hy_obj.file_name], 'r') as outfile:
        hy_obj.topo = json.load(outfile)

def calc_topo_coeffs(actors,topo_dict):

    if topo_dict['type'] == 'precomputed':
        print("Using precomputed topographic coefficients.")
        _ = ray.get([a.do.remote(load_topo_precomputed,topo_dict) for a in actors])

    else:
        print("Calculating topographic coefficients.")

        _ = ray.get([a.gen_mask.remote(mask_create,'calc_topo',
                                       topo_dict['calc_mask']) for a in actors])

        if topo_dict['type'] == 'scs+c':
            _ = ray.get([a.do.remote(calc_scsc_coeffs,topo_dict) for a in actors])

        elif topo_dict['type'] == 'scs':
            _ = ray.get([a.do.remote(calc_scs_coeffs,topo_dict) for a in actors])

        elif topo_dict['type'] == 'mod_minneart':
            _ = ray.get([a.do.remote(calc_modminn_coeffs,topo_dict) for a in actors])

        elif topo_dict['type'] == 'cosine':
            _ = ray.get([a.do.remote(calc_cosine_coeffs,topo_dict) for a in actors])

        elif topo_dict['type'] == 'c':
            _ = ray.get([a.do.remote(calc_c_coeffs,topo_dict) for a in actors])

    _ = ray.get([a.do.remote(lambda x: x.corrections.append('topo')) for a in actors])



