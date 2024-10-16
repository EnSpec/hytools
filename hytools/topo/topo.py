# -*- coding: utf-8 -*-

"""
HyTools:  Hyperspectral image processing library
Copyright (C) 2021 University of Wisconsin

Authors: Adam Chlus, Zhiwei Ye, Philip Townsend.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

Topographic correction

"""
import json
import numpy as np
import ray
from .modminn import apply_modminn,calc_modminn_coeffs
from .scsc import apply_scsc,calc_scsc_coeffs, calc_scsc_coeffs_group
from .cosine import apply_cosine,calc_cosine_coeffs
from .c import apply_c,calc_c_coeffs, calc_c_coeffs_group
from .scs import apply_scs,calc_scs_coeffs
from ..masks import mask_create
from ..misc import set_topo

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

def get_topo_sample_mask(hy_obj,topo_dict):

    sample_ratio = float(topo_dict["sample_perc"])   
    
    subsample_mask = np.copy(hy_obj.mask['calc_topo'])

    idx = np.array(np.where(subsample_mask!=0)).T

    idxRand= idx[np.random.choice(range(len(idx)),int(len(idx)*(1-sample_ratio)), replace = False)].T

    subsample_mask[idxRand[0],idxRand[1]] = 0
    subsample_mask = subsample_mask.astype(np.int8)

    hy_obj.ancillary['sample_mask']=subsample_mask

def calc_topo_coeffs(actors,topo_dict,actor_group_list=None,group_tag_list=None):
#def calc_topo_coeffs(actors,actor_group_list,topo_dict,group_tag_list):
    if topo_dict['type'] == 'precomputed':
        print("Using precomputed topographic coefficients.")
        _ = ray.get([a.do.remote(load_topo_precomputed,topo_dict) for a in actors]) # actors

        #_ = ray.get([a.do.remote(lambda x: x.corrections.append('topo')) for a in actors])

    else:
        print("Calculating topographic coefficients.")

        _ = ray.get([a.do.remote(set_topo,topo_dict) for a in actors])

        _ = ray.get([a.gen_mask.remote(mask_create,'calc_topo',
                                       topo_dict['calc_mask']) for a in actors])

        if (actor_group_list is None) or (topo_dict['type'] in ['scs','mod_minneart','cosine']):
        # no grouping
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

            #_ = ray.get([a.do.remote(lambda x: x.corrections.append('topo')) for a in actors])

        else:
            
            _ = ray.get([a.do.remote(get_topo_sample_mask,topo_dict) for a in actors])

            for group_order, sub_actors in enumerate(actor_group_list):

                #return 0

                if topo_dict['type'] == 'scs+c':
                    calc_scsc_coeffs_group(sub_actors,topo_dict,group_tag_list[group_order])

                elif topo_dict['type'] == 'c':
                    calc_c_coeffs_group(sub_actors,topo_dict,group_tag_list[group_order])

    _ = ray.get([a.do.remote(lambda x: x.corrections.append('topo')) for a in actors])


def calc_topo_coeffs_single(hy_obj,topo_dict):

    if topo_dict['type'] == 'precomputed':
        print("Using precomputed topographic coefficients.")
        load_topo_precomputed(hy_obj,topo_dict)

    else:
        print("Calculating topographic coefficients.")

        hy_obj.gen_mask(mask_create,'calc_topo',topo_dict['calc_mask'])


        if topo_dict['type'] == 'scs+c':
            calc_scsc_coeffs(hy_obj,topo_dict)

        elif topo_dict['type'] == 'scs':
            calc_scs_coeffs(hy_obj,topo_dict)

        elif topo_dict['type'] == 'mod_minneart':
            calc_modminn_coeffs(hy_obj,topo_dict)

        elif topo_dict['type'] == 'cosine':
            calc_cosine_coeffs(hy_obj,topo_dict)

        elif topo_dict['type'] == 'c':
            calc_c_coeffs(hy_obj,topo_dict)
    
    hy_obj.corrections.append('topo')

