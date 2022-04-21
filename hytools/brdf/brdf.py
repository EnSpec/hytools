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
    hy_obj.brdf['norm_solar_zn'] = float(solar_zn)
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

        solar_zn_samples = ray.get([a.do.remote(set_solar_zn) for a in actors])

        #Calculate and assign normalization angles
        if 'norm_solar_zn' in brdf_dict.keys():
            if brdf_dict['norm_solar_zn'] == 'scene':
                norm_angle = float(np.mean(solar_zn_samples))
                _ =  ray.get([a.do.remote(update_brdf,{'key':'norm_solar_zn',
                           'value': norm_angle}) for a in actors])
                print("Scene average solar zenith angle : %s degrees" % round(np.degrees(norm_angle),3))
            elif isinstance(brdf_dict['norm_solar_zn'],float):
                norm_angle = brdf_dict['norm_solar_zn']
                _ =  ray.get([a.do.remote(update_brdf,{'key':'norm_solar_zn',
                           'value': norm_angle}) for a in actors])
            elif brdf_dict['norm_solar_zn'] == 'line':
                print('Using solar zenith line mean')
            else:
                print('Unrecognized solar zenith angle normalization, using line mean')
        else:
            print("solar_zn normalization angle not set, using line mean.")


        for angle in ['sensor_az','sensor_zn','solar_az']:
            if "norm_%s" % angle in brdf_dict.keys():
                norm_angle = brdf_dict["norm_%s" % angle]
            else:
                print("%s normalization angle not set, using 0." % angle )
                norm_angle = 0
            _ =  ray.get([a.do.remote(update_brdf,{'key':"norm_%s" % angle,
                                                       'value': norm_angle}) for a in actors])

        print("Calculating BRDF coefficients")
        if brdf_dict['type']== 'universal':
            universal_brdf(actors,config_dict)
        elif brdf_dict['type'] == 'flex':
            flex_brdf(actors,config_dict)
        elif brdf_dict['type'] == 'local':
            print('Local/class BRDF correction....under development')

    _ = ray.get([a.do.remote(lambda x: x.corrections.append('brdf')) for a in actors])
