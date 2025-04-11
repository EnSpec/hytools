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
import h5py
from .universal import universal_brdf,apply_universal
from .flex import flex_brdf,apply_flex,ndvi_stratify, get_kernel_samples, ndvi_bins, get_band_samples
from ..masks import mask_create
from ..misc import set_brdf, update_brdf, progbar

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

def ndvi_stratify_samples(combine_dict):
    '''Create NDVI bin stratification mask
    '''

    ndvi = combine_dict["ndi_samples"]
    class_mask = np.zeros(ndvi.shape)

    for bin_num in combine_dict['brdf_dict']['bins']:
        start,end =  combine_dict['brdf_dict']['bins'][bin_num]
        class_mask[(ndvi > start) & (ndvi <= end)] = bin_num

    class_mask = class_mask.astype(np.int8)
    combine_dict['ndvi_classes'] = class_mask

def get_topo_var_samples_pre(hy_obj):
    '''Get variables for group_topo correction, run after ndvi_stratify()
    '''
    slope = hy_obj.get_anc('slope')
    cosine_i = hy_obj.cosine_i()
    sample_ind = (hy_obj.ancillary['ndvi_classes'] !=0)

    return slope[sample_ind], cosine_i[sample_ind]

def calc_flex_single_post(combine_data_dict,brdf_dict,load_reflectance_mode):

    combine_data_dict["brdf_dict"] = brdf_dict
    bad_bands = combine_data_dict['bad_bands']
    # Determine bin dimensions and create class mask
    if brdf_dict['bin_type'] == 'dynamic':
        bins = ndvi_bins(combine_data_dict["ndi_samples"],brdf_dict)
        #Update number of bins
        #print(bins)
        combine_data_dict["brdf_dict"]['num_bins']=len(bins)  #hy_obj.brdf['num_bins'] = len(bins)
    else:
        bins = brdf_dict['bins']   

    combine_data_dict['brdf_dict']['bins'] = {k:v for (k,v) in enumerate(bins,start=1)}

    ndvi_stratify_samples(combine_data_dict)

    coeffs = {}
    good_band_count=0
    for band_num,band in enumerate(bad_bands):
        if ~band:
            coeffs[band_num] = {}

            if load_reflectance_mode==0:
                band_samples = combine_data_dict["reflectance_samples"][:,good_band_count]   #ray.get([a.do.remote(get_band_samples,
                                     #{'band_num':band_num}) for a in actors])
            else:
                combine_refl = []
                for h5name in combine_data_dict["reflectance_samples"]:
                    h5_obj = h5py.File(h5name, "r")
                    sub_refl_samples = h5_obj["reflectance_samples"][()][:,good_band_count]
                    combine_refl += [sub_refl_samples]
                    h5_obj.close()
                band_samples = np.concatenate(combine_refl,axis=0)

            band_coeffs= []
            for bin_num in combine_data_dict['brdf_dict']['bins']:

                bin_mask =  (combine_data_dict["ndvi_classes"]== bin_num)

                X = np.concatenate([combine_data_dict["kernels_samples"],np.ones((bin_mask.shape[0],1))],axis=1)[bin_mask]    #kernel_samples[:,:3][bin_mask]
                y = band_samples[bin_mask]
                band_coeffs.append(np.linalg.lstsq(X, y,rcond=-1)[0].flatten().tolist())
            coeffs[band_num]  = band_coeffs
            progbar(np.sum(~bad_bands[:band_num+1]),np.sum(~bad_bands))
            good_band_count+=1

    print('\n')

    combine_data_dict["brdf_dict"]['coeffs'] = coeffs



def calc_flex_single_pre(hy_obj,brdf_dict):
    ''' get samples of a single image for future BRDF coefficients estimation
    '''
    hy_obj.brdf['coeffs'] ={}

    # Determine bin dimensions and create class mask
    if hy_obj.brdf['bin_type'] == 'dynamic':
        bins = ndvi_bins(hy_obj.ndi()[hy_obj.mask['no_data']],brdf_dict)
        #Update number of bins
        hy_obj.brdf['num_bins'] = len(bins)
    else:
        bins = brdf_dict['bins']

    hy_obj.brdf['bins'] = {k:v for (k,v) in enumerate(bins,start=1)}
    ndvi_stratify(hy_obj)
    kernel_samples= get_kernel_samples(hy_obj)

    # Loop each band
    refl_samples_list = []
    used_band = []
    for band_num,band in enumerate(hy_obj.bad_bands):
        if ~band:
            band_samples = hy_obj.do(get_band_samples, {'band_num':band_num})
            refl_samples_list+=[band_samples[:,None]]
            used_band+=[hy_obj.wavelengths[band_num]]

    refl_samples = np.concatenate(refl_samples_list,axis=1)

    slope_samples, cos_i_samples = get_topo_var_samples_pre(hy_obj)  # slope and cosine_i

    return kernel_samples[:,:2], refl_samples, used_band, slope_samples, cos_i_samples


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

def calc_brdf_coeffs_pre(hy_obj,config_dict):

    brdf_dict = config_dict['brdf']

    if brdf_dict['type'] == 'precomputed':
        print("Using precomputed BRDF coefficients")
        load_brdf_precomputed(hy_obj,config_dict['brdf'])

    else:
        # Set BRDF dict

        set_brdf(hy_obj,brdf_dict)
        set_solar_zn_0 = set_solar_zn(hy_obj)

        # Create masks used for calculating coefficients
        hy_obj.gen_mask(mask_create,'calc_brdf',brdf_dict['calc_mask'])

        kernel_samples, reflectance_samples, used_band, slope_samples, cos_i_samples = calc_flex_single_pre(hy_obj,brdf_dict)

    hy_obj.corrections.append('brdf')

    return {
        "set_solar_zn":set_solar_zn_0,
        #"ndvi":hy_obj.ndi(),
        "kernel_samples":kernel_samples,
        "reflectance_samples":reflectance_samples,
        "used_band":used_band,
        "slope_samples":slope_samples,
        "cos_i_samples":cos_i_samples,
    }










