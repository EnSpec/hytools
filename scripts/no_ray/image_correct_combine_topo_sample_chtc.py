

import json
import os
import warnings
import sys
import numpy as np
from scipy.optimize import nnls

import h5py

import hytools as ht
from hytools.io.envi import *
from hytools.masks import mask_create

from hytools.topo.c import calc_c

warnings.filterwarnings("ignore")
np.seterr(divide='ignore', invalid='ignore')

def main():

    config_file = sys.argv[1]
    sample_folder = sys.argv[2]
    topo_subgroup_id = sys.argv[3]

    with open(config_file, 'r') as outfile:
        config_dict = json.load(outfile)

    images = []

    if config_dict["topo"]["subgrouped"]:
        topo_dict = config_dict['topo']
        subgroup = topo_dict["subgroup"]

        sample_h5_list=[]
        for each_img_name in subgroup.keys():
            each_h5_name = f"{sample_folder}/{os.path.splitext(os.path.basename(each_img_name))[0]}_pretopo_sample.h5"
            if subgroup[each_img_name]==topo_subgroup_id:
                sample_h5_list+=[each_h5_name]
                images+=[each_img_name]

        if len(sample_h5_list)==0:
            print(f"Cannot find subgroup '{topo_subgroup_id}', exit.")
            return

        sample_dict = load_sample_h5(sample_h5_list)

        calc_topo_single_post(sample_dict,topo_dict) #,update topo coeffs

        apply_topo_scsc(sample_dict) # update reflectance

        if config_dict['export']['coeffs'] and len(config_dict["corrections"]) > 0:
            print("Exporting correction coefficients.")
            export_coeffs_topo(sample_dict,config_dict['export'],images,sample_h5_list)

        for image_order, imagename in enumerate(images):
            export_h5(imagename,config_dict['export'],sample_dict[sample_h5_list[image_order]])

    else:
        print("No subgroup is defined, exit.")


def load_sample_h5(h5_file_list):

    h5_all_dict = {}

    bad_bands=None #get from the 1st image

    for i_order, h5name in enumerate(h5_file_list):
        h5_obj = h5py.File(h5name, "r")
        wavelist = h5_obj["wavelengths"][()]
        full_image_wavelist = h5_obj["image_wavelengths"][()]
        set_solar_zn = h5_obj["kernels_samples"].attrs['set_solar_zn']
        refl_samples = h5_obj["reflectance_samples"][()]
        kernel_samples = h5_obj["kernels_samples"][()]
        slope_samples = h5_obj["slope_samples"][()]
        cosine_i_samples = h5_obj["cosine_i_samples"][()]
        bad_bands = h5_obj["bad_bands"][()]

        h5_obj.close()

        sample_nir=refl_samples[:,get_wave(850,wavelist)]
        sample_red=refl_samples[:,get_wave(660,wavelist)]
        sample_ndi = (sample_nir-sample_red)/(sample_nir+sample_red)

        h5_all_dict[h5name] =  {
          "kernels_samples":kernel_samples,
          "reflectance_samples":refl_samples,
          "ndi_samples":sample_ndi,
          "bad_bands":bad_bands,
          "full_image_wavelist":full_image_wavelist,
          "set_solar_zn":set_solar_zn,
          "wavelist":wavelist,
          "slope_samples":slope_samples,
          "cosine_i_samples":cosine_i_samples,
          "topo_dict":None,
        }

    return h5_all_dict


def export_coeffs_topo(data_dict,export_dict,images,h5_list):
    '''Export correction coefficients to file.
    '''

    for img_order, image in enumerate(images):
        coeff_file = export_dict['output_dir']
        coeff_file += os.path.splitext(os.path.basename(image))[0]
        coeff_file += "_%s_coeffs_%s_chtc.json" % ("topo",export_dict["suffix"])

        with open(coeff_file, 'w') as outfile:
            corr_dict = data_dict[h5_list[img_order]]['topo_dict']
            json.dump(corr_dict,outfile)

def get_wave(wave,wavelengths):
    """Return the band image corresponding to the input wavelength.
    If not an exact match the closest wavelength will be returned.

    Args:
        wave (float): Wavelength in image units.
        wavelengths (list): Wavelength list


    Returns:
        band index 0-based.

    """

    if (wave  > wavelengths.max()) | (wave  < wavelengths.min()):
        print("Input wavelength outside wavelength range!")
        band_ind = None
    else:
        band_ind = np.argmin(np.abs(wavelengths - wave))
    return band_ind

def export_h5(imagename,export_dict,obj_dict):

    out_filename = f"{export_dict['output_dir']}{os.path.splitext(os.path.basename(imagename))[0]}_prebrdf_sample.h5"

    h5_obj = h5py.File(out_filename, "w")
    h5_obj.attrs['Image Name']=f"{os.path.splitext(os.path.basename(imagename))[0]}"

    dset1 = h5_obj.create_dataset("kernels_samples", data=obj_dict["kernels_samples"])
    dset2 = h5_obj.create_dataset("reflectance_samples", data=obj_dict["reflectance_samples"])
    dset3 = h5_obj.create_dataset("wavelengths", data=obj_dict["wavelist"])
    dset1.attrs['set_solar_zn']=obj_dict['set_solar_zn']
    dset1.attrs['kernels_names']='["Volume","Geometry"]'
    dset1.attrs['Solar Zenith Unit']="Radians"

    dset4 = h5_obj.create_dataset("image_wavelengths", data=obj_dict["full_image_wavelist"])
    dset5 = h5_obj.create_dataset("bad_bands", data=obj_dict["bad_bands"])

    h5_obj.close()
    print(f"{out_filename} saved.")

def calc_topo_single_post(sample_dict,topo_dict):
    combine_refl = []
    combine_cos_i = []
    combine_slope = []
    ndi_list = []

    for h5_name in sample_dict.keys():
        sub_dict = sample_dict[h5_name]
        ndi_list+=[sub_dict["ndi_samples"]]
        combine_refl+=[sub_dict["reflectance_samples"]]
        combine_cos_i+=[sub_dict["cosine_i_samples"]]
        combine_slope+=[sub_dict["slope_samples"]]
        bad_bands = sub_dict["bad_bands"]

    combine_refl=np.concatenate(combine_refl,axis=0)
    combine_cos_i=np.concatenate(combine_cos_i,axis=0)
    combine_slope=np.concatenate(combine_slope,axis=0)
    ndi_list=np.concatenate(ndi_list,axis=0)

    mask = np.ones(ndi_list.shape).astype(bool)

    mask &= (ndi_list >= float(topo_dict['calc_mask'][0][1]['min'])) & (ndi_list <= float(topo_dict['calc_mask'][0][1]['max']))
    mask &= (combine_slope >= float(topo_dict['calc_mask'][1][1]['min'])) & (combine_slope <= float(topo_dict['calc_mask'][1][1]['max']))
    mask &= (combine_cos_i >= float(topo_dict['calc_mask'][2][1]['min'])) & (combine_cos_i <= float(topo_dict['calc_mask'][2][1]['max']))

    feasible_sample_count=np.count_nonzero(mask)

    if feasible_sample_count>10:
        used_reflectance_samples = combine_refl[mask==1,:]
        used_cos_i = combine_cos_i[mask==1]

        topo_dict['coeffs'] = {}

        band_cursor=0
        for band_num,band in enumerate(bad_bands):
            if ~band:
                topo_dict['coeffs'][band_num] = calc_c(used_reflectance_samples[:,band_cursor],used_cos_i,
                                                       fit_type=topo_dict['c_fit_type'])
                band_cursor+=1
    else:
        topo_dict['coeffs'] = {}
        band_cursor=0
        for band_num,band in enumerate(bad_bands):
            if ~band:
                topo_dict['coeffs'][band_num] = 100000.0
                band_cursor+=1

    for h5_name in sample_dict.keys():
        sub_dict = sample_dict[h5_name]
        sub_dict["topo_dict"] = topo_dict

def apply_topo_scsc(sample_dict):
    for h5_name in sample_dict.keys():        
        data_dict = sample_dict[h5_name]
        topo_dict = data_dict["topo_dict"]
        slope_samples = data_dict['slope_samples']
        c1 = np.cos(slope_samples) * np.cos(data_dict['set_solar_zn'])
        cosine_i = data_dict['cosine_i_samples']
        ndi_list = data_dict['ndi_samples']
        refl_samples = data_dict['reflectance_samples']

        C_arr = np.array(list(topo_dict['coeffs'].values()))

        mask = np.ones(ndi_list.shape).astype(bool)

        # mask order in the config matters here
        mask &= (ndi_list >= float(topo_dict['apply_mask'][0][1]['min'])) & (ndi_list <= float(topo_dict['apply_mask'][0][1]['max']))
        mask &= (slope_samples >= float(topo_dict['apply_mask'][1][1]['min'])) & (slope_samples <= float(topo_dict['apply_mask'][1][1]['max']))
        mask &= (cosine_i >= float(topo_dict['apply_mask'][2][1]['min'])) & (cosine_i <= float(topo_dict['apply_mask'][2][1]['max']))

        for band_order in range(refl_samples.shape[1]):
            band = np.copy(refl_samples[:,band_order])
            correction_factor = (c1 + C_arr[band_order])/(cosine_i + C_arr[band_order])
            band[mask] = band[mask]*correction_factor[mask]
            refl_samples[:,band_order]=band


if __name__== "__main__":
    main()
