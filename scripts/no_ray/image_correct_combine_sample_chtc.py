
import json
import os
import warnings
import sys
import numpy as np

import h5py

import hytools as ht
from hytools.io.envi import *
from hytools.brdf import calc_flex_single_post   
from hytools.glint import set_glint_parameters
from hytools.masks import mask_create

warnings.filterwarnings("ignore")
np.seterr(divide='ignore', invalid='ignore')

def main():

    config_file = sys.argv[1]
    sample_folder = sys.argv[2]

    with open(config_file, 'r') as outfile:
        config_dict = json.load(outfile)

    images = config_dict["input_files"]
    brdf_dict = config_dict['brdf']

    sample_h5_list = [ f"{sample_folder}/{os.path.splitext(os.path.basename(image))[0]}_prebrdf_sample.h5" for image in images]

    sample_dict = load_sample_h5(sample_h5_list)

    if isinstance(brdf_dict['solar_zn_type'],str):
        if brdf_dict['solar_zn_type'] == 'scene':
            brdf_dict["solar_zn_norm_radians"]=float(sample_dict['mean_solar_zn'])
            print("Scene average solar zenith angle : %s degrees" % round(np.degrees(brdf_dict["solar_zn_norm_radians"]),3))
        elif isinstance(brdf_dict['solar_zn_type'],float):
            brdf_dict["solar_zn_norm_radians"]=brdf_dict['solar_zn_type']
        else:
            print('Unrecognized solar zenith angle normalization')


    calc_flex_single_post(sample_dict,brdf_dict) #,Ht_Obj.bad_bands

    if config_dict['export']['coeffs'] and len(config_dict["corrections"]) > 0:
        print("Exporting correction coefficients.")
        export_coeffs_brdf(sample_dict,config_dict['export'],images)



def load_sample_h5(h5_file_list):
    '''Load information from H5 files, and return a dictionary with all the info needed.
    '''

    combine_refl = []
    combine_kernel = []
    solar_zn_list = []
    ndi_list = []

    bad_bands=None #get from the 1st image

    for i_order, h5name in enumerate(h5_file_list):
        h5_obj = h5py.File(h5name, "r")
        wavelist = h5_obj["wavelengths"][()]
        set_solar_zn = h5_obj["kernels_samples"].attrs['set_solar_zn']
        refl_samples = h5_obj["reflectance_samples"][()]
        kernel_samples = h5_obj["kernels_samples"][()]

        if i_order==0:
            bad_bands = h5_obj["bad_bands"][()] 

        h5_obj.close()

        sample_nir=refl_samples[:,get_wave(850,wavelist)]
        sample_red=refl_samples[:,get_wave(660,wavelist)]
        sample_ndi = (sample_nir-sample_red)/(sample_nir+sample_red)
        ndi_list+=[sample_ndi]

        combine_refl+=[refl_samples]
        solar_zn_list+=[set_solar_zn]
        combine_kernel+=[kernel_samples]

    return {
        "kernels_samples":np.concatenate(combine_kernel,axis=0),
        "reflectance_samples":np.concatenate(combine_refl,axis=0),
        "ndi_samples":np.concatenate(ndi_list),
        "mean_solar_zn": np.array(solar_zn_list).mean(),
        "bad_bands":bad_bands,
    }


def export_coeffs_brdf(data_dict,export_dict,images):
    '''Export correction coefficients to file.
    '''

    for image in images:
        coeff_file = export_dict['output_dir']
        coeff_file += os.path.splitext(os.path.basename(image))[0]
        coeff_file += "_%s_coeffs_%s_chtc.json" % ("brdf",export_dict["suffix"])

        with open(coeff_file, 'w') as outfile:
            corr_dict = data_dict['brdf_dict']  #hy_obj.brdf
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

if __name__== "__main__":
    main()
