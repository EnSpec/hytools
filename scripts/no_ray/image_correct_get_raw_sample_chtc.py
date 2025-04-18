

import json
import os
import warnings
import sys
import numpy as np

import h5py

import hytools as ht
from hytools.io.envi import *
from hytools.topo import calc_topo_coeffs_single
from hytools.brdf import calc_brdf_coeffs_pre
from hytools.glint import set_glint_parameters
from hytools.masks import mask_create

warnings.filterwarnings("ignore")
np.seterr(divide='ignore', invalid='ignore')

def main():

    config_file = sys.argv[1]
    flightline_index = int(sys.argv[2])
    group_topo_bool = bool(int(sys.argv[3]))
    print("For Group Topo (no correction applied to samples)",group_topo_bool)

    with open(config_file, 'r') as outfile:
        config_dict = json.load(outfile)

    image = config_dict["input_files"][flightline_index]

    Ht_Obj = ht.HyTools()

    if config_dict['file_type'] == 'envi':
        anc_files = config_dict["anc_files"]
        Ht_Obj.read_file(image,config_dict['file_type'],anc_files[image])

    elif config_dict['file_type'] == 'neon':
        Ht_Obj.read_file(image,config_dict['file_type'])

    elif config_dict['file_type'] == 'ncav' or config_dict['file_type'] == 'emit':
        anc_files = config_dict["anc_files"]
        Ht_Obj.read_file(image,config_dict['file_type'],anc_files[image])

    Ht_Obj.create_bad_bands(config_dict['bad_bands'])

    non_water_count = np.count_nonzero(Ht_Obj.ndi()[Ht_Obj.mask['no_data']]>0.001)
    if non_water_count<50:
        print("Not enough ground pixels... exit")
        return

    for correction in config_dict["corrections"]:
        if correction =='topo':
            if group_topo_bool is False: # get single line topo coeffs 
                calc_topo_coeffs_single(Ht_Obj,config_dict['topo'])
        elif correction == 'brdf':
            data_dict=calc_brdf_coeffs_pre(Ht_Obj,config_dict)
            print(f"{data_dict['kernel_samples'].shape[0]} pixels extracted.")
            export_h5(Ht_Obj,config_dict['export'],data_dict,group_topo_bool)


    if (config_dict['export']['coeffs'] and len(config_dict["corrections"]) > 0) and group_topo_bool is False:
        print("Exporting correction coefficients.")
        export_coeffs_topo(Ht_Obj,config_dict['export'])


def export_h5(hy_obj,export_dict,obj_dict,topo_bool):

    if topo_bool==False:
        out_filename = f"{export_dict['output_dir']}{os.path.splitext(os.path.basename(hy_obj.file_name))[0]}_prebrdf_sample.h5"
    else:
        out_filename = f"{export_dict['output_dir']}{os.path.splitext(os.path.basename(hy_obj.file_name))[0]}_pretopo_sample.h5"

    h5_obj = h5py.File(out_filename, "w")

    dset1 = h5_obj.create_dataset("kernels_samples", data=obj_dict["kernel_samples"]) #chunks=(50, 50), compression="gzip"
    dset2 = h5_obj.create_dataset("reflectance_samples", data=obj_dict["reflectance_samples"])
    dset3 = h5_obj.create_dataset("wavelengths", data=obj_dict["used_band"])
    dset1.attrs['set_solar_zn']=obj_dict['set_solar_zn']
    dset1.attrs['kernels_names']='["Volume","Geometry"]'
    dset1.attrs['Solar Zenith Unit']="Radians"

    dset4 = h5_obj.create_dataset("image_wavelengths", data=np.array(hy_obj.wavelengths))
    dset5 = h5_obj.create_dataset("bad_bands", data=np.array(hy_obj.bad_bands))

    dset6 = h5_obj.create_dataset("slope_samples", data=obj_dict["slope_samples"])
    dset6.attrs['Slope Unit']="Radians"
    dset7 = h5_obj.create_dataset("cosine_i_samples", data=obj_dict["cos_i_samples"])    

    h5_obj.close()
    print(f"{out_filename} saved.")


def export_coeffs_topo(hy_obj,export_dict):
    '''Export correction coefficients to file.
    '''
    for correction in hy_obj.corrections:
        if not correction == 'topo':
            continue

        coeff_file = export_dict['output_dir']
        coeff_file += os.path.splitext(os.path.basename(hy_obj.file_name))[0]
        coeff_file += "_%s_coeffs_%s_chtc.json" % (correction,export_dict["suffix"])

        with open(coeff_file, 'w') as outfile:
            if correction == 'topo':
                corr_dict = hy_obj.topo
            elif correction == 'glint':
                continue
            else:
                corr_dict = hy_obj.brdf
            json.dump(corr_dict,outfile)

if __name__== "__main__":
    main()