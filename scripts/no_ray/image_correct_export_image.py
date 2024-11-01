
import json
import os
import warnings
import sys
import numpy as np

import hytools as ht
from hytools.io.envi import *
from hytools.topo import load_topo_precomputed
from hytools.brdf import load_brdf_precomputed
from hytools.glint import set_glint_parameters_single
from hytools.masks import mask_create

warnings.filterwarnings("ignore")
np.seterr(divide='ignore', invalid='ignore')

def main():

    config_file = sys.argv[1]
    image_order = int(sys.argv[2])

    with open(config_file, 'r') as outfile:
        config_dict = json.load(outfile)

    image = config_dict["input_files"][image_order]

    actor = ht.HyTools()

    if config_dict['file_type'] == 'envi':
        anc_files = config_dict["anc_files"]
        actor.read_file(image,config_dict['file_type'],anc_files[image])

    elif config_dict['file_type'] == 'neon':
        actor.read_file(image,config_dict['file_type'])

    actor.create_bad_bands(config_dict['bad_bands'])

    for correction in config_dict["corrections"]:
        if correction =='topo':
            if config_dict['topo']['type'] == 'precomputed':
                print("Using precomputed topographic coefficients.")
                load_topo_precomputed(actor,config_dict['topo'])
                actor.corrections.append('topo')
            else:
                print('Only precomputed topographic coefficients are accepted. Quit.')
                return
        elif correction == 'brdf':
            if config_dict['brdf']['type'] == 'precomputed':
                print("Using precomputed BRDF coefficients.")
                load_brdf_precomputed(actor,config_dict['brdf'])
                actor.corrections.append('brdf')
            else:
                print('Only precomputed BRDF coefficients are accepted. Quit.')
                return

        elif correction == 'glint':
            set_glint_parameters_single(actor, config_dict)

    if config_dict['export']['image']:
        print("Exporting corrected image.")
        apply_corrections_single(actor,config_dict)



def apply_corrections_single(hy_obj,config_dict):
    '''Apply correction to image and export
        to file.
    '''

    header_dict = hy_obj.get_header()
    header_dict['data ignore value'] = hy_obj.no_data
    header_dict['data type'] = 4

    output_name = config_dict['export']['output_dir']
    output_name += os.path.splitext(os.path.basename(hy_obj.file_name))[0]
    output_name +=  "_%s" % config_dict['export']["suffix"]
    
    #Export all wavelengths
    if len(config_dict['export']['subset_waves']) == 0:

        if config_dict["resample"] == True:
            hy_obj.resampler = config_dict['resampler']
            waves= hy_obj.resampler['out_waves']
        else:
            waves = hy_obj.wavelengths

        header_dict['bands'] = len(waves)
        header_dict['wavelength'] = waves

        writer = WriteENVI(output_name,header_dict)
        iterator = hy_obj.iterate(by='line', corrections=hy_obj.corrections,
                                  resample=config_dict['resample'])
        while not iterator.complete:
            line = iterator.read_next()
            writer.write_line(line,iterator.current_line)
        writer.close()

    #Export subset of wavelengths
    else:
        waves = config_dict['export']['subset_waves']
        bands = [hy_obj.wave_to_band(x) for x in waves]
        waves = [round(hy_obj.wavelengths[x],2) for x in bands]
        header_dict['bands'] = len(bands)
        header_dict['wavelength'] = waves

        writer = WriteENVI(output_name,header_dict)
        for b,band_num in enumerate(bands):
            band = hy_obj.get_band(band_num,
                                   corrections=hy_obj.corrections)
            writer.write_band(band, b)
        writer.close()
    
    #Export masks
    if (config_dict['export']['masks']) and (len(config_dict["corrections"]) > 0):
        masks = []
        mask_names = []
        
        for correction in config_dict["corrections"]:
            for mask_type in getattr(hy_obj,correction)['apply_mask']:
                mask_names.append(correction + '_' + mask_type[0])
                masks.append(mask_create(hy_obj, [mask_type]))

        header_dict['data type'] = 1
        header_dict['bands'] = len(masks)
        header_dict['band names'] = mask_names
        header_dict['samples'] = hy_obj.columns
        header_dict['lines'] = hy_obj.lines
        header_dict['wavelength'] = []
        header_dict['fwhm'] = []
        header_dict['wavelength units'] = ''
        header_dict['data ignore value'] = 255


        output_name = config_dict['export']['output_dir']
        output_name += os.path.splitext(os.path.basename(hy_obj.file_name))[0]
        output_name +=  "_%s_mask" % config_dict['export']["suffix"]

        writer = WriteENVI(output_name,header_dict)

        for band_num,mask in enumerate(masks):
            mask = mask.astype(int)
            mask[~hy_obj.mask['no_data']] = 255
            writer.write_band(mask,band_num)

        del masks


if __name__== "__main__":
    main()