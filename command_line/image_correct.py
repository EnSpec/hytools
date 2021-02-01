import glob
import ray
import numpy as np
import json
import os
import warnings
import sys
import hytools as ht
from hytools.io.envi import *
from hytools.topo.topo import topo_coeffs
from hytools.brdf.brdf import brdf_coeffs

warnings.filterwarnings("ignore")
np.seterr(divide='ignore', invalid='ignore')

def main():

    config_file = sys.argv[1]

    with open(config_file, 'r') as outfile:
        config_dict = json.load(outfile)

    images= config_dict["input_files"]

    if ray.is_initialized():
        ray.shutdown()
    print("Using %s CPUs." % config_dict['num_cpus'])
    ray.init(num_cpus = config_dict['num_cpus'])

    HyTools = ray.remote(ht.HyTools)
    actors = [HyTools.remote() for image in images]

    if config_dict['file_type'] == 'envi':
        anc_files = config_dict["anc_files"]
        _ = ray.get([a.read_file.remote(image,config_dict['file_type'],
                                        anc_files[image]) for a,image in zip(actors,images)])
    elif config_dict['file_type'] == 'neon':
        _ = ray.get([a.read_file.remote(image,config_dict['file_type']) for a,image in zip(actors,images)])

    #Here is where the outlier detection should probably happen.

    _ = ray.get([a.create_bad_bands.remote(config_dict['bad_bands']) for a in actors])

    for correction in config_dict["corrections"]:
        if correction =='topo':
            topo_coeffs(actors,config_dict['topo'])
        elif correction == 'brdf':
            brdf_coeffs(actors,config_dict)

    if config_dict['export']['coeffs'] and len(config_dict["corrections"]) > 0:
        print("Exporting correction coefficients.")
        _ = ray.get([a.do.remote(export_coeffs,config_dict['export']) for a in actors])

    if config_dict['export']['image']:
        print("Exporting corrected images.")
        _ = ray.get([a.do.remote(apply_corrections,config_dict) for a in actors])

    ray.shutdown()

def export_coeffs(hy_obj,export_dict):
    '''Export correction coefficients to file.
    '''
    for correction in hy_obj.corrections:
        coeff_file = export_dict['output_dir']
        coeff_file += os.path.splitext(os.path.basename(hy_obj.file_name))[0]
        coeff_file += "_%s_coeffs_%s.json" % (correction,export_dict["suffix"])

        with open(coeff_file, 'w') as outfile:
            if correction == 'topo':
                corr_dict = hy_obj.topo
            else:
                corr_dict = hy_obj.brdf
            json.dump(corr_dict,outfile)

def apply_corrections(hy_obj,config_dict):
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
        iterator = hy_obj.iterate(by = 'line', corrections = hy_obj.corrections,
                                  resample = config_dict['resample'])
        while not iterator.complete:
            line = iterator.read_next()
            writer.write_line(line,iterator.current_line)
            progbar(iterator.current_line,hy_obj.lines)

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
                                   corrections = hy_obj.corrections)
            writer.write_band(band, b)
        writer.close()

if __name__== "__main__":
    main()
