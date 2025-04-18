
import json
import os
import warnings
import sys
import ray
import numpy as np

import time

import hytools as ht
from hytools.io.envi import *
from hytools.topo import calc_topo_coeffs
from hytools.brdf import calc_brdf_coeffs
from hytools.masks import mask_create


warnings.filterwarnings("ignore")
np.seterr(divide='ignore', invalid='ignore')

def main():

    time_start = time.perf_counter()

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
    elif config_dict['file_type'] == 'emit' or config_dict['file_type'] == 'ncav':
        anc_files = config_dict["anc_files"]
        if bool(config_dict["glt_files"]):
            glt_files = config_dict["glt_files"]
            _ = ray.get([a.read_file.remote(image,config_dict['file_type'],
                                        anc_path=anc_files[image],glt_path=glt_files[image]) for a,image in zip(actors,images)])
        else:    
            _ = ray.get([a.read_file.remote(image,config_dict['file_type'],
                                        anc_path=anc_files[image]) for a,image in zip(actors,images)])

    _ = ray.get([a.create_bad_bands.remote(config_dict['bad_bands']) for a in actors])

    for correction in config_dict["corrections"]:
        if correction =='topo':
            time_topo_start = time.perf_counter()
            #set_trace()
            calc_topo_coeffs(actors,config_dict['topo'])
            time_topo_end = time.perf_counter() 
            print("TOPO Time: {} sec.".format(time_topo_end - time_topo_start))
        elif correction == 'brdf':
            time_brdf_start = time.perf_counter()
            calc_brdf_coeffs(actors,config_dict)
            time_brdf_end = time.perf_counter() 
            print("BRDF Time: {} sec.".format(time_brdf_end - time_brdf_start))


    if config_dict['export']['coeffs'] and len(config_dict["corrections"]) > 0:
        print("Exporting correction coefficients.")
        _ = ray.get([a.do.remote(export_coeffs,config_dict['export']) for a in actors])


    time_export_start = time.perf_counter()
    if config_dict['export']['image']:
        print("Exporting corrected images.")
        _ = ray.get([a.do.remote(apply_corrections,config_dict) for a in actors])
    time_export_end = time.perf_counter()
    print("{} Export Time: {} sec.".format(images[0],time_export_end - time_export_start))

    ray.shutdown()

    time_end = time.perf_counter()
    print("Total Time: {} sec.".format(time_end - time_start))

def export_coeffs(hy_obj,export_dict):
    '''Export correction coefficients to file.
    '''
    for correction in hy_obj.corrections:
        if correction=='unsmooth':
            continue

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
    
    use_glt_output_bool=False
    if 'use_glt' in config_dict['export']:
        use_glt_output_bool = config_dict['export']['use_glt']
        if use_glt_output_bool==True:
            header_dict = hy_obj.get_header(warp_glt=True)
        else:
            header_dict = hy_obj.get_header()    
    else:
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
        header_dict['fwhm'] = hy_obj.fwhm

        writer = WriteENVI(output_name,header_dict)
        iterator = hy_obj.iterate(by = 'line', corrections = hy_obj.corrections,
                                  resample = config_dict['resample'])
        
        if use_glt_output_bool==False:
            while not iterator.complete:
                line = iterator.read_next()
                writer.write_line(line,iterator.current_line)
        else:
            while not iterator.complete:
                line = iterator.read_next()
                line_ind=np.where(hy_obj.glt_y==iterator.current_line + 1)
                redundant_line = line[hy_obj.glt_x[line_ind]-1]
                writer.write_line_glt(redundant_line[:,:2],line_ind[0],line_ind[1])
        
        writer.close()

    #Export subset of wavelengths
    else:
        waves = config_dict['export']['subset_waves']
        bands = [hy_obj.wave_to_band(x) for x in waves]
        waves = [round(hy_obj.wavelengths[x],2) for x in bands]
        header_dict['bands'] = len(bands)
        header_dict['wavelength'] = waves
        header_dict['fwhm'] = [hy_obj.fwhm[x] for x in bands]

        writer = WriteENVI(output_name,header_dict)
        if use_glt_output_bool==False:
            for b,band_num in enumerate(bands):
                band = hy_obj.get_band(band_num,
                                    corrections = hy_obj.corrections)
                writer.write_band(band, b)
        else:
            for b,band_num in enumerate(bands):
                band = hy_obj.get_band(band_num,
                                    corrections = hy_obj.corrections)
                writer.write_band_glt(band,b, (hy_obj.glt_y[hy_obj.fill_mask]-1,hy_obj.glt_x[hy_obj.fill_mask]-1),hy_obj.fill_mask)
        writer.close()

    #Export masks
    # does not work for precomputed json coeffs
    if (config_dict['export']['masks']) and (len(config_dict["corrections"]) > 0):
        masks = []
        mask_names = []

        for correction in config_dict["corrections"]:
            if correction=='glint':
                continue

            for mask_type in config_dict[correction]['apply_mask']:
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

        if use_glt_output_bool==False:
            for band_num,mask in enumerate(masks):
                mask =mask.astype(int)
                mask[~hy_obj.mask['no_data']] = 255
                writer.write_band(mask,band_num)
        else:
            for band_num,mask in enumerate(masks):
                mask =mask.astype(int)
                mask[~hy_obj.mask['no_data']] = 255
                writer.write_band_glt(mask,band_num, (hy_obj.glt_y[hy_obj.fill_mask]-1,hy_obj.glt_x[hy_obj.fill_mask]-1),hy_obj.fill_mask)
        del masks

if __name__== "__main__":
    main()
