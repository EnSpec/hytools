'''tanager2envi.py
'''
import argparse
import os
import sys
import ray
import numpy as np

import hytools as ht
from hytools.io.envi import WriteENVI

def main():
    '''This command line tool exports Tanager Reflectance HDF imaging spectroscopy data
    to an ENVI formatted binary file, with the option of also exporting
    ancillary data following formatting used by NASA JPL for AVIRIS
    observables. The script utilizes ray to export images in parallel.

    '''
    parser = argparse.ArgumentParser(description = "Convert Tanager1 Ortho HDF to ENVI format")
    parser.add_argument('images',help="Input image pathnames", nargs='*')
    parser.add_argument('output_dir',help="Output directory", type = str)
    parser.add_argument("-rfl", help="Output reflectance image", required=False, action='store_true')
    parser.add_argument("-anc", help="Output ancillary", required=False, action='store_true')
    parser.add_argument("-mask", help="Output existing mask", required=False, action='store_true')

    args = parser.parse_args()

    if not args.output_dir.endswith("/"):
        args.output_dir+="/"

    if ray.is_initialized():
        ray.shutdown()
    ray.init(num_cpus = len(args.images))

    hytool = ray.remote(ht.HyTools)
    actors = [hytool.remote() for image in args.images]
    _ = ray.get([a.read_file.remote(image,'tanager') for a,image in zip(actors,args.images)])

    def tanager_to_envi(hy_obj):
        basemame = os.path.basename(os.path.splitext(hy_obj.file_name)[0])
        print("Exporting %s " % basemame)
        output_name = args.output_dir+ basemame
        writer = WriteENVI(output_name,hy_obj.get_header())

        iterator = hy_obj.iterate(by = 'chunk',
                                  chunk_size = (int(np.ceil(hy_obj.lines/16)),int(np.ceil(hy_obj.columns/16))))
        pixels_processed = 0
        while not iterator.complete:
            chunk = iterator.read_next()
            pixels_processed += chunk.shape[0]*chunk.shape[1]
            writer.write_chunk(chunk,iterator.current_line,iterator.current_column)
            if iterator.complete:
                writer.close()

    def export_anc(hy_obj):
        anc_header = hy_obj.get_header()
        anc_header['bands'] = 5
        anc_header['band_names'] = [  'to-sensor azimuth',
                                    'to-sensor zenith','to-sun azimuth',
                                      'to-sun zenith','UTC time']
        anc_header['wavelength units'] = np.nan
        anc_header['wavelength'] = np.nan
        anc_header['data type'] = 4

        output_name = args.output_dir+ os.path.basename(os.path.splitext(hy_obj.file_name)[0])
        writer = WriteENVI(output_name + "_ancillary", anc_header)
        #writer.write_band(hy_obj.get_anc("path_length",radians = False),1)
        writer.write_band(hy_obj.get_anc("sensor_az",radians = False),0)
        writer.write_band(hy_obj.get_anc("sensor_zn",radians = False),1)
        writer.write_band(hy_obj.get_anc("solar_az",radians = False),2)
        writer.write_band(hy_obj.get_anc("solar_zn",radians = False),3)
        #writer.write_band(hy_obj.get_anc("phase placeholder"),5)
        #writer.write_band(hy_obj.get_anc("slope",radians = False),6)
        #writer.write_band(hy_obj.get_anc("aspect",radians = False),7)
        #writer.write_band(hy_obj.cosine_i(),8)
        writer.write_band(hy_obj.get_anc("UTC time",radians = False),4)
        writer.close()

    def export_msk(hy_obj):
        anc_header = hy_obj.get_header()
        anc_header['bands'] = 3
        anc_header['band_names'] = ['beta_cirrus_mask',
                                    'beta_cloud_mask',
                                    'nodata_pixels']

        anc_header['data type'] = 2

        output_name = args.output_dir+ os.path.basename(os.path.splitext(hy_obj.file_name)[0])
        writer = WriteENVI(output_name + "_mask", anc_header)
        writer.write_band(hy_obj.get_anc("beta_cirrus_mask",radians = False),0)
        writer.write_band(hy_obj.get_anc("beta_cloud_mask",radians = False),1)
        writer.write_band(hy_obj.get_anc("nodata_pixels",radians = False),2)
        writer.close()

    if args.rfl:
        print("\nExporting reflectance data")
        _ = ray.get([a.do.remote(tanager_to_envi) for a in actors])

    if args.anc:
        print("\nExporting ancillary data")
        _ = ray.get([a.do.remote(export_anc) for a in actors])

    if args.mask:
        print("\nExporting mask data")
        _ = ray.get([a.do.remote(export_msk) for a in actors])

    print("Export complete.")


if __name__== "__main__":
    main()
