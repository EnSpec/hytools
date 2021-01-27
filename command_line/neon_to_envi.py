''' This script exports NEON AOP HDF imaging spectroscopy data
to and ENVI formated binary file, with option of also exporting
ancillary dataset following formatting used by NASA JPL for AVIRIS
observables.

TODO: Add phase and UTC time to ancillary output

'''
import argparse
import warnings
import os
import numpy as np
import hytools as ht
from hytools.io.envi import WriteENVI
from hytools.misc import progbar

warnings.filterwarnings("ignore")

def main():

    parser = argparse.ArgumentParser(description = "Convert NEON AOP H5 to ENVI format")
    parser.add_argument("-inp", help="Input image pathname",required=True, type = str)
    parser.add_argument("-out", help="Output directory",required=True, type = str)
    parser.add_argument("-anc", help="Ouput ancillary", required=False, action='store_true')

    args = parser.parse_args()

    if not args.out.endswith("/"):
        args.out+="/"

    hy_obj =ht.HyTools()
    hy_obj.read_file(args.inp,'neon')

    output_name = args.out+ os.path.basename(os.path.splitext(args.inp)[0])
    writer = WriteENVI(output_name,hy_obj.get_header())

    print("\nExporting reflectance data.")

    iterator = hy_obj.iterate(by = 'chunk')
    pixels_processed = 0
    while not iterator.complete:
        chunk = iterator.read_next()
        pixels_processed += chunk.shape[0]*chunk.shape[1]
        progbar(pixels_processed, hy_obj.columns*hy_obj.lines)
        writer.write_chunk(chunk,iterator.current_line,iterator.current_column)
        if iterator.complete:
            writer.close()

    if args.anc:
        print("\nExporting ancillary data")
        anc_header = hy_obj.get_header()
        anc_header['bands'] = 10
        anc_header['band_names'] = ['path length', 'to-sensor azimuth', 'to-sensor zenith','to-sun azimuth',
                                      'to-sun zenith','phase', 'slope', 'aspect', 'cosine i','UTC time']
        anc_header['wavelength units'] = np.nan
        anc_header['wavelength'] = np.nan
        anc_header['data type'] = 4

        writer = WriteENVI(output_name + "_ancillary", anc_header)
        writer.write_band(hy_obj.get_anc("path_length"),0)
        writer.write_band(hy_obj.get_anc("sensor_az"),1)
        writer.write_band(hy_obj.get_anc("sensor_zn"),2)
        writer.write_band(hy_obj.get_anc("solar_az"),3)
        writer.write_band(hy_obj.get_anc("solar_zn"),4)
        #writer.write_band(hy_obj.get_anc("phase placeholder"),5)
        writer.write_band(hy_obj.get_anc("slope"),6)
        writer.write_band(hy_obj.get_anc("aspect"),7)
        writer.write_band(hy_obj.cosine_i(),8)
        #writer.write_band('UTC time placeholder',9)
        writer.close()

    print("Export complete.")

if __name__== "__main__":
    main()
