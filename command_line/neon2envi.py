'''neon2envi.py

TODO: Add phase and UTC time to ancillary output
TODO: Implement progress bar like from this example:
      https://docs.ray.io/en/master/auto_examples/progress_bar.html

'''
import argparse
import os
import ray
import numpy as np
import hytools as ht
from hytools.io.envi import WriteENVI

def main():
    '''This command line tool exports NEON AOP HDF imaging spectroscopy data
    to an ENVI formated binary file, with the option of also exporting
    ancillary data following formatting used by NASA JPL for AVIRIS
    observables. The script utilizes ray to export images in parralel.

    '''
    parser = argparse.ArgumentParser(description = "Convert NEON AOP H5 to ENVI format")
    parser.add_argument('images',help="Input image pathnames", nargs='*')
    parser.add_argument('output_dir',help="Output directory", type = str)
    parser.add_argument("-anc", help="Output ancillary", required=False, action='store_true')

    args = parser.parse_args()

    if not args.output_dir.endswith("/"):
        args.output_dir+="/"

    if ray.is_initialized():
        ray.shutdown()
    ray.init(num_cpus = len(args.images))

    hytool = ray.remote(ht.HyTools)
    actors = [hytool.remote() for image in args.images]
    _ = ray.get([a.read_file.remote(image,'neon') for a,image in zip(actors,args.images)])

    def neon_to_envi(hy_obj):
        basemame = os.path.basename(os.path.splitext(hy_obj.file_name)[0])
        print("Exporting %s " % basemame)
        output_name = args.output_dir+ basemame
        writer = WriteENVI(output_name,hy_obj.get_header())
        iterator = hy_obj.iterate(by = 'chunk')
        pixels_processed = 0
        while not iterator.complete:
            chunk = iterator.read_next()
            pixels_processed += chunk.shape[0]*chunk.shape[1]
            writer.write_chunk(chunk,iterator.current_line,iterator.current_column)
            if iterator.complete:
                writer.close()

    def export_anc(hy_obj):
        anc_header = hy_obj.get_header()
        anc_header['bands'] = 10
        anc_header['band_names'] = ['path length', 'to-sensor azimuth',
                                    'to-sensor zenith','to-sun azimuth',
                                      'to-sun zenith','phase', 'slope',
                                      'aspect', 'cosine i','UTC time']
        anc_header['wavelength units'] = np.nan
        anc_header['wavelength'] = np.nan
        anc_header['data type'] = 4

        output_name = args.output_dir+ os.path.basename(os.path.splitext(hy_obj.file_name)[0])
        writer = WriteENVI(output_name + "_ancillary", anc_header)
        writer.write_band(hy_obj.get_anc("path_length"),0)
        writer.write_band(hy_obj.get_anc("sensor_az",radians = False),1)
        writer.write_band(hy_obj.get_anc("sensor_zn",radians = False),2)
        writer.write_band(hy_obj.get_anc("solar_az",radians = False),3)
        writer.write_band(hy_obj.get_anc("solar_zn",radians = False),4)
        #writer.write_band(hy_obj.get_anc("phase placeholder"),5)
        writer.write_band(hy_obj.get_anc("slope",radians = False),6)
        writer.write_band(hy_obj.get_anc("aspect",radians = False),7)
        writer.write_band(hy_obj.cosine_i(),8)
        #writer.write_band('UTC time placeholder',9)
        writer.close()

    _ = ray.get([a.do.remote(neon_to_envi) for a in actors])

    if args.anc:
        print("\nExporting ancillary data")
        _ = ray.get([a.do.remote(export_anc) for a in actors])

    print("Export complete.")


if __name__== "__main__":
    main()
