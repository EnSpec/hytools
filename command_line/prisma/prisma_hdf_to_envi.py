'''
TODO: Docstrings
'''
import os
import argparse
import warnings
import h5py
from hytools.io.envi import WriteENVI,envi_header_dict

warnings.filterwarnings("ignore")

def progbar(curr, total, full_progbar=100):
    '''Progress bar
    '''
    frac = curr/total
    filled_progbar = round(frac*full_progbar)
    print('\r', '#'*filled_progbar + '-'*(full_progbar-filled_progbar), '[{:>7.2%}]'.format(frac), end='')

def main():
    '''
    :return: DESCRIPTION
    :rtype: TYPE
    '''
    parser = argparse.ArgumentParser(description = "Convert PRISMA H5 to merged VNIR-SWIR ENVI format")
    parser.add_argument("--img", help="Input image pathname",required=True, type = str)
    parser.add_argument("--out", help="Output directory",required=True, type = str)

    args = parser.parse_args()

    if not args.out.endswith("/"):
        args.out+="/"

    # Load HDF
    hdf_obj = h5py.File(args.img,'r')
    base_key = [key for key in hdf_obj['HDFEOS']["SWATHS"].keys() if 'HCO' in key][0]
    vnir_data = hdf_obj['HDFEOS']["SWATHS"][base_key]['Data Fields']['VNIR_Cube']
    swir_data =  hdf_obj['HDFEOS']["SWATHS"][base_key]['Data Fields']['SWIR_Cube']
    upper_x = hdf_obj.attrs.get('Product_ULcorner_easting')
    upper_y = hdf_obj.attrs.get('Product_ULcorner_northing')

    resolution =30
    projection =hdf_obj.attrs.get('Projection_Name').decode('UTF-8')
    zone =hdf_obj.attrs.get('Projection_Id').decode('UTF-8')
    direction = 'N'
    if hdf_obj.attrs.get('Product_center_lat') < 0:
        direction = 'S'#May not be accurate in scenes that straddle the equator....
    map_info = [projection, 1, 1,upper_x ,
                       upper_y,resolution,resolution,
                       zone, direction, 'WGS-84','units=Meters']

    output_name = args.out+ os.path.basename(os.path.splitext(args.img)[0])

    vnir_waves =  hdf_obj.attrs.get('List_Cw_Vnir')[6:].tolist()
    swir_waves =  hdf_obj.attrs.get('List_Cw_Swir')[:-3].tolist()
    vnir_waves.reverse()
    swir_waves.reverse()

    vnir_fwhm =  hdf_obj.attrs.get('List_Fwhm_Vnir')[6:].tolist()
    swir_fwhm   =hdf_obj.attrs.get('List_Fwhm_Swir')[:-3].tolist()
    vnir_fwhm.reverse()
    swir_fwhm.reverse()

    #Specific output file path
    header_dict = envi_header_dict()
    header_dict['lines']= vnir_data.shape[0]
    header_dict['samples']= vnir_data.shape[2]
    header_dict['bands']= len(vnir_waves + swir_waves)
    header_dict['wavelength']= vnir_waves + swir_waves
    header_dict['fwhm']= vnir_fwhm + swir_fwhm
    header_dict['interleave']= 'bsq'
    header_dict['data type'] = 12
    no_data = 2**16 -1
    header_dict['data ignore value'] = no_data
    header_dict['wavelength units'] = "nanometers"
    header_dict['map info'] = map_info

    #Create mask to change no data from 0 to -9999
    mask = (vnir_data[:,20,:]==0)&(vnir_data[:,10,:]==0)&(swir_data[:,100,:]==0)

    #Write bands to file
    writer = WriteENVI(output_name,header_dict)
    band_num =0

    for i in range(65,5,-1):
        band =  vnir_data[:,i,:]
        band[mask] = no_data
        writer.write_band(band,band_num)
        band_num+=1
        progbar(band_num, header_dict['bands'])

    for i in range(169,-1,-1):
        band =  swir_data[:,i,:]
        band[mask] = no_data
        writer.write_band(band,band_num)
        band_num+=1
        progbar(band_num, header_dict['bands'])

    # TODO: Add export of observables
    writer.close()
    print("\nExport complete.")

if __name__== "__main__":
    main()
