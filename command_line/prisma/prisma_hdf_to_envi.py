import argparse,warnings
import numpy as np,os
import hytools as ht
from hytools.io.envi import writeENVI,ENVI_header_from_hdf

warnings.filterwarnings("ignore")

def progbar(curr, total, full_progbar=100):
    frac = curr/total
    filled_progbar = round(frac*full_progbar)
    print('\r', '#'*filled_progbar + '-'*(full_progbar-filled_progbar), '[{:>7.2%}]'.format(frac), end='')

def main():
    parser = argparse.ArgumentParser(description = "Convert PRISMA H5 to merged VNIR-SWIR ENVI format")
    parser.add_argument("--img", help="Input image pathname",required=True, type = str)
    parser.add_argument("--out", help="Output directory",required=True, type = str)

    args = parser.parse_args()

    if not args.out.endswith("/"):
        args.out+="/"

    #Load data objects memory
    hyObj = ht.openPRISMA(args.img)
    print(hyObj.file_name)
    hyObj.load_data()
    
    output_name = args.out+ os.path.basename(os.path.splitext(args.img)[0]) 
    
    vnir_waves = hyObj.wavelengths_vnir[6:].tolist()
    swir_waves = hyObj.wavelengths_swir[:-3].tolist()
    vnir_waves.reverse()
    swir_waves.reverse()
    
    vnir_fwhm = hyObj.fwhm_vnir[6:].tolist()
    swir_fwhm = hyObj.fwhm_swir[:-3].tolist()
    vnir_fwhm.reverse()
    swir_fwhm.reverse()
    
    #Specific output file path
    header_dict = ENVI_header_from_hdf(hyObj)
    header_dict['bands']= len(vnir_waves + swir_waves)
    header_dict['wavelength']= vnir_waves + swir_waves
    header_dict['fwhm']= vnir_fwhm + swir_fwhm
    header_dict['interleave']= 'bsq'
    header_dict['data type'] = 12
    no_data = 2**16 -1
    header_dict['data ignore value'] = no_data
    
    #Create mask to change no data from 0 to -9999
    mask = (hyObj.vnir_data[:,20,:]==0)&(hyObj.vnir_data[:,10,:]==0)&(hyObj.swir_data[:,100,:]==0)
     
    #Write bands to file
    writer = writeENVI(output_name,header_dict)
    band_num =0
    
    for i in range(65,5,-1):
        band =  hyObj.vnir_data[:,i,:]
        band[mask] = no_data
        writer.write_band(band,band_num)
        band_num+=1
        progbar(band_num, header_dict['bands'])
    
    for i in range(169,-1,-1):
        band =  hyObj.swir_data[:,i,:]
        band[mask] = no_data
        writer.write_band(band,band_num)
        band_num+=1
        progbar(band_num, header_dict['bands'])

    # TODO: Add export of observables

    print("\nExport complete.")

if __name__== "__main__":
    main()



