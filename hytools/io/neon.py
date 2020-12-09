#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NEON AOP HDF opener
"""
import h5py



def open_neon(hy_obj, no_data = -9999):
    """Load and parse NEON formated HDF image into a HyTools file object.

    Args:
        src_file (str): pathname of input HDF file.
        no_data (float, optional): No data value. Defaults to -9999.
        load_obs (bool, optional): Map observables to memory. Defaults to False.

    Returns:
        HyTools file object: Populated HyTools file object.

    """

    hdf_obj = h5py.File(hy_obj.file_name,'r')
    hy_obj.base_key = list(hdf_obj.keys())[0]
    metadata = hdf_obj[hy_obj.base_key]["Reflectance"]["Metadata"]
    data = hdf_obj[hy_obj.base_key]["Reflectance"]["Reflectance_Data"]

    hy_obj.projection = metadata['Coordinate_System']['Coordinate_System_String'][()].decode("utf-8")
    hy_obj.map_info = metadata['Coordinate_System']['Map_Info'][()].decode("utf-8").split(',')
    hy_obj.transform = (float(hy_obj.map_info [3]),float(hy_obj.map_info [1]),0,float(hy_obj.map_info [4]),0,-float(hy_obj.map_info [2]))
    hy_obj.fwhm =  metadata['Spectral_Data']['FWHM'][()]
    hy_obj.wavelengths = metadata['Spectral_Data']['Wavelength'][()]
    hy_obj.wavelength_units = metadata['Spectral_Data']['Wavelength'].attrs['Units']
    hy_obj.lines = data.shape[0]
    hy_obj.columns = data.shape[1]
    hy_obj.bands = data.shape[2]
    hy_obj.no_data = no_data
    hy_obj.observables = {'path_length': ['Ancillary_Imagery','Path_Length'],
                        'sensor_az': ['to-sensor_Azimuth_Angle'],
                        'sensor_zn': ['to-sensor_Zenith_Angle'],
                        'solar_az': ['Logs','Solar_Azimuth_Angle'],
                        'solar_zn': ['Logs','Solar_Zenith_Angle'],
                        'slope': ['Ancillary_Imagery','Slope'],
                        'aspect':['Ancillary_Imagery','Aspect']}

    return hy_obj
