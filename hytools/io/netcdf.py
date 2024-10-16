# -*- coding: utf-8 -*-
"""
HyTools:  Hyperspectral image processing library
Copyright (C) 2021 University of Wisconsin

Authors: Adam Chlus, Zhiwei Ye, Philip Townsend.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

NASA NetCDF opener
"""
import os
import h5py
import numpy as np
from .envi import parse_envi_header

unit_dict = {'nm':'nanometers'}

def open_netcdf(hy_obj, anc_path = {}, glt_path = {}):
    """Load and parse NASA formatted NetCDF image into a HyTools file object.

    Args:


    Returns:
        HyTools file object: Populated HyTools file object.

    """

    nc4_obj = h5py.File(hy_obj.file_name,'r')
    hy_obj.base_key = list(nc4_obj.keys())[0]
    
    metadata = nc4_obj.attrs  
    data = nc4_obj["reflectance"]   

    hy_obj.fwhm =  nc4_obj['sensor_band_parameters']['fwhm'][()]
    hy_obj.wavelengths = nc4_obj['sensor_band_parameters']['wavelengths'][()]
    hy_obj.wavelength_units = unit_dict[nc4_obj['sensor_band_parameters']['wavelengths'].attrs['units'].decode("utf-8")]
    hy_obj.lines = data.shape[0]
    hy_obj.columns = data.shape[1]


    hy_obj.bands = data.shape[2]
    hy_obj.bad_bands =  np.array(1-nc4_obj['sensor_band_parameters']['good_wavelengths'][()]).astype(np.bool)  #np.array([False for band in range(hy_obj.bands)])
    hy_obj.no_data = data.attrs['_FillValue'][0]

    hy_obj.anc_path = anc_path

    if bool(glt_path)==True:
        hy_obj.glt_path = glt_path

        glt_header_file = os.path.splitext(glt_path[list(glt_path.keys())[0]][0])[0] + ".hdr"
        glt_header=parse_envi_header(glt_header_file)

        hy_obj.map_info = glt_header["map info"]
        hy_obj.lines_glt = glt_header["lines"]
        hy_obj.columns_glt = glt_header["samples"]  

    else:    
        hy_obj.glt_path = { "glt_x": ["location","glt_x"],
                            "glt_y": ["location","glt_y"]}
        hy_obj.projection = metadata['spatial_ref'].decode("utf-8")
        geotransform = nc4_obj.attrs['geotransform'][()]
        hy_obj.map_info = ['Geographic Lat/Lon','1','1',str(geotransform[0]),str(geotransform[3]),str(geotransform[1]),str(-geotransform[5]),'WGS-84']
        hy_obj.transform = tuple(metadata['geotransform'][()])
        glt_x = nc4_obj['location']['glt_x']
    
        hy_obj.lines_glt = glt_x.shape[0]
        hy_obj.columns_glt = glt_x.shape[1]

    return hy_obj
