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

Tanager Reflectance HDF opener
"""
import h5py
import numpy as np

unit_dict = {'nm':'nanometers'}
utm_south_dict = {False:'North',True:'South'}

def open_tanager(hy_obj, anc_path = {}, no_data = -9999,):
    """Load and parse NEON formated HDF image into a HyTools file object.

    Args:
        src_file (str): pathname of input HDF file.
        no_data (float, optional): No data value. Defaults to -9999.

    Returns:
        HyTools file object: Populated HyTools file object.

    """

    hdf_obj = h5py.File(hy_obj.file_name,'r')
    hy_obj.base_key = list(hdf_obj.keys())[0]
    metadata = hdf_obj["HDFEOS INFORMATION"]["StructMetadata.0"][()].split()   #["Reflectance"]["Metadata"]
    list_metadata = [m.decode('utf-8') for m in metadata]
    data = hdf_obj[hy_obj.base_key]["GRIDS"]["HYP"]["Data Fields"]["surface_reflectance"]

    hy_obj.lines = data.shape[1]
    hy_obj.columns = data.shape[2]
    hy_obj.bands = data.shape[0]
    #print(hy_obj.lines,hy_obj.columns,hy_obj.bands)

    ulc=[i for i in list_metadata if "UpperLeftPointMtrs" in i][0]
    ulx = float(ulc.split('=(')[-1].replace(')','').split(',')[0])
    uly = float(ulc.split('=(')[-1].replace(')','').split(',')[1])
    lrc=[i for i in list_metadata if "LowerRightMtrs" in i][0]
    lrx = float(lrc.split('=(')[-1].replace(')','').split(',')[0])
    lry = float(lrc.split('=(')[-1].replace(')','').split(',')[1])

    utm_zone_tag = [i for i in list_metadata if "ZoneCode" in i][0]
    utm_zone_code = int(utm_zone_tag.split('=')[-1])
    bool_south = utm_zone_code < 0

    #assuming no rotation
    xRes = (lrx - ulx) / hy_obj.columns
    yRes = (uly - lry) / hy_obj.lines

    hy_obj.transform = (ulx,xRes,0,uly,0,-yRes)
    hy_obj.map_info = ['UTM','1.0','1.0',
                               str(ulx),str(uly),
                               str(xRes),str(yRes),
                               str(np.abs(utm_zone_code)),utm_south_dict[bool_south],
                               'WGS-84',
                               'units=Meters']

    hy_obj.fwhm =  data.attrs['fwhm'][()]
    hy_obj.wavelengths = data.attrs['wavelengths'][()]
    hy_obj.wavelength_units = unit_dict[data.attrs['wavelengths_units']]
    hy_obj.bad_bands = data.attrs['good_wavelengths'][()]

    hy_obj.no_data = no_data
    hy_obj.anc_path = {'path_length': ['sensor_to_ground_path_length'],
                        'sensor_az': ['sensor_azimuth'],
                        'sensor_zn': ['sensor_zenith'],
                        'solar_az': ['sun_azimuth'],
                        'solar_zn': ['sun_zenith'],
                        'aod': ['aerosol_optical_depth'],
                        'beta_cirrus_mask': ['beta_cirrus_mask'],
                        'beta_cloud_mask': ['beta_cloud_mask'],
                        'nodata_pixels': ['nodata_pixels'],
                        'UTC time':['time'],
                        'water_vapor': ['column_water_vapor']}

    if anc_path is not None:
        #override all paths
        for key in anc_path:
            hy_obj.anc_path[key] = anc_path[key]

    return hy_obj
