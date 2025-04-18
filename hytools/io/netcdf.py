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
import h5netcdf
import numpy as np
from .envi import parse_envi_header, WriteENVI

unit_dict = {'nm':'nanometers'}
utm_zone_dict = {'N':'North','S':'South'}

def open_netcdf(hy_obj, sensor,anc_path = {}, glt_path = {}):
    """Load and parse NASA formatted NetCDF AVIRIS/EMIT image into a HyTools file object.

    Args:
        HyTools file object: Populated HyTools file object.
        sensor (str): sensor name for reading, either 'emit' (EMIT) or 'ncav' (AVIRIS)
        anc_path (dict): Dictionary with pathnames and band numbers of ancillary datasets.
        glt_path (list): Dictionary with pathnames and band numbers of external GLT datasets.
    Returns:
        HyTools file object: Populated HyTools file object.

    """

    nc4_obj = h5py.File(hy_obj.file_name,'r')
    hy_obj.base_key = list(nc4_obj.keys())[0]

    metadata = nc4_obj.attrs
    if sensor=='AV':
        data = nc4_obj["reflectance"]["reflectance"]
        hy_obj.fwhm =  nc4_obj['reflectance']['fwhm'][()]
        hy_obj.wavelengths = nc4_obj['reflectance']['wavelength'][()]
        hy_obj.wavelength_units = unit_dict[nc4_obj['reflectance']['wavelength'].attrs['units'].decode("utf-8")]

        hy_obj.lines = data.shape[1]
        hy_obj.columns = data.shape[2]
        hy_obj.bands = data.shape[0]

    elif sensor == 'EMIT':
        data = nc4_obj["reflectance"]
        hy_obj.fwhm =  nc4_obj['sensor_band_parameters']['fwhm'][()]
        hy_obj.wavelengths = nc4_obj['sensor_band_parameters']['wavelengths'][()]
        hy_obj.wavelength_units = unit_dict[nc4_obj['sensor_band_parameters']['wavelengths'].attrs['units'].decode("utf-8")]
        hy_obj.lines = data.shape[0]
        hy_obj.columns = data.shape[1]
        hy_obj.bands = data.shape[2]

        hy_obj.bad_bands =  np.array(1-nc4_obj['sensor_band_parameters']['good_wavelengths'][()]).astype(np.bool)

    hy_obj.no_data = data.attrs['_FillValue'][0]

    hy_obj.anc_path = anc_path

    if bool(glt_path)==True:
        hy_obj.glt_path = glt_path

        glt_header_file = os.path.splitext(glt_path[list(glt_path.keys())[0]][0])[0] + ".hdr"
        glt_header=parse_envi_header(glt_header_file)
        hy_obj.map_info = glt_header["map info"]
        hy_obj.lines_glt = glt_header["lines"]
        hy_obj.columns_glt = glt_header["samples"]

        hy_obj.transform = (float(glt_header["map info"][3]),float(glt_header["map info"][5]),0,
                            float(glt_header["map info"][4]),0,-float(glt_header["map info"][6]))

        if "coordinate system string" in glt_header:
            hy_obj.projection = glt_header["coordinate system string"]
        else:
            hy_obj.projection = ''

    else:
        if sensor == 'EMIT':

            hy_obj.glt_path = { "glt_x": ["location","glt_x"],
                                "glt_y": ["location","glt_y"]}
            hy_obj.projection = metadata['spatial_ref'].decode("utf-8")
            geotransform = nc4_obj.attrs['geotransform'][()]
            hy_obj.map_info = ['Geographic Lat/Lon','1','1',
                               str(geotransform[0]),str(geotransform[3]),
                               str(geotransform[1]),str(-geotransform[5]),
                               'WGS-84']
            hy_obj.transform = tuple(metadata['geotransform'][()])
            glt_x = nc4_obj['location']['glt_x']

            hy_obj.lines_glt = glt_x.shape[0]
            hy_obj.columns_glt = glt_x.shape[1]

        elif sensor == 'AV':
            hy_obj.projection = nc4_obj['transverse_mercator'].attrs['spatial_ref'].decode("utf-8")
            geotransform = [ float(x) for x in nc4_obj['transverse_mercator'].attrs['GeoTransform'].decode("utf-8").split(' ')]
            utm_zone_tag=((hy_obj.projection).split('UTM zone ')[1]).split('",GEOGCS')[0]
            hy_obj.map_info = ['UTM','1','1',
                               str(geotransform[0]),str(geotransform[3]),
                               str(geotransform[1]),str(-geotransform[5]),
                               utm_zone_tag[:-1],utm_zone_dict[utm_zone_tag[-1]],'WGS-84']
            hy_obj.transform = tuple(geotransform)

            hy_obj.lines_glt = hy_obj.lines
            hy_obj.columns_glt = hy_obj.columns
            hy_obj.glt_path = { "glt_x": ["geolocation_lookup_table","sample"],
                                "glt_y": ["geolocation_lookup_table","line"]}            

    return hy_obj

def set_wavelength_meta(nc4_obj,header_dict,glt_bool):
    file_type = (header_dict['file_type']).lower()

    if file_type=="ncav" or (file_type=="emit" and glt_bool==True):
        gp=nc4_obj.create_group("reflectance")
        wavelength_var=nc4_obj.create_variable("/reflectance/wavelength",("wavelength",),
                                               data=header_dict['wavelength'],
                                               dtype=np.float32)
        fwhm_var = nc4_obj.create_variable("/reflectance/fwhm",("wavelength",),
                                           data=header_dict['fwhm'],
                                           dtype=np.float32)
    elif file_type=="emit":
        if glt_bool: # handled in above codes
            pass
        else: # do not warp with GLT
            nc4_obj.dimensions["bands"]=header_dict['bands']
            wavelength_var=nc4_obj.create_variable("/sensor_band_parameters/wavelengths",("bands",),
                                                   data=np.array(header_dict['wavelength']),
                                                   dtype=np.float32)
            fwhm_var = nc4_obj.create_variable("/sensor_band_parameters/fwhm", ("bands",),
                                               data=header_dict['fwhm'],
                                               dtype=np.float32)


def write_netcdf_refl_meta(nc4_obj,header_dict,glt_bool):
    set_wavelength_meta(nc4_obj,header_dict,glt_bool)
    write_netcdf_meta(nc4_obj,header_dict,glt_bool)

class WriteNetCDF(WriteENVI):
    """Iterator class for writing to a NetCDF data file.
        The class inherites all the write functionss from WriteENVI: write pixel, line, band, chunk, etc. 
    """
    def __init__(self,output_name, header_dict, attr_dict, glt_bool, type_tag, band_name=None):
        """
        Args:
            output_name (str): Pathname of output ENVI data file.
            header_dict (dict): Dictionary containing ENVI header information.

        Returns:
            None.

        """

        if type_tag=="reflectance": # for reflectance
            self.header_dict = header_dict
            self.output_name = output_name
            self.file_type = header_dict['file_type'].lower()

            self.nc4_obj = h5netcdf.File(output_name, "w")

            write_netcdf_refl_meta(self.nc4_obj,header_dict,glt_bool)
            if self.file_type in ["ncav","envi"]:
                self.interleave = "bsq"
                self.data = self.nc4_obj.create_variable("/reflectance/reflectance",
                                                         ("wavelength","northing","easting"),
                                                         np.float32,
                                                         chunks=(2,256,256),
                                                         compression='gzip')
                self.data.attrs["grid_mapping"] =  "projection"
            elif self.file_type == "emit":
                if glt_bool:
                    self.interleave = "bsq"
                    self.data = self.nc4_obj.create_variable("/reflectance/reflectance",
                                                             ("wavelength","northing","easting"),
                                                             np.float32,
                                                             chunks=(1,256,256),
                                                             compression='gzip')
                    self.data.attrs["grid_mapping"] =  "projection"
                else:
                    self.interleave = "bip"
                    self.data = self.nc4_obj.create_variable("reflectance",
                                                             ("downtrack","crosstrack","bands"),
                                                             np.float32,
                                                             chunks=(256,256,2),
                                                             compression='gzip')

            self.data.attrs["_FillValue"]=-9999.0
            self.external_nc_attrs(attr_dict)
        elif type_tag=="mask":  # for masks
            self.interleave = "bsq"
            self.header_dict = header_dict
            self.file_type = header_dict['file_type'].lower()
            self.nc4_obj = h5netcdf.File(output_name, "r+")

            if self.file_type in ["ncav","envi"]:
                self.data = self.nc4_obj.create_variable(f"/masks/{band_name}",
                                                         ("northing","easting"),
                                                         np.uint8,
                                                         chunks=(256,256),
                                                         compression='gzip')
                self.data.attrs["grid_mapping"] =  "projection"
            elif self.file_type == "emit":
                if glt_bool:
                    self.data = self.nc4_obj.create_variable(f"/masks/{band_name}",
                                                             ("northing","easting"),
                                                             np.uint8,
                                                             chunks=(256,256),
                                                             compression='gzip')
                    self.data.attrs["grid_mapping"] =  "projection"
                else:
                    self.data = self.nc4_obj.create_variable(f"/masks/{band_name}",
                                                             ("downtrack","crosstrack"),
                                                             np.uint8,
                                                             chunks=(256,256),
                                                             compression='gzip')
            self.data.attrs["_FillValue"]=255
            self.external_nc_attrs(attr_dict)
        elif type_tag=="trait":
            self.interleave = "bsq"
            self.file_type = header_dict['file_type'].lower()
            self.nc4_obj = h5netcdf.File(output_name, "w")

            self.nc4_obj.dimensions["bands"]=2

            self.interleave = "bsq"

            write_netcdf_meta(self.nc4_obj,header_dict,glt_bool)
            if self.file_type in ["ncav","envi"]:
                self.data = self.nc4_obj.create_variable(f"/{band_name}/stack",
                                                         ("bands","northing","easting"),
                                                         np.float32,
                                                         chunks=(1,256,256),
                                                         compression='gzip')
                self.data.attrs["grid_mapping"] =  "projection"

            elif self.file_type == "emit":
                if glt_bool:
                    self.data = self.nc4_obj.create_variable(f"/{band_name}/stack",
                                                             ("bands","northing","easting"),
                                                             np.float32,
                                                             chunks=(1,256,256),
                                                             compression='gzip')
                    self.data.attrs["grid_mapping"] =  "projection"
                else:
                    self.data = self.nc4_obj.create_variable(f"/{band_name}/stack",
                                                             ("bands","downtrack","crosstrack"),
                                                             np.float32,
                                                             chunks=(1,256,256),
                                                             compression='gzip')

            self.data.attrs["band_names"] =  header_dict["band names"][:2]
            self.data.attrs["_FillValue"]=-9999.0

    def write_mask_band(self,band):
        self.data[:,:]  = band

    def write_mask_band_glt(self,band,glt_indices,fill_mask):
        tmp_band = np.ones(fill_mask.shape)*self.header_dict['data ignore value']
        tmp_band[fill_mask] = band[glt_indices]
        tmp_band[~fill_mask] = 255

        self.data[:,:] = tmp_band


    def write_glt_dataset(self,glt_x_arr,glt_y_arr,dim_x_name="ortho_x",dim_y_name="ortho_y"):
        var_glt_x = self.nc4_obj.create_variable("/location/glt_x",(dim_y_name,dim_x_name),
                                                 data=glt_x_arr,
                                                 dtype=np.int32,
                                                 chunks=(256,256),
                                                 compression='gzip')
        var_glt_y = self.nc4_obj.create_variable("/location/glt_y",(dim_y_name,dim_x_name),
                                                 data=glt_y_arr,
                                                 dtype=np.int32,
                                                 chunks=(256,256),
                                                 compression='gzip')

        var_glt_x.attrs["grid_mapping"] =  "projection"
        var_glt_y.attrs["grid_mapping"] =  "projection"

        var_glt_x.attrs["_FillValue"]=0
        var_glt_y.attrs["_FillValue"]=0


    def write_netcdf_band_glt(self,band,index,glt_indices,fill_mask):
        """
        Args:
            band (numpy.ndarray): Band array (lines,columns).
            index (int): Zero-based band index.
            glt_indices (numpy.ndarray,numpy.ndarray): Zero-based tuple indices.

        Returns:
            None.

        """

        tmp_band = np.ones(fill_mask.shape)*(-9999)
        tmp_band[fill_mask] = band[glt_indices]
        tmp_band[~fill_mask] = -9999

        if self.interleave == "bip":
            self.data[:,:,index]=tmp_band
        elif self.interleave == "bil":
            self.data[:,index,:]=tmp_band
        elif self.interleave == "bsq":
            self.data[index,:,:]=tmp_band

    def external_nc_attrs(self,attr_dict):

        if attr_dict is None:
            return

        for key in attr_dict:
            split_key = key.split('/')
            if len(split_key[0])==0:
                split_key.pop(0)
            if len(split_key)>1:
                group_path = '/'+'/'.join(split_key[:-1])
                self.nc4_obj[group_path].attrs[split_key[-1]]=str(attr_dict[key])
            else:
                self.nc4_obj.attrs[key]=str(attr_dict[key])


    def close(self):
        """Delete 
        """
        self.nc4_obj.close()

def write_netcdf_meta(nc4_obj,header_dict,glt_bool):

    file_type = (header_dict['file_type']).lower()

    if file_type=="ncav" or (file_type=="emit" and glt_bool==True):
        transform=header_dict['transform']

        nc4_obj.dimensions["northing"]=header_dict['lines'] #dim0
        nc4_obj.dimensions["easting"]=header_dict['samples'] #dim1

        tm_var = nc4_obj.create_variable("/projection",data=np.array([0]),dtype=np.uint8)
        tm_var.attrs["GeoTransform"]=' '.join([str(x) for x in header_dict['transform']])
        tm_var.attrs["crs_wkt"]=header_dict['projection']
        tm_var.attrs["spatial_ref"]=header_dict['projection']


    elif file_type=="emit":
        if glt_bool: # handled in above codes
            pass
        else: # do not warp with GLT
            loc_gp=nc4_obj.create_group("location")
            nc4_obj.dimensions["downtrack"]=header_dict['lines'] #dim0
            nc4_obj.dimensions["crosstrack"]=header_dict['samples'] #dim1

            nc4_obj.dimensions["ortho_y"]=header_dict['lines_glt']
            nc4_obj.dimensions["ortho_x"]=header_dict['samples_glt']

            nc4_obj.attrs["geotransform"]=' '.join([str(x) for x in header_dict['transform']])
            nc4_obj.attrs["spatial_ref"]=header_dict['projection']
            nc4_obj.attrs["spatialResolution"]=np.sqrt(header_dict['transform'][1]**2+header_dict['transform'][2]**2)

            tm_var = nc4_obj.create_variable("/projection",data=np.array([0]),dtype=np.uint8)
            tm_var.attrs["GeoTransform"]=' '.join([str(x) for x in header_dict['transform']])
            tm_var.attrs["crs_wkt"]=header_dict['projection']
            tm_var.attrs["spatial_ref"]=header_dict['projection']
