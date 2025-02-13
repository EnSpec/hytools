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

Base

TODO: Add corrections to ndi()

"""
import os
import json
import numpy as np
import h5py
import warnings
import sys
from .io.envi import envi_read_band,envi_read_pixels
from .io.envi import envi_read_line,envi_read_column,envi_read_chunk
from .io.envi import open_envi,parse_envi_header,envi_header_from_neon,envi_header_from_nc
from .io.neon import open_neon
from .io.netcdf import open_netcdf
from .brdf import apply_brdf_correct
from .glint import apply_glint_correct
from .brdf.kernels import calc_volume_kernel,calc_geom_kernel
from .topo import calc_cosine_i,apply_topo_correct
from .transform.resampling import *

warnings.filterwarnings("ignore")

class HyTools:
    """HyTools file object"""

    def __init__(self):
        """Constructor method
        """
        self.anc_path = {}
        self.ancillary = {}
        self.bad_bands = []
        self.bands = None
        self.base_key = None
        self.base_name = None
        self.brdf = {'type': None}
        self.glint= {'type': None}
        self.byte_order = None
        self.columns = None
        self.columns_glt = None
        self.corrections = []
        self.crs = None
        self.data = None
        self.dtype = None
        self.endianness = None
        self.file_name = None
        self.file_type = None
        self.fill_mask = None
        self.fwhm = []
        self.glt_path = {}
        self.glt_x = None
        self.glt_y = None
        self.hdf_obj  = None
        self.interleave = None
        self.lines = None
        self.lines_glt = None
        self.map_info = None
        self.mask = {}
        self.nc4_obj  = None
        self.no_data = None
        self.offset = 0
        self.projection = None
        self.resampler = {'type': None}
        self.shape = None
        self.topo = {'type': None}
        self.ulx = None
        self.uly = None
        self.wavelength_units = None
        self.wavelengths = []

    def read_file(self,file_name,file_type = 'envi',anc_path = None, ext = False, glt_path = None):
        self.file_name = file_name
        self.file_type = file_type

        if file_type == 'envi':
            open_envi(self,anc_path,ext)
        elif file_type == "neon":
            open_neon(self)
        elif file_type == "emit":
            open_netcdf(self,'EMIT',anc_path,glt_path)
        elif file_type == "ncav":
            open_netcdf(self,'AV',anc_path,glt_path)

        else:
            print("Unrecognized file type.")

        # Create a no data mask
        if self.bands>11:
          self.mask['no_data'] = self.get_wave(660) > 0.5*self.no_data
        else:
          self.mask['no_data'] = self.get_band(0) > 0.5*self.no_data

        #Match mask with ancillary mask
        if anc_path:

            if file_type == 'envi':
                ancillary = HyTools()
                ancillary.read_file(self.anc_path['solar_zn'][0],'envi')
                if not np.array_equal(self.mask['no_data'],ancillary.mask['no_data']):
                    print('Reflectance and ancillary no data extents do not match, combining no data masks.')
                    self.mask['no_data'] &= ancillary.mask['no_data']
                ancillary.close_data()
                del ancillary
            elif file_type == 'emit' and not self.anc_path['slope'][0].endswith('.nc'):
                ancillary = HyTools()
                ancillary.read_file(self.anc_path['slope'][0],'envi')
                if not np.array_equal(self.mask['no_data'],ancillary.mask['no_data']):
                    print('Reflectance and ancillary no data extents do not match, combining no data masks.')
                    self.mask['no_data'] &= ancillary.mask['no_data']
                ancillary.close_data()
                del ancillary

        self.base_name = os.path.basename(os.path.splitext(self.file_name)[0])

    def create_bad_bands(self,bad_regions):
        """Create bad bands mask, Good: True, bad : False.

        Args:
            bad_regions (list of lists): start and end values of wavelength
            regions considered bad. Wavelengths should be in the same units as
            data units. ex: [[350,400].....[2450,2500]].

        Returns:
            None.

        """

        bad_bands = []
        for wavelength in self.wavelengths:
            bad=False
            for start,end in bad_regions:
                bad = ((wavelength >= start) & (wavelength <=end)) or bad
            bad_bands.append(bad)
        self.bad_bands = np.array(bad_bands)


    def load_data(self, mode = 'r'):
        """Load data object to memory.

        Args:
            mode (str, optional): File read mode. Defaults to 'r'.
            offset (int, optional): Offset in bytes. Defaults to 0.

        Returns:
            None.

        """

        if self.file_type  == "envi":
            self.data = np.memmap(self.file_name,dtype = self.dtype, mode=mode,
                                  shape = self.shape,offset=self.offset)
        elif self.file_type  == "neon":
            self.hdf_obj = h5py.File(self.file_name,'r')
            self.data = self.hdf_obj[self.base_key]["Reflectance"]["Reflectance_Data"]
        elif self.file_type  == "emit":
            self.nc4_obj = h5py.File(self.file_name,'r')
            self.data = self.nc4_obj["reflectance"]
            self.glt_x = self.load_glt('glt_x')
            self.glt_y = self.load_glt('glt_y')
            self.fill_mask = (self.glt_x>0)
        elif self.file_type  == "ncav":
            self.nc4_obj = h5py.File(self.file_name,'r')
            self.data = self.nc4_obj["reflectance"]["reflectance"]


    def close_data(self):
        """Close data object.

        """
        if self.file_type  == "envi":
            del self.data
        elif self.file_type  == "neon":
            self.hdf_obj.close()
            self.hdf_obj = None
        elif self.file_type  == "emit" or self.file_type == "ncav":
            self.nc4_obj.close()
            self.nc4_obj = None
        self.data = None


    def iterate(self,by,chunk_size= (100,100),corrections = [],resample=False):
        """Create data Iterator.

        Args:
            by (str): Dimension along which to iterate: "line","column","band","chunk".
            chunk_size (tuple, optional): Two dimensional chunk size (Y,X).
                                          Applies only when "chunk" selected.
                                          Defaults to (100,100).

        Returns:
            Iterator class object: Data Iterator.

        """

        return Iterator(self,by,chunk_size,corrections =corrections,resample=resample)

    def wave_to_band(self,wave):
        """Return band index corresponding to input wavelength. Return closest band if
           not an exact match.

        Args:
            wave (float): Wavelength of band to be retrieved in image wavelength units.

        Returns:
            int: Band index.

        """

        if (wave  > self.wavelengths.max()) | (wave  < self.wavelengths.min()):
            print("Input wavelength outside image range!")
            band_num = None
        else:
            band_num = np.argmin(np.abs(self.wavelengths - wave))
        return band_num

    def get_band(self,index,corrections= [], mask =None):
        """
        Args:
            index (int): Zero-indexed band index.
            mask (str): Return masked values using named mask.
            corrections(list): Corrections to apply, will be applied in
            order listed.

        Returns:
            numpy.ndarray: A 2D (lines x columns) array or 1D if masked.

        """

        self.load_data()
        if self.file_type == "neon":
            band =  self.data[:,:,index]
        elif self.file_type == "emit":
            band =  self.data[:,:,index]
        elif self.file_type == "ncav":
            band =  self.data[index,:,:]
        elif self.file_type == "envi":
            band = envi_read_band(self.data,index,self.interleave)
            if self.endianness != sys.byteorder:
                band = band.byteswap()
        self.close_data()

        band = self.correct(band,'band',index,corrections)

        if mask:
            band = band[self.mask[mask]]

        return band


    def get_wave(self,wave,corrections= [],mask =None):
        """Return the band image corresponding to the input wavelength.
        If not an exact match the closest wavelength will be returned.

        Args:
            wave (float): Wavelength in image units.
            mask (str): Return masked values using named mask.


        Returns:
            numpy.ndarray: Band image array (line,columns).

        """

        if (wave  > self.wavelengths.max()) | (wave  < self.wavelengths.min()):
            print("Input wavelength outside wavelength range!")
            band = None
        else:
            band_num = np.argmin(np.abs(self.wavelengths - wave))
            band = self.get_band(band_num,corrections= corrections, mask=mask)
        return band

    def get_pixels(self,lines,columns,corrections= [],resample = False):
        """
        Args:
            lines (list): List of zero-indexed line indices.
            columns (list): List of zero-indexed column indices.

        Returns:
            numpy.ndarray: Pixel array (pixels,bands).

        """

        self.load_data()
        if self.file_type == "neon" or self.file_type == "emit":
            pixels = []
            for line,column in zip(lines,columns):
                pixels.append(self.data[line,column,:])
            pixels = np.array(pixels)
        elif self.file_type == "ncav":
            pixels = self.data[:,lines,columns]
        elif self.file_type == "envi":
            pixels = envi_read_pixels(self.data,lines,columns,self.interleave)
            if self.endianness != sys.byteorder:
                pixels = pixels.byteswap()
        self.close_data()

        pixels = self.correct(pixels,'pixels',
                         [lines,columns],corrections)

        if resample:
            pixels = pixels[np.newaxis,:,~self.bad_bands]
            pixels = apply_resampler(self,pixels)[0,:,:]

        return pixels

    def get_line(self,index, corrections= [],resample = False):
        """
        Args:
            index (int): Zero-indexed line index.

        Returns:
            numpy.ndarray: Line array (columns, bands).

        """

        self.load_data()
        if self.file_type == "neon" or self.file_type == "emit":
            line = self.data[index,:,:]
        elif self.file_type == "ncav":
            line = np.moveaxis(self.data[:,index,:],0,1)
        elif self.file_type == "envi":
            line = envi_read_line(self.data,index,self.interleave)
            if self.endianness != sys.byteorder:
                line = line.byteswap()
        self.close_data()

        line = self.correct(line,'line',index,corrections)

        if resample:
            line = line[np.newaxis,:,~self.bad_bands]
            line = apply_resampler(self,line)[0,:,:]

        return line

    def get_column(self,index,corrections = [],resample = False):
        """
        Args:
            index (int): Zero-indexed column index.

        Returns:
            numpy.ndarray: Column array (lines, bands).

        """

        self.load_data()
        if self.file_type == "neon" or self.file_type == "emit":
            column = self.data[:,index,:]
        elif self.file_type == "ncav":
            column = np.moveaxis(self.data[:,:,index],0,1)
        elif self.file_type == "envi":
            column = envi_read_column(self.data,index,self.interleave)
            if self.endianness != sys.byteorder:
                column = column.byteswap()
        self.close_data()

        column = self.correct(column,'column',index,corrections)

        if resample:
            column = column[:,np.newaxis,~self.bad_bands]
            column = apply_resampler(self,column)[:,0,:]

        return column

    def get_chunk(self,col_start,col_end,line_start,line_end, corrections= [],resample = False):
        """
        Args:
            col_start (int): Chunk starting column.
            col_end (int): Noninclusive chunk ending column index.
            line_start (int): Chunk starting line.
            line_end (int): Noninclusive chunk ending line index.
            corrections(list): Corrections to apply, will be applied in
            order listed.
            resample (bool): Resample wavelengths. Defaults to False.

        Returns:
            numpy.ndarray: Chunk array (line_end-line_start,col_end-col_start,bands).

        """

        self.load_data()
        if self.file_type == "neon" or self.file_type == "emit":
            chunk = self.data[line_start:line_end,col_start:col_end,:]
        elif self.file_type == "ncav":
            chunk = np.moveaxis(self.data[:,line_start:line_end,col_start:col_end],0,-1)
        elif self.file_type == "envi":
            chunk =  envi_read_chunk(self.data,col_start,col_end,
                                     line_start,line_end,self.interleave)
            if self.endianness != sys.byteorder:
                chunk = chunk.byteswap()
        self.close_data()

        chunk = self.correct(chunk,'chunk',
                        [col_start,col_end,line_start,line_end],
                        corrections)
        if resample:
            chunk = apply_resampler(self,chunk[:,:,~self.bad_bands])
        return chunk

    def correct(self,data,dimension,index,corrections):
        for correction in corrections:
            if correction == 'topo':
                data = apply_topo_correct(self,data,dimension,index)
            elif correction == 'brdf':
                data = apply_brdf_correct(self,data,dimension,index)
            elif correction == 'glint':
                data = apply_glint_correct(self,data,dimension,index)
        return data

    def get_anc(self,anc,radians = True,mask = None):
        """Read ancillary datasets to memory.

        Args:
            anc (str): Ancillary dataset name.
            radians (bool, optional): Convert angular measures to radians. Defaults to True.

        Returns:
            anc_data (numpy.ndarray)

        """

        angular_anc = ['slope','sensor_az','sensor_zn','aspect','solar_zn','solar_az']

        if self.file_type == "envi":
            ancillary = HyTools()
            ancillary.read_file(self.anc_path[anc][0],'envi')
            ancillary.load_data()
            anc_data = np.copy(ancillary.get_band(self.anc_path[anc][1]))
            if ancillary.endianness != sys.byteorder:
                anc_data = anc_data.byteswap()
            ancillary.close_data()

        elif self.file_type == "neon":
            keys = self.anc_path[anc]

            if len(keys)==2 and isinstance(keys[-1],int):
                ancillary = HyTools()
                ancillary.read_file(self.anc_path[anc][0],'envi')
                ancillary.load_data()
                anc_data = np.copy(ancillary.get_band(self.anc_path[anc][1]))
                if ancillary.endianness != sys.byteorder:
                    anc_data = anc_data.byteswap()
                ancillary.close_data()
                
            else:
                hdf_obj = h5py.File(self.file_name,'r')

                metadata = hdf_obj[self.base_key]["Reflectance"]["Metadata"]
                for key in keys:
                    metadata = metadata[key]
                anc_data = metadata[()]
                hdf_obj.close()

            #Make solar geometry into 2D array
            if anc in ['solar_zn','solar_az']:
                anc_data = np.ones((self.lines, self.columns)) * anc_data

        elif self.file_type == "emit" or self.file_type == "ncav":
            if bool(self.anc_path)==False:
                return None

            else:    
                if (self.anc_path[anc][0]).endswith('nc'):
                    nc4_anc_obj = h5py.File(self.anc_path[anc][0],'r')
                    
                    if self.file_type == "emit":
                        anc_data = nc4_anc_obj['obs'][()][:,:,self.anc_path[anc][1]]
                    elif self.file_type == "ncav":
                        anc_data_raw = nc4_anc_obj['observation_parameters'][self.anc_path[anc][1]][()]
                        obs_glt_x = np.abs(nc4_anc_obj['geolocation_lookup_table']['sample'][()]) # some values in the GLT are negative for unknown reason
                        obs_glt_y = np.abs(nc4_anc_obj['geolocation_lookup_table']['line'][()])
                        anc_data = np.zeros(obs_glt_x.shape)
                        anc_data[obs_glt_x<=0] = nc4_anc_obj['observation_parameters'][self.anc_path[anc][1]].attrs['_FillValue'][0]  # -9999
                        data_mask_to_fill = (obs_glt_x>0)
                        anc_data[data_mask_to_fill] = anc_data_raw[obs_glt_y[data_mask_to_fill].astype(int)-1,obs_glt_x[data_mask_to_fill].astype(int)-1]
                        
                    nc4_anc_obj.close()
                else:
                    ancillary = HyTools()
                    ancillary.read_file(self.anc_path[anc][0],'envi')
                    ancillary.load_data()
                    anc_data = np.copy(ancillary.get_band(self.anc_path[anc][1]))
                    if ancillary.endianness != sys.byteorder:
                        anc_data = anc_data.byteswap()
                    ancillary.close_data()

        if radians and (anc in angular_anc):
            anc_data= np.radians(anc_data)


        if mask:
            anc_data = anc_data[self.mask[mask]]

        return anc_data

    def load_anc(self,anc,radians = True):
        self.ancillary[anc] = self.get_anc(self,anc,radians)


    def load_glt(self,glt):
        # check if GLT inside nc is used
        if self.glt_path[glt][0]=='location':
            glt_data = self.nc4_obj[self.glt_path[glt][0]][self.glt_path[glt][1]][()]
        else:
            glt_img = HyTools()
            glt_img.read_file(self.glt_path[glt][0],'envi')
            glt_img.load_data()
            glt_data = np.copy(glt_img.get_band(self.glt_path[glt][1]))
            if glt_img.endianness != sys.byteorder:
                glt_data = glt_data.byteswap()
            glt_img.close_data()

        return glt_data


    def volume_kernel(self,kernel):
        """Calculate volume scattering kernel.
        """

        return calc_volume_kernel(self.get_anc('solar_az'), self.get_anc('solar_zn'),
                                  self.get_anc('sensor_az'), self.get_anc('sensor_zn'),
                                               kernel)

    def geom_kernel(self,kernel,b_r=1.,h_b =2.):
        """Calculate volume scattering kernel.
        """

        return calc_geom_kernel(self.get_anc('solar_az'),self.get_anc('solar_zn'),
                                self.get_anc('sensor_az'),self.get_anc('sensor_zn'),
                                kernel,b_r=b_r,h_b =h_b)

    def cosine_i(self):
        """ Calculate the cosine of the solar incidence angle. Assumes
        path to required ancillary datasets have been specified.

        Returns:
            cos_i numpy.ndarray: Cosine of solar incidence angle.

        """

        cos_i = calc_cosine_i(self.get_anc('solar_zn'), self.get_anc('solar_az'),
                          self.get_anc('aspect') ,self.get_anc('slope'))
        return cos_i

    def ndi(self,wave1= 850,wave2 = 660,mask = None):
        """ Calculate normalized difference index.
            Defaults to NDVI. Assumes input wavelengths are in
            nanometers

        Args:
            wave1 (int,float): Wavelength of first band. Defaults to 850.
            wave2 (int,float): Wavelength of second band. Defaults to 660.
            mask (bool): Mask data

        Returns:
            ndi numpy.ndarray:

        """

        wave1 = self.get_wave(wave1)
        wave2 = self.get_wave(wave2)
        ndi = (wave1-wave2)/(wave1+wave2)

        if mask:
            ndi = ndi[self.mask[mask]]
        return ndi


    def set_mask(self,mask,name):
        """Generate mask using masking function which takes a HyTools object as
        an argument.
        """
        self.mask[name] = mask

    def gen_mask(self,masker,name,args = None):
        """Generate mask using masking function which takes a HyTools object as
        an argument.
        """
        if args:
            self.mask[name] = masker(self,args)
        else:
            self.mask[name] = masker(self)

    def do(self,function,args = None):
        """Run a function and return the results.

        """
        if args:
            return function(self, args)
        else:
            return function(self)


    def get_header(self,warp_glt = False):
        """ Return header dictionary

        """
        if self.file_type == "neon":
            header_dict = envi_header_from_neon(self)
        elif self.file_type == "emit" or self.file_type == "ncav":
            header_dict = envi_header_from_nc(self,warp_glt = warp_glt)    
        elif self.file_type == "envi":
            header_dict = parse_envi_header(self.header_file)
        return header_dict


    def load_coeffs(self, coeff_file,kind):
        with open(coeff_file, 'r') as outfile:
            if kind == 'brdf':
                self.brdf = json.load(outfile, cls =Decoder)
            elif  kind == 'topo':
                self.topo = json.load(outfile, cls =Decoder)


class Iterator:
    """Iterator class
    """

    def __init__(self,hy_obj,by,chunk_size = None,corrections = [],resample = False):
        """
        Args:
            hy_obj (Hytools object): Populated Hytools file object.
            by (str): Iterator slice dimension: "line", "column", "band"",chunk".
            chunk_size (tuple, optional): Chunk size. Defaults to None.

        Iterator cannot be pickled when reading HDF files.

        Returns:
            None.

        """

        self.chunk_size= chunk_size
        self.by = by
        self.current_column = -1
        self.current_line = -1
        self.current_band = -1
        self.complete = False
        self.hy_obj = hy_obj
        self.resample = resample
        self.corrections = corrections


    def read_next(self):
        """ Return next line/column/band/chunk.
        """

        if self.by == "line":
            self.current_line +=1
            if self.current_line == self.hy_obj.lines-1:
                self.complete = True
            subset = self.hy_obj.get_line(self.current_line,
                                            corrections =self.corrections,
                                            resample = self.resample)
        elif self.by == "column":
            self.current_column +=1
            if self.current_column == self.hy_obj.columns-1:
                self.complete = True
            subset = self.hy_obj.get_column(self.current_column,
                                            corrections =self.corrections,
                                            resample = self.resample)
        elif self.by == "band":
            self.current_band +=1
            if self.current_band == self.hy_obj.bands-1:
                self.complete = True
            subset = self.hy_obj.get_band(self.current_band,
                                            corrections =self.corrections)

        elif self.by == "chunk":
            if self.current_column == -1:
                self.current_column +=1
                self.current_line +=1
            else:
                self.current_column += self.chunk_size[1]
            if self.current_column >= self.hy_obj.columns:
                self.current_column = 0
                self.current_line += self.chunk_size[0]

            y_start = self.current_line
            y_end = self.current_line + self.chunk_size[0]
            if y_end >= self.hy_obj.lines:
                y_end = self.hy_obj.lines
            x_start = self.current_column
            x_end = self.current_column + self.chunk_size[1]
            if x_end >= self.hy_obj.columns:
                x_end = self.hy_obj.columns

            if (y_end == self.hy_obj.lines) and (x_end == self.hy_obj.columns):
                self.complete = True

            subset = self.hy_obj.get_chunk(x_start,x_end, y_start,y_end,
                                            corrections =self.corrections,
                                            resample = self.resample)

        elif self.by == "glt_line":
            self.current_line +=1
            if self.current_line == self.hy_obj.lines-1:
                self.complete = True
            valid_mask=self.hy_obj.fill_mask[self.current_line,:]
            
            valid_subset = self.hy_obj.get_pixels(
                                            self.hy_obj.glt_y[self.current_line,valid_mask]-1,self.hy_obj.glt_x[self.current_line,valid_mask]-1,
                                            corrections = self.corrections,
                                            resample = self.resample)

            subset = np.full((self.hy_obj.columns_glt,valid_subset.shape[1]),-9999).astype(np.float32)
            subset[valid_mask,:] = valid_subset

        return subset

    def reset(self):
        """Reset counters.
        """
        self.current_column = -1
        self.current_line = -1
        self.current_band = -1
        self.complete = False


class Decoder(json.JSONDecoder):
    def decode(self, s):
        result = super().decode(s)  # result = super(Decoder, self).decode(s) for Python 2.x
        return self._decode(result)

    def _decode(self, o):
        if isinstance(o, str):
            try:
                return int(o)
            except ValueError:
                return o
        elif isinstance(o, dict):
            return {k: self._decode(v) for k, v in o.items()}
        elif isinstance(o, list):
            return [self._decode(v) for v in o]
        else:
            return o
