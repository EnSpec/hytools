# -*- coding: utf-8 -*-
"""
Base
"""
import numpy as np
import h5py
import ray
from .io.envi import envi_read_band,envi_read_pixels
from .io.envi import envi_read_line,envi_read_column,envi_read_chunk
from .io.envi import open_envi
from .io.neon import open_neon
from .correction.brdf import calc_volume_kernel,calc_geom_kernel

@ray.remote
class HyTools:
    """HyTools file object"""

    def __init__(self):
        """Constructor method
        """
        self.file_type = None
        self.interleave = None
        self.file_name = None
        self.shape = None
        self.lines = None
        self.columns = None
        self.bands = None
        self.wavelengths = None
        self.fwhm = []
        self.bad_bands = []
        self.no_data = None
        self.map_info = None
        self.crs = None
        self.ulx = None
        self.uly = None
        self.dtype = None
        self.data = None
        self.header_dict = None
        self.projection = None
        self.byte_order = None
        self.wavelength_units = None
        self.hdf_obj  = None
        self.offset = 0
        self.base_key = None
        self.observables = None
        self.mask = None


    def read_file(self,file_name,file_type):

        self.file_name = file_name
        self.file_type = file_type

        if file_type == 'envi':
            open_envi(self)
        elif file_type == "neon":
            open_neon(self)
        elif file_type == "prisma":
            print("PRISMA not yet supported.")
        else:
            print("Unrecognized file type.")

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

    def close_data(self):
        """Close data object.

        """
        if self.file_type  == "envi":
            del self.data
        elif self.file_type  == "neon":
            self.hdf_obj.close()
            self.hdf_obj = None
        self.data = None


    def iterate(self,by,chunk_size= (100,100)):
        """Create data Iterator.

        Args:
            by (str): Dimension along which to iterate: "line","column","band","chunk".
            chunk_size (tuple, optional): Two dimensional chunk size (Y,X).
                                          Applies only when "chunk" selected.
                                          Defaults to (100,100).

        Returns:
            Iterator class object: Data Iterator.

        """

        return Iterator(self,by,chunk_size)

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

    def get_band(self,index,mask_values =False):
        """
        Args:
            index (int): Zero-indexed band index.
            mask_values (bool): Return only masked values,
            requires that a mask is set.

        Returns:
            numpy.ndarray: A 2D (lines x columns) array or 1D if masked.

        """

        self.load_data()
        if self.file_type == "neon":
            band =  self.data[:,:,index]
        elif self.file_type == "envi":
            band = envi_read_band(self.data,index,self.interleave)
        self.close_data()

        if mask_values and isinstance(self.mask, np.ndarray):
            band = band[self.mask]

        return band


    def get_wave(self,wave,mask_values =False):
        """Return the band image corresponding to the input wavelength.
        If not an exact match the closest wavelength will be returned.

        Args:
            wave (float): DESCRIPTION.
            mask_values (bool): Return only masked values,
            requires that a mask is set.

        Returns:
            numpy.ndarray: Band image array (line,columns).

        """

        if (wave  > self.wavelengths.max()) | (wave  < self.wavelengths.min()):
            print("Input wavelength outside wavelength range!")
            band = None
        else:
            band_num = np.argmin(np.abs(self.wavelengths - wave))
            band = self.get_band(band_num,mask_values)
        return band

    def get_pixels(self,lines,columns):
        """
        Args:
            lines (list): List of zero-indexed line indices.
            columns (list): List of zero-indexed column indices.

        Returns:
            numpy.ndarray: Pixel array (pixels,bands).

        """

        self.load_data()
        if self.file_type == "neon":
            pixels = []
            for line,column in zip(lines,columns):
                pixels.append(self.data[line,column,:])
            pixels = np.array(pixels)
        elif self.file_type == "envi":
            pixels = envi_read_pixels(self.data,lines,columns,self.interleave)
        self.close_data()
        return pixels

    def get_line(self,index):
        """
        Args:
            index (int): Zero-indexed line index.

        Returns:
            numpy.ndarray: Line array (columns, bands).

        """

        self.load_data()
        if self.file_type == "neon":
            line = self.data[index,:,:]
        elif self.file_type == "envi":
            line = envi_read_line(self.data,index,self.interleave)
        self.close_data()
        return line

    def get_column(self,index):
        """
        Args:
            index (int): Zero-indexed column index.

        Returns:
            numpy.ndarray: Column array (lines, bands).

        """

        self.load_data()
        if self.file_type == "neon":
            column = self.data[:,index,:]
        elif self.file_type == "envi":
            column = envi_read_column(self.data,index,self.interleave)
        self.close_data()
        return column

    def get_chunk(self,col_start,col_end,line_start,line_end):
        """
        Args:
            col_start (int): Chunk starting column.
            col_end (int): Noninclusive chunk ending column index.
            line_start (int): Chunk starting line.
            line_end (int): Noninclusive chunk ending line index.

        Returns:
            numpy.ndarray: Chunk array (line_end-line_start,col_end-col_start,bands).

        """

        self.load_data()
        if self.file_type == "neon":
            chunk = self.data[line_start:line_end,col_start:col_end,:]
        elif self.file_type == "envi":
            chunk =  envi_read_chunk(self.data,col_start,col_end,
                                     line_start,line_end,self.interleave)
        self.close_data()
        return chunk


    def get_obs(self,obs, radians = True):
        """ Read observable dataset to memory.

        Args:
            obs (str): Observable name.
            radians (bool, optional): Convert angular measures to radians. Defaults to True.

        Returns:
            obs_data (numpy.ndarray)

        """

        angular_obs = ['slope','sensor_az','sensor_zn','aspect','solar_zn','solar_az']

        if self.file_type == "envi":
            observables = open_envi(self.observables[obs][0])
            observables.load_data()
            obs_data = np.copy(observables.get_band(self.observables[obs][1]))
            observables.close_data()

        else:
            hdf_obj = h5py.File(self.file_name,'r')
            metadata = hdf_obj[self.base_key]["Reflectance"]["Metadata"]
            keys = self.observables[obs]
            for key in keys:
                metadata = metadata[key]
            obs_data = metadata[()]

            #Make solar geometry into 2D array
            if obs in ['solar_zn','solar_az']:
                obs_data = np.ones((self.lines, self.columns)) * obs_data
            hdf_obj.close()

        if radians and (obs in angular_obs):
            obs_data= np.radians(obs_data)

        return obs_data


    def volume_kernel(self,kernel):
        """Calculate volume scattering kernel.
        """

        return calc_volume_kernel(self.get_obs('solar_az'), self.get_obs('solar_zn'),
                                  self.get_obs('sensor_az'), self.get_obs('sensor_zn'),
                                               kernel)


    def geom_kernel(self,kernel,b_r=10.,h_b =2.):
        """Calculate volume scattering kernel.
        """

        return calc_geom_kernel(self.get_obs('solar_az'),self.get_obs('solar_zn'),
                                self.get_obs('sensor_az'),self.get_obs('sensor_zn'),
                                kernel,b_r=b_r,h_b =h_b)


    def set_mask(self,mask):
        """Generate mask using masking function which takes a HyTools object as
        an argument.
        """
        self.mask = mask

    def gen_mask(self,masker):
        """Generate mask using masking function which takes a HyTools object as
        an argument.
        """
        self.mask = masker(self)

    def do(self,function):
        """Run a function and return the results
        """
        return function(self)


class Iterator:
    """Iterator class
    """

    def __init__(self,hy_obj,by,chunk_size = None):
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
        self.hy_obj.load_data()

    def read_next(self):
        """ Return next line/column/band/chunk.
        """

        if self.by == "line":
            self.current_line +=1
            if self.current_line == self.hy_obj.lines-1:
                self.complete = True
                subset = None
            if self.hy_obj.file_type == "neon":
                subset =  self.hy_obj.data[self.current_band,:,:]
            else:
                subset =  envi_read_line(self.hy_obj.data,self.current_line,
                                         self.hy_obj.interleave)

        elif self.by == "column":
            self.current_column +=1
            if self.current_column == self.hy_obj.columns-1:
                self.complete = True
            if self.hy_obj.file_type == "neon":
                subset =  self.hy_obj.data[:,self.current_band,:]
            else:
                subset =  envi_read_column(self.hy_obj.data,self.current_column,
                                           self.hy_obj.interleave)

        elif self.by == "band":
            self.current_band +=1
            if self.current_band == self.hy_obj.bands-1:
                self.complete = True
            if self.hy_obj.file_type == "neon":
                subset =  self.hy_obj.data[:,:,self.current_band]
            else:
                subset =  envi_read_band(self.hy_obj.data,self.current_band,
                                         self.hy_obj.interleave)

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

            if self.hy_obj.file_type == "neon":
                subset = self.hy_obj.data[y_start:y_end,x_start:x_end,:]
            else:
                subset =  envi_read_chunk(self.hy_obj.data,x_start,x_end,
                                          y_start,y_end,self.hy_obj.interleave)
            if (y_end == self.hy_obj.lines) and (x_end == self.hy_obj.columns):
                self.complete = True

        if self.complete:
            self.hy_obj.close_data()

        return subset

    def reset(self):
        """Reset counters.
        """
        self.current_column = -1
        self.current_line = -1
        self.current_band = -1
        self.complete = False




