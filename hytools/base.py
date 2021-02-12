"""
Base

TODO: Add corrections to ndi()

"""
import os
import json
import numpy as np
import h5py
import warnings
from .io.envi import envi_read_band,envi_read_pixels
from .io.envi import envi_read_line,envi_read_column,envi_read_chunk
from .io.envi import open_envi,parse_envi_header,envi_header_from_neon
from .io.neon import open_neon
from .brdf import apply_brdf_correct
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
        self.byte_order = None
        self.columns = None
        self.corrections = []
        self.crs = None
        self.data = None
        self.dtype = None
        self.file_name = None
        self.file_type = None
        self.fwhm = []
        self.hdf_obj  = None
        self.interleave = None
        self.lines = None
        self.map_info = None
        self.mask = {}
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

    def read_file(self,file_name,file_type,anc_path = None):
        self.file_name = file_name
        self.file_type = file_type

        if file_type == 'envi':
            open_envi(self,anc_path)
        elif file_type == "neon":
            open_neon(self)
        else:
            print("Unrecognized file type.")

        # Create a no data mask
        self.mask['no_data'] = self.get_band(0) != self.no_data
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

    def close_data(self):
        """Close data object.

        """
        if self.file_type  == "envi":
            del self.data
        elif self.file_type  == "neon":
            self.hdf_obj.close()
            self.hdf_obj = None
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
        elif self.file_type == "envi":
            band = envi_read_band(self.data,index,self.interleave)
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
        if self.file_type == "neon":
            pixels = []
            for line,column in zip(lines,columns):
                pixels.append(self.data[line,column,:])
            pixels = np.array(pixels)
        elif self.file_type == "envi":
            pixels = envi_read_pixels(self.data,lines,columns,self.interleave)
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
        if self.file_type == "neon":
            line = self.data[index,:,:]
        elif self.file_type == "envi":
            line = envi_read_line(self.data,index,self.interleave)
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
        if self.file_type == "neon":
            column = self.data[:,index,:]
        elif self.file_type == "envi":
            column = envi_read_column(self.data,index,self.interleave)
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
        if self.file_type == "neon":
            chunk = self.data[line_start:line_end,col_start:col_end,:]
        elif self.file_type == "envi":
            chunk =  envi_read_chunk(self.data,col_start,col_end,
                                     line_start,line_end,self.interleave)
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
                data =  apply_topo_correct(self,data,dimension,index)
            elif correction == 'brdf':
                data =  apply_brdf_correct(self,data,dimension,index)
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
            ancillary.close_data()

        else:
            hdf_obj = h5py.File(self.file_name,'r')
            metadata = hdf_obj[self.base_key]["Reflectance"]["Metadata"]
            keys = self.anc_path[anc]
            for key in keys:
                metadata = metadata[key]
            anc_data = metadata[()]

            #Make solar geometry into 2D array
            if anc in ['solar_zn','solar_az']:
                anc_data = np.ones((self.lines, self.columns)) * anc_data
            hdf_obj.close()

        if radians and (anc in angular_anc):
            anc_data= np.radians(anc_data)


        if mask:
            anc_data = anc_data[self.mask[mask]]

        return anc_data

    def load_anc(self,anc,radians = True):
        self.ancillary[anc] = self.get_anc(self,anc,radians)

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


    def get_header(self):
        """ Return header dictionary

        """
        if self.file_type == "neon":
            header_dict = envi_header_from_neon(self)
        elif self.file_type == "envi":
            header_file = os.path.splitext(self.file_name)[0] + ".hdr"
            header_dict = parse_envi_header(header_file)
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