"""
 Functions for reading NEON formatted HDF files
"""
import h5py,os
import numpy as np
from .base import HyTools

class HyToolsNEON(object):
    """HyTools NEON class container"""
    
    def __init__(self):
              
        self.file_type = "NEON"
        self.interleave = np.nan
        self.file_name = np.nan
        self.shape = np.nan
        self.lines = np.nan
        self.columns = np.nan
        self.bands = np.nan
        self.wavelengths = np.nan
        self.fwhm = []
        self.bad_bands = []
        self.no_data = np.nan
        self.map_info = np.nan
        self.crs = np.nan
        self.ulX = np.nan
        self.ulY = np.nan
        self.dtype = np.nan
        self.data = np.nan
        self.header_dict = np.nan
        self.solar_zn = []
        self.solar_az = []
        self.sensor_zn = []
        self.sensor_az = []
        
    def create_bad_bands(self,bad_regions):
        """Create bad bands list based upon spectrum regions. Good: 1, bad : 0.
        
        Parameters
        ----------
        bad_regions : list
            List of lists containing start and end values of wavelength regions considered bad.
            ex: [[350,400].....[2450,2500]] Wavelengths should be in the same units as
            data units. Assumes wavelength attribute is populated.
        """
        bad_bands = []
        
        for wavelength in self.wavelengths:
            bad =True
            for start,end in bad_regions:
                if (wavelength >= start) and (wavelength <=end):
                    bad = False
            bad_bands.append(bad)             
        self.bad_bands = np.array(bad_bands)
    
    def load_data(self):
        """Load data object to memory.
        
        Parameters
        ----------
        mode: str 
            File read mode, default: read-only
            
        """
        self.hdfObj = h5py.File(self.file_name,'r')
        base_key = list(self.hdfObj.keys())[0]
        self.data = self.hdfObj[base_key]["Reflectance"]["Reflectance_Data"] 
                
    def close_data(self):
        """Close data object.
        """
        self.hdfObj.close()
         
    def iterate(self,by,chunk_size= (100,100)):    
        """Return data iterator.
        
        Parameters
        ----------
        by: str
            Dimension along which to iterate. 
            Lines,columns,bands or chunks.
        chunk_size : shape (columns , rows)
            Size of chunks to iterate over, only applicable when
            by == chunks.
        
        Returns
        -------
        iterator: hyTools iterator class
            
        """
        iterator = iterHDF(self.data,by,chunk_size)

        return iterator     

    
    def get_band(self,band):
        """Return the i-th band of the image.

        Parameters
        ----------
        band: int
                Zero-indexed band index
        Returns
        -------
        band : np.array (lines, columns)
        """
        
        band = hdf_read_band(self.data,band)
        return band
    
    
    def get_wave(self,wave):
        """Return the band image corresponding to the input wavelength, 
        if not an exact match the closest wavelength will be returned.

        Parameters
        ----------
        wave: int
                Wavelength of band to be gotten.
        Returns
        -------
        band : np.array (lines, columns)
        """
        
        # Perform wavelength unit conversion if nescessary
        if self.wavelength_units == "micrometers" and wave > 3:
            wave/= 1000
        if self.wavelength_units == "nanometers" and wave < 3:
            wave*= 1000

        if wave in self.wavelengths:
            band = np.argwhere(wave == self.wavelengths)[0][0]
        elif (wave  > self.wavelengths.max()) | (wave  < self.wavelengths.min()):
            print("Input wavelength outside image range!")
            return
        else: 
            band = np.argmin(np.abs(self.wavelengths - wave))
        
        band = hdf_read_band(self.data,band)

        return band
    
            
    def wave_to_band(self,wave):
        """Return band number corresponding to input wavelength. Return closest band if
           not an exact match. 
          
           wave : float/int
                  Wavelength 

        """
        # Perform wavelength unit conversion if nescessary
        if self.wavelength_units == "micrometers" and wave > 3:
            wave/= 1000
        if self.wavelength_units == "nanometers" and wave < 3:
            wave*= 1000

        if wave in self.wavelengths:
            band = np.argwhere(wave == self.wavelengths)[0][0]
        elif (wave  > self.wavelengths.max()) | (wave  < self.wavelengths.min()):
            print("Input wavelength outside image range!")
            band = np.nan
        else: 
            band = np.argmin(np.abs(self.wavelengths - wave))
        return band
          
    def get_line(self,line):        
        """Return the i-th band of the image.
        
        Parameters
        ----------
        band: int
                Zero-indexed band index

        Returns
        -------
        line : np.array (columns, bands)
        """
        
        line = hdf_read_line(self.data,line)
        return line
            
    def get_column(self,column):
        """Return the i-th column of the image.
       
        Parameters
        ----------
        column: int
                Zero-indexed column index

        Returns
        -------
        column : np.array (lines, bands)
        """
        column = hdf_read_column(self.data,column)
        return column
        
    def get_chunk(self,x_start,x_end,y_start,y_end):
        """Return chunk from image.

        Parameters
        ----------
        x_start : int
        x_end : int
        y_start : int
        y_end : int
            
        Returns
        -------
        chunk : np.array (y_end-y_start,x_end-x_start, bands)
        """
        chunk = hdf_read_chunk(self.data,x_start,x_end,y_start,y_end)
        return chunk


def openNEON(srcFile, no_data = -9999,load_obs = False):
    """Load and parse NEON formated HDF image into a HyTools data object
        
    Parameters
    ----------
    srcFile : str
        pathname of input HDF file

    no_data: int
        No data value
    """

    if not os.path.isfile(srcFile):
        print("File not found.")
        return
    
    hyObj = HyToolsNEON()

    # Load metadata and populate HyTools object
    hdfObj = h5py.File(srcFile,'r')
    base_key = list(hdfObj.keys())[0]
    metadata = hdfObj[base_key]["Reflectance"]["Metadata"]
    data = hdfObj[base_key]["Reflectance"]["Reflectance_Data"] 
    hyObj.projection = metadata['Coordinate_System']['Coordinate_System_String'][()].decode("utf-8")
    hyObj.map_info = metadata['Coordinate_System']['Map_Info'][()].decode("utf-8").split(',')
    hyObj.transform = (float(hyObj.map_info [3]),float(hyObj.map_info [1]),0,float(hyObj.map_info [4]),0,-float(hyObj.map_info [2]))
    hyObj.fwhm =  metadata['Spectral_Data']['FWHM'][()]
    hyObj.wavelengths = metadata['Spectral_Data']['Wavelength'][()]
    
    #If wavelengths units are not specified guess
    try:
        hyObj.wavelength_units = metadata['Spectral_Data']['Wavelength'].attrs['Units']
    except:
        if hyObj.wavelengths.min() >100:
            hyObj.wavelength_units = "nanometers"
        else:
            hyObj.wavelength_units = "micrometers"
            
            
    hyObj.lines = data.shape[0]
    hyObj.columns = data.shape[1]
    hyObj.bands = data.shape[2]
    hyObj.no_data = no_data
    hyObj.file_name = srcFile
 
    # Map observables to memory
    if load_obs: 
        hyObj.solar_zn = np.ones((hyObj.lines, hyObj.columns)) * np.radians(metadata['Logs']['Solar_Zenith_Angle'][()])
        hyObj.solar_az = np.ones((hyObj.lines, hyObj.columns)) * np.radians(metadata['Logs']['Solar_Azimuth_Angle'][()])
        hyObj.sensor_zn = np.radians(metadata['to-sensor_Zenith_Angle'][:,:])
        hyObj.sensor_az = np.radians(metadata['to-sensor_Azimuth_Angle'][:,:])
        hyObj.slope = np.radians(metadata['Ancillary_Imagery']['Slope'][()])
        hyObj.aspect =  np.radians(metadata['Ancillary_Imagery']['Aspect'][()])
        
    hdfObj.close()
    return hyObj


class iterHDF(object):
    """Iterator class for reading HDF data file.
    
    """

    def __init__(self,data,by, chunk_size = None):
        """
        
        Parameters
        ----------
        data : memmap object
    
        by: iterator slice lines, columns, bands or chunks
        
        chunk_size: y,x chunks size
            
            
        """
        self.chunk_size= chunk_size
        self.by = by
        self.current_column = -1
        self.current_line = -1
        self.current_band = -1
        self.data = data
        self.complete = False
        self.lines,self.columns,self.bands = data.shape

    def read_next(self):
        """ Return next line/column/band/chunk.
        
        """
        if self.by == "line":
            self.current_line +=1
            if self.current_line == self.lines-1:
                self.complete = True
                subset = np.nan
            subset =  hdf_read_line(self.data,self.current_line)

        elif self.by == "column":
            self.current_column +=1
            if self.current_column == self.columns-1:
                self.complete = True
            subset =  hdf_read_column(self.data,self.current_column)

        elif self.by == "band":
            self.current_band +=1
            if self.current_band == self.bands-1:
                self.complete = True
            subset =  hdf_read_band(self.data,self.current_band)

        elif self.by == "chunk":
            
            if self.current_column == -1:
                self.current_column +=1
                self.current_line +=1
            else:
                self.current_column += self.chunk_size[1]
                
            if self.current_column >= self.columns:
                self.current_column = 0
                self.current_line += self.chunk_size[0]

            # Set array indices for current chunk and update current line and column.
            y_start = self.current_line
            y_end = self.current_line + self.chunk_size[0]  
            if y_end >= self.lines:
                y_end = self.lines 
            x_start = self.current_column 
            x_end = self.current_column + self.chunk_size[1]
            if x_end >= self.columns:
                x_end = self.columns 

            subset =  hdf_read_chunk(self.data,x_start,x_end,y_start,y_end)
            if (y_end == self.lines) and (x_end == self.columns):
                self.complete = True
         
        return subset
        
    def reset(self):
        """Reset counters.
        """
        self.current_column = -1
        self.current_line = -1
        self.current_band = -1
        self.complete = False
        
        
def hdf_read_line(dataArray,line):
    """ Read line from hdf file.
    
    """
    lineSubset = dataArray[line,:,:] 
    return lineSubset

def hdf_read_column(dataArray,column):
    """ Read column from hdf file.
    
    """
    columnSubset = dataArray[:,column,:] 
    return columnSubset 
    
def hdf_read_band(dataArray,band):
    """ Read band from hdf file.
    
    """
    bandSubset =  dataArray[:,:,band] 
    return bandSubset

    
def hdf_read_chunk(dataArray,x_start,x_end,y_start,y_end):
    """ Read chunk from hdf file.
    """
    chunkSubset = dataArray[y_start:y_end,x_start:x_end,:]
    return chunkSubset