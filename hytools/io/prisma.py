"""
 Functions for reading PRISMA formatted HDF files
"""
import h5py,os
import numpy as np


'''
UNDER CONSTRUCTION

class HyToolsPRISMA(object):
    """HyTools PRISMA class container"""
    
    def __init__(self):
              
        self.file_type = "PRISMA"
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
        self.mask = np.nan
        
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
    
    def load_data(self, mode = 'r', offset = 0):
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


        
    def set_mask(self,mask):
        """Set mask for image analysis.
          
          mask: m x n numpy array 
               A boolean mask to exclude pixels from analysis, shape should be the same
               as the number of line and columns of the image.

        """
        
        if mask.shape == (self.lines,self.columns):
            self.mask = mask
        else:
            print("Error: Shape of mask does not match shape of image.")



def openPRISMA(srcFile, no_data = -9999,load_obs = False):
    """Load and parse PRISMA formated HDF image into a HyTools data object
        
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
    
    hyObj = HyToolsPRISMA()

    return hyObj

"""