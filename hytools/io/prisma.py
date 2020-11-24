"""
 Functions for reading PRISMA formatted HDF files
"""
import h5py,os
import numpy as np


class HyToolsPRISMA(object):
    """HyTools PRISMA class container"""
    
    def __init__(self):
              
        self.file_type = "PRISMA"
        self.interleave = np.nan
        self.file_name = np.nan
        self.bad_bands = []
        self.no_data = np.nan
        self.map_info = np.nan
        self.crs = np.nan
        self.ulX = np.nan
        self.ulY = np.nan
        self.dtype = np.nan
        self.data = np.nan
        self.header_dict = np.nan
    
    def load_data(self):
        """Load data object to memory.
        
        Parameters
        ----------
        mode: str 
            File read mode, default: read-only
            
        """
        self.hdfObj = h5py.File(self.file_name,'r')
        base_key = [key for key in self.hdfObj ['HDFEOS']["SWATHS"].keys() if 'HCO' in key][0]  
        self.vnir_data = self.hdfObj ['HDFEOS']["SWATHS"][base_key]['Data Fields']['VNIR_Cube'] 
        self.swir_data = self.hdfObj ['HDFEOS']["SWATHS"][base_key]['Data Fields']['SWIR_Cube'] 
                    
    def close_data(self):
        """Close data object.
        """
        self.hdfObj.close()


def openPRISMA(srcFile):
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

    # Load metadata and populate HyTools object
    hdfObj = h5py.File(srcFile,'r')
    
    base_key = [key for key in hdfObj['HDFEOS']["SWATHS"].keys() if 'HCO' in key][0]  
    vnir_data = hdfObj['HDFEOS']["SWATHS"][base_key]['Data Fields']['VNIR_Cube'] 
    swir_data = hdfObj['HDFEOS']["SWATHS"][base_key]['Data Fields']['SWIR_Cube'] 
    
    resolution =30
    projection =hdfObj.attrs.get('Projection_Name').decode('UTF-8')
    zone =hdfObj.attrs.get('Projection_Id').decode('UTF-8')
    direction = 'N'
    if hdfObj.attrs.get('Product_center_lat') < 0:
        direction = 'S'#May not be accurate in scenes that straddle the equator....
        
    corner_easting = hdfObj.attrs.get('Product_ULcorner_easting')
    corner_northing = hdfObj.attrs.get('Product_ULcorner_northing')    
    
    map_info_string = [projection, 1, 1,corner_easting ,
                       corner_northing,resolution,resolution,
                       zone, direction, 'WGS-84','units=Meters']
    
    hyObj.map_info = [projection, 1, 1,corner_easting ,corner_northing,resolution,resolution,zone, direction, 'WGS-84','units=Meters']
    hyObj.transform = (float(hyObj.map_info [3]),float(hyObj.map_info [1]),0,float(hyObj.map_info [4]),0,-float(hyObj.map_info [2]))
    
    hyObj.fwhm_vnir =   hdfObj.attrs.get('List_Fwhm_Vnir')
    hyObj.wavelengths_vnir = hdfObj.attrs.get('List_Fwhm_Vnir')

    hyObj.fwhm_swir =   hdfObj.attrs.get('List_Fwhm_Swir')
    hyObj.wavelengths_swir = hdfObj.attrs.get('List_Fwhm_Swir')

    hyObj.wavelength_units = "nanometers"
  
    hyObj.lines_vnir = vnir_data.shape[0]
    hyObj.column_vnir = vnir_data.shape[1]
    hyObj.bands_vnir = vnir_data.shape[2]
    
    hyObj.lines_swir = swir_data.shape[0]
    hyObj.column_swir = swir_data.shape[1]
    hyObj.bands_swir = swir_data.shape[2]

    hyObj.file_name = srcFile

    hdfObj.close()
    return hyObj

















