"""
 Functions for reading ENVI formatted binary files

See https://www.l3harrisgeospatial.com/docs/enviheaderfiles.html for formatting details
"""
import numpy as np
import os

# ENVI datatype conversion dictionary
dtypeDict = {1:np.uint8,
             2:np.int16,
             3:np.int32,
             4:np.float32,
             5:np.float64,
             12:np.uint16,
             13:np.uint32,
             14:np.int64,
             15:np.uint64}

class HyToolsENVI(object):
    """HyTools  class object"""
    
    def __init__(self):
              
        self.file_type = "ENVI"
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
        
        self.data = np.memmap(self.file_name,dtype = self.dtype, mode=mode, shape = self.shape,offset=offset)

# Needs testing            
#            if self.byte_order == 1:
#                self.data = self.data.byteswap()
#                self.byte_order = 0
#                self.header_dict['byte order'] = 0

                
    def close_data(self):
        """Close data object.
        """

        del self.data

         
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

        iterator = iterENVI(self.data,by,self.interleave,chunk_size)
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
        
        band = envi_read_band(self.data,band,self.interleave)
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
        
        band = envi_read_band(self.data,band,self.interleave)
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

        line = envi_read_line(self.data,line,self.interleave)
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

        column = envi_read_column(self.data,column,self.interleave)
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

        chunk =   envi_read_chunk(self.data,x_start,x_end,y_start,y_end,self.interleave)
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

    
    # Move to ENVI...
    # def load_obs(self,observables):
    #     """
    #     Load observables to memory.
        
    #     """
    #     if self.file_type == "ENVI":
    #         observables = openENVI(observables)
    #         observables.load_data()
    #         self.sensor_az = np.radians(observables.get_band(1))
    #         self.sensor_zn = np.radians(observables.get_band(2))
    #         self.solar_az = np.radians(observables.get_band(3))
    #         self.solar_zn = np.radians(observables.get_band(4))
    #         self.slope = np.radians(observables.get_band(6))
    #         self.aspect = np.radians(observables.get_band(7))
    #         observables.close_data()
                
                
                

def openENVI(srcFile):
    """Load and parse ENVI image header into a HyTools data object
    
    Parameters
    ----------
    srcFile : str
        Pathname of input ENVI data file, header assumed to be located in 
        same directory
        
    Returns
    -------
    Populated HyTools data object

    """
    
    if not os.path.isfile(os.path.splitext(srcFile)[0] + ".hdr"):
        print("ERROR: Header file not found.")
        return None

    hyObj = HyToolsENVI()

    # Load header into dictionary
    header_dict = parse_ENVI_header(os.path.splitext(srcFile)[0] + ".hdr")

    # Assign HyTools object attributes
    hyObj.lines =  header_dict["lines"]
    hyObj.columns =  header_dict["samples"]
    hyObj.bands =   header_dict["bands"]
    hyObj.interleave =  header_dict["interleave"]
    hyObj.fwhm =  header_dict["fwhm"]
    hyObj.wavelengths = header_dict["wavelength"]
    hyObj.wavelength_units = header_dict["wavelength units"]
    hyObj.dtype = dtypeDict[header_dict["data type"]]
    hyObj.no_data = header_dict['data ignore value']
    hyObj.map_info = header_dict['map info']
    hyObj.byte_order = header_dict['byte order']
    hyObj.header_dict =  header_dict 

    hyObj.file_name = srcFile
    
    if type(header_dict['bbl']) == np.ndarray:
        hyObj.bad_bands = np.array([x==1 for x in header_dict['bbl']])

    if header_dict["interleave"] == 'bip':    
        hyObj.shape = (hyObj.lines, hyObj.columns, hyObj.bands)
    elif header_dict["interleave"] == 'bil':    
        hyObj.shape = (hyObj.lines, hyObj.bands, hyObj.columns) 
    elif header_dict["interleave"] == 'bsq':
        hyObj.shape = (hyObj.bands, hyObj.lines, hyObj.columns)
    else:
        print("ERROR: Unrecognized interleave type.")
        hyObj = None
        
    #Convert all units to nanometers
    if hyObj.wavelength_units == "micrometers":
       hyObj.wavelength_units ="nanometers" 
       hyObj.wavelengths*=10**3
       hyObj.fwhm*= 10**3
       
    if hyObj.wavelength_units not in ["nanometers","micrometers"]:
        try:
            if hyObj.wavelengths.min() <100:
                hyObj.wavelengths*=10**3
                hyObj.fwhm*= 10**3        
            hyObj.wavelength_units = "nanometers"
        except:
            hyObj.wavelength_units = "unknown"
            
    # If no_data value is not specified guess using image corners.   
    if np.isnan(hyObj.no_data):  
        print("No data value specified, guessing.")
        hyObj.load_data()
        ul = hyObj.data[0,0,0]
        ur = hyObj.data[0,-1,0]
        ll = hyObj.data[-1,0,0]
        lr = hyObj.data[-1,-1,0]
        counts = {v: k for k, v in Counter([ul,ur,ll,lr]).items()}
        hyObj.no_data = counts[max(counts.keys())]
        hyObj.close_data()
        
  
    del header_dict
    return hyObj    


class writeENVI(object):
    """Iterator class for writing to an ENVI data file.
    
    """
    
    
    def __init__(self,output_name,headerDict):
        """
        Parameters
        ----------
        srcFile : str
            Pathname of output ENVI data file
        
        headerDict : dict
            Dictionary containing ENVI header information
        Returns
        -------
        Populated hyTools data object
            
        """
    
        self.interleave = headerDict['interleave']
        self.headerDict = headerDict
        self.output_name =output_name
        dtype = dtypeDict[headerDict["data type"]]
        lines = headerDict['lines']
        columns = headerDict['samples']
        bands = headerDict['bands']
        
        # Create numpy mem map file on disk
        if self.interleave == "bip":
            self.data = np.memmap(output_name,dtype = dtype, mode='w+', shape = (lines,columns,bands))
        elif self.interleave == "bil":
            self.data = np.memmap(output_name,dtype = dtype, mode='w+', shape =(lines,bands,columns))
        elif self.interleave == "bsq":
            self.data = np.memmap(output_name,dtype = dtype, mode='w+',shape =(bands,lines,columns))    
        write_ENVI_header(self.output_name,self.headerDict)    
            
    def write_line(self,dataArray,line):
        """ Write line to ENVI file.
        
        """
               
        if self.interleave == "bip":
            self.data[line,:,:] = dataArray
    
        elif self.interleave == "bil":
            self.data[line,:,:] = dataArray
        
        elif self.interleave == "bsq":
            self.data[:,line,:] = dataArray
                
    def write_column(self,dataArray,column):
        """ Write column to ENVI file.
        
        """
           
        if self.interleave == "bip":
            self.data[:,column,:]  = dataArray
        elif self.interleave == "bil":
            self.data[:,:,column] = dataArray     
        elif self.interleave == "bsq":
            self.data[:,:,column] = dataArray
                    
    def write_band(self,dataArray,band):
        """ Write band to ENVI file.
        
        """
        if self.interleave == "bip":
            self.data[:,:,band]  = dataArray
        elif self.interleave == "bil":
            self.data[:,band,:] = dataArray
        elif self.interleave == "bsq":
            self.data[band,:,:]= dataArray
            
    
        
    def write_chunk(self,dataArray,line,column):
        """ Write chunk to ENVI file.
        
        """
    
        x_start = column 
        x_end = column + dataArray.shape[1]
        y_start = line
        y_end = line + dataArray.shape[0]
    
        if self.interleave == "bip":
            self.data[y_start:y_end,x_start:x_end,:] = dataArray
        elif self.interleave == "bil":
            self.data[y_start:y_end,:,x_start:x_end] = np.moveaxis(dataArray,-1,-2)
        elif self.interleave == "bsq":
            self.data[:,y_start:y_end,x_start:x_end] = np.moveaxis(dataArray,-1,0)
        
    def close(self):
        del self.data
                        
        

def ENVI_header_from_hdf(hyObj, interleave = 'bil'):
    """Create an ENVI header dictionary from HDF metadata
    """

    headerDict = {}

    headerDict["ENVI description"] = "{}"
    headerDict["samples"] = hyObj.columns
    headerDict["lines"]   = hyObj.lines
    headerDict["bands"]   = hyObj.bands
    headerDict["header offset"] = 0
    headerDict["file type"] = "ENVI Standard"
    headerDict["data type"] = 2
    headerDict["interleave"] = interleave
    headerDict["sensor type"] = ""
    headerDict["byte order"] = 0
    headerDict["map info"] = hyObj.map_info
    headerDict["coordinate system string"] = hyObj.projection
    headerDict["wavelength units"] = hyObj.wavelength_units
    headerDict["data ignore value"] =hyObj.no_data
    headerDict["wavelength"] =hyObj.wavelengths

    return headerDict
    

def write_ENVI_header(output_name,headerDict):
    """Parse ENVI header into dictionary
    """

    headerFile = open(output_name + ".hdr",'w+')
    headerFile.write("ENVI\n")
    
    for key in headerDict.keys():
        value = headerDict[key]
        # Convert list to comma seperated strings
        if type(value) == list or type(value) == np.ndarray:
            value = "{%s}" % ",".join(map(str, value))
        else:
            value = str(value)
        
        # Skip entires with nan as value
        if value != 'nan':
            headerFile.write("%s = %s\n" % (key,value))
    
    headerFile.close()
 


def empty_ENVI_header_dict():
    # Dictionary of all types
    headerDict = {"acquisition time": np.nan,
                 "band names":np.nan, 
                 "bands": np.nan, 
                 "bbl": np.nan,
                 "byte order": np.nan,
                 "class lookup": np.nan,
                 "class names": np.nan,
                 "classes": np.nan,
                 "cloud cover": np.nan,
                 "complex function": np.nan,
                 "coordinate system string": np.nan,
                 "correction factors": np.nan,
                 "data gain values": np.nan,
                 "data ignore value": np.nan,
                 "data offset values": np.nan,
                 "data reflectance gain values": np.nan,
                 "data reflectance offset values": np.nan,
                 "data type": np.nan,
                 "default bands": np.nan,
                 "default stretch": np.nan,
                 "dem band": np.nan,
                 "dem file": np.nan,
                 "description": np.nan,
                 "envi description":np.nan,
                 "file type": np.nan,
                 "fwhm": np.nan,
                 "geo points": np.nan,
                 "header offset": np.nan,
                 "interleave": np.nan,
                 "lines": np.nan,
                 "map info": np.nan,
                 "pixel size": np.nan,
                 "projection info": np.nan,
                 "read procedures": np.nan,
                 "reflectance scale factor": np.nan,
                 "rpc info": np.nan,
                 "samples":np.nan,
                 "security tag": np.nan,
                 "sensor type": np.nan,
                 "smoothing factors": np.nan,
                 "solar irradiance": np.nan,
                 "spectra names": np.nan,
                 "sun azimuth": np.nan,
                 "sun elevation": np.nan,
                 "wavelength": np.nan,
                 "wavelength units": np.nan,
                 "x start": np.nan,
                 "y start": np.nan,
                 "z plot average": np.nan,
                 "z plot range": np.nan,
                 "z plot titles": np.nan}
    return headerDict



class iterENVI(object):
    """Iterator class for reading ENVI data file.
    
    """

    def __init__(self,data,by,interleave, chunk_size = None):
        """
        
        Parameters
        ----------
        data : memmap object
    
        by: iterator slice lines, columns, bands or chunks
        
        chunk_size: y,x chunks size
            
            
        """
        self.interleave = interleave
        self.chunk_size= chunk_size
        self.by = by
        self.current_column = -1
        self.current_line = -1
        self.current_band = -1
        self.data = data
        self.complete = False
        
        if interleave == "bip":
            self.lines,self.columns,self.bands = data.shape
        elif interleave == "bil":
            self.lines,self.bands,self.columns = data.shape
        elif interleave == "bsq":
            self.bands,self.lines,self.columns = data.shape
        else:
            print("ERROR: Iterator unit not recognized.")   
    
    def read_next(self):
        """ Return next line/column/band/chunk.
        
        """
        if self.by == "line":
            self.current_line +=1
            if self.current_line == self.lines-1:
                self.complete = True
                subset = np.nan
            subset =  envi_read_line(self.data,self.current_line,self.interleave)

        elif self.by == "column":
            self.current_column +=1
            if self.current_column == self.columns-1:
                self.complete = True
            subset =  envi_read_column(self.data,self.current_column,self.interleave)

        elif self.by == "band":
            self.current_band +=1
            if self.current_band == self.bands-1:
                self.complete = True
            subset =  envi_read_band(self.data,self.current_band,self.interleave)

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

            subset =  envi_read_chunk(self.data,x_start,x_end,y_start,y_end,self.interleave)
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
        
        
def envi_read_line(dataArray,line,interleave):
    """ Read line from ENVI file.
    
    """
           
    if interleave == "bip":
        lineSubset = dataArray[line,:,:] 

    elif interleave == "bil":
        lineSubset = dataArray[line,:,:] 
    
    elif interleave == "bsq":
        lineSubset = dataArray[:,line,:]
        
    return lineSubset

def envi_read_column(dataArray,column,interleave):
    """ Read column from ENVI file.
    
    """
       
    if interleave == "bip":
        columnSubset = dataArray[:,column,:] 
    elif interleave == "bil":
        columnSubset = dataArray[:,:,column]     
    elif interleave == "bsq":
       columnSubset =  dataArray[:,:,column]
        
    return columnSubset 
    
def envi_read_band(dataArray,band,interleave):
    """ Read band from ENVI file.
    
    """
       
    if interleave == "bip":
        bandSubset =  dataArray[:,:,band] 
    elif interleave == "bil":
        bandSubset = dataArray[:,band,:] 
    elif interleave == "bsq":
        bandSubset = dataArray[band,:,:]
        
    return bandSubset

    
def envi_read_chunk(dataArray,x_start,x_end,y_start,y_end,interleave):
    """ Read chunk from ENVI file.
    """

    if interleave == "bip":
        chunkSubset = dataArray[y_start:y_end,x_start:x_end,:]
    elif interleave == "bil":
        chunkSubset = np.moveaxis(dataArray[y_start:y_end,:,x_start:x_end],-1,-2)
    elif interleave == "bsq":
        chunkSubset = np.moveaxis(dataArray[:,y_start:y_end,x_start:x_end],0,-1)
    
    return chunkSubset


     
def parse_ENVI_header(hdrFile):
    """Parse ENVI header into dictionary
    """

    # Dictionary of all types
    fieldDict = {"acquisition time": "str",
                 "band names":"list_str", 
                 "bands": "int", 
                 "bbl": "list_float",
                 "byte order": "int",
                 "class lookup": "str",
                 "class names": "str",
                 "classes": "int",
                 "cloud cover": "float",
                 "complex function": "str",
                 "coordinate system string": "str",
                 "correction factors": "list_float",
                 "data gain values": "list_float",
                 "data ignore value": "float",
                 "data offset values": "list_float",
                 "data reflectance gain values": "list_float",
                 "data reflectance offset values": "list_float",
                 "data type": "int",
                 "default bands": "list_float",
                 "default stretch": "str",
                 "dem band": "int",
                 "dem file": "str",
                 "description": "str",
                 "envi description":"str",
                 "file type": "str",
                 "fwhm": "list_float",
                 "geo points": "list_float",
                 "header offset": "int",
                 "interleave": "str",
                 "lines": "int",
                 "map info": "list_str",
                 "pixel size": "list_str",
                 "projection info": "str",
                 "read procedures": "str",
                 "reflectance scale factor": "float",
                 "rpc info": "str",
                 "samples":"int",
                 "security tag": "str",
                 "sensor type": "str",
                 "smoothing factors": "list_float",
                 "solar irradiance": "float",
                 "spectra names": "list_str",
                 "sun azimuth": "float",
                 "sun elevation": "float",
                 "wavelength": "list_float",
                 "wavelength units": "str",
                 "x start": "float",
                 "y start": "float",
                 "z plot average": "str",
                 "z plot range": "str",
                 "z plot titles": "str"}

    headerDict = {}

    headerFile = open(hdrFile,'r')
    line = headerFile.readline()
      
    while line :
        if "=" in line:
            key,value = line.rstrip().split("=",1)
            # Add field not in ENVI default list
            if key.strip() not in fieldDict.keys():
                fieldDict[key.strip()] = "str"
            
            valType = fieldDict[key.strip()]
            
            if "{" in value and not "}" in value: 
                while "}" not in line:
                    line = headerFile.readline()
                    value+=line
    
            if '{}' in value: 
                value = np.nan
            elif valType == "list_float":
                value= np.array([float(x) for x in value.translate(str.maketrans("\n{}","   ")).split(",")])
            elif valType == "list_int":
                value= np.array([int(x) for x in value.translate(str.maketrans("\n{}","   ")).split(",")])
            elif valType == "list_str":
                value= [x.strip() for x in value.translate(str.maketrans("\n{}","   ")).split(",")]
            elif valType == "int":
                value = int(value.translate(str.maketrans("\n{}","   ")))
            elif valType == "float":
                value = float(value.translate(str.maketrans("\n{}","   ")))
            elif valType == "str":
                value = value.translate(str.maketrans("\n{}","   ")).strip().lower()

            headerDict[key.strip()] = value
                            
        line = headerFile.readline()
    
    # Fill unused fields with nans
    for key in fieldDict.keys():
        if key not in headerDict.keys():
            headerDict[key] = np.nan
    
    headerFile.close()
    return headerDict
        