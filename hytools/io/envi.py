"""
 Functions for reading ENVI formatted binary files
See https://www.l3harrisgeospatial.com/docs/ENVIHeaderFiles.html for formatting details

TODO: Implement opening of ENVI files with different byte order
"""
import numpy as np

# ENVI datatype conversion dictionary
dtype_dict = {1:np.uint8,
             2:np.int16,
             3:np.int32,
             4:np.float32,
             5:np.float64,
             12:np.uint16,
             13:np.uint32,
             14:np.int64,
             15:np.uint64}

# Dictionary of all ENVI header fields
field_dict = {"acquisition time": "str",
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




class WriteENVI:
    """Iterator class for writing to an ENVI data file.

    """
    def __init__(self,output_name,header_dict):
        """Constructor method

        :param output_name: Pathname of output ENVI data file
        :type output_name: str
        :param header_dict: Dictionary containing ENVI header information
        :type header_dict: dict
        """

        self.interleave = header_dict['interleave']
        self.header_dict = header_dict
        self.output_name =output_name
        dtype = dtype_dict[header_dict["data type"]]
        lines = header_dict['lines']
        columns = header_dict['samples']
        bands = header_dict['bands']

        if self.interleave == "bip":
            self.data = np.memmap(output_name,dtype = dtype,
                                  mode='w+', shape = (lines,columns,bands))
        elif self.interleave == "bil":
            self.data = np.memmap(output_name,dtype = dtype,
                                  mode='w+', shape =(lines,bands,columns))
        elif self.interleave == "bsq":
            self.data = np.memmap(output_name,dtype = dtype,
                                  mode='w+',shape =(bands,lines,columns))
        write_envi_header(self.output_name,self.header_dict)

    def write_line(self,line,index):
        """Write line to ENVI file.

        :param line: Line array (columns,bands)
        :type line: numpy.ndarray
        :param index: Zero-based line index
        :type index: int
        """

        if self.interleave == "bip":
            self.data[index,:,:] = line

        elif self.interleave == "bil":
            self.data[index,:,:] = line

        elif self.interleave == "bsq":
            self.data[:,index,:] = line

    def write_column(self,column,index):
        """Write column to ENVI file.

        :param column: Column array (lines,bands)
        :type column: numpy.ndarray
        :param index: Zero-based column index
        :type index: int
        """

        if self.interleave == "bip":
            self.data[:,index,:]  = column
        elif self.interleave == "bil":
            self.data[:,:,index] = column
        elif self.interleave == "bsq":
            self.data[:,:,index] = column

    def write_band(self,band,index):
        """Write band to ENVI file.

        :param band: Band array (lines,columns)
        :type band: numpy.ndarray
        :param index: Zero-based band index
        :type index: int
        """
        if self.interleave == "bip":
            self.data[:,:,index]  = band
        elif self.interleave == "bil":
            self.data[:,index,:] = band
        elif self.interleave == "bsq":
            self.data[index,:,:]= band

    def write_chunk(self,chunk,line_index,column_index):
        """Write chunk to ENVI file.

        :param line: Chunks array (chunk lines,chunk columns,bands)
        :type line: numpy.ndarray
        :param line_index: Zero-based upper line index
        :type line_index: int
        :param column_index: Zero-based left column index
        :type column_index: int
        """
        x_start = column_index
        x_end = column_index + chunk.shape[1]
        y_start = line_index
        y_end = line_index + chunk.shape[0]

        if self.interleave == "bip":
            self.data[y_start:y_end,x_start:x_end,:] = chunk
        elif self.interleave == "bil":
            self.data[y_start:y_end,:,x_start:x_end] = np.moveaxis(chunk,-1,-2)
        elif self.interleave == "bsq":
            self.data[:,y_start:y_end,x_start:x_end] = np.moveaxis(chunk,-1,0)

    def close(self):
        """Delete numpy memmap
        """
        del self.data



def envi_header_from_hdf(hy_obj, interleave = 'bil'):
    """Create an ENVI header dictionary from HDF metadata

    :param hy_obj: Populated HyTools class object
    :type hy_obj: HyTools class object
    :param interleave: Date interleave type, defaults to 'bil'
    :type interleave: str, optional
    :return: ENVI header dictionary
    :rtype: dict
    """

    header_dict = {}
    header_dict["ENVI description"] = "{}"
    header_dict["samples"] = hy_obj.columns
    header_dict["lines"]   = hy_obj.lines
    header_dict["bands"]   = hy_obj.bands
    header_dict["header offset"] = 0
    header_dict["file type"] = "ENVI Standard"
    header_dict["data type"] = 2
    header_dict["interleave"] = interleave
    header_dict["sensor type"] = ""
    header_dict["byte order"] = 0
    header_dict["map info"] = hy_obj.map_info
    header_dict["coordinate system string"] = hy_obj.projection
    header_dict["wavelength units"] = hy_obj.wavelength_units
    header_dict["data ignore value"] =hy_obj.no_data
    header_dict["wavelength"] =hy_obj.wavelengths
    return header_dict


def write_envi_header(output_name,header_dict):
    """Write ENVI header to disk.

    :param output_name: Header file pathname
    :type output_name: str
    :param header_dict: ENVI header dictionary
    :type header_dict: dict
    """

    header_file = open(output_name + ".hdr",'w+')
    header_file.write("ENVI\n")

    for key in header_dict.keys():
        value = header_dict[key]
        # Convert list to comma seperated strings
        if isinstance(value,(list,np.ndarray)):
            value = "{%s}" % ",".join(map(str, value))
        else:
            value = str(value)
        # Skip entires with nan as value
        if value != 'None':
            header_file.write("%s = %s\n" % (key,value))
    header_file.close()



def envi_header_dict():
    """Returns empty dictionary with all ENVI header data fields.

    :return: Empty ENVI header dictionary
    :rtype: dict
    """
    return {key:None for (key,value) in field_dict.items()}


def envi_read_line(data,index,interleave):
    """Read line from ENVI file.

    :param data: hyTools data method
    :type data: hyTools data method
    :param index: Zero based line index
    :type index: int
    :param interleave: Data interleave type
    :type index: str
    :return: Line array (columns,bands)
    :rtype: numpy.ndarray
    """
    if interleave == "bip":
        line = data[index,:,:]
    elif interleave == "bil":
        line = data[index,:,:]
    elif interleave == "bsq":
        line = data[:,index,:]
    return line

def envi_read_column(data,index,interleave):
    """Read column from ENVI file.

    :param data: hyTools data method
    :type data: hyTools data method
    :param index: Zero based line index
    :type index: int
    :param interleave: Data interleave type
    :type index: str
    :return: Column array (lines,bands)
    :rtype: numpy.ndarray
    """
    if interleave == "bip":
        column = data[:,index,:]
    elif interleave == "bil":
        column = data[:,:,index]
    elif interleave == "bsq":
        column =  data[:,:,index]
    return column

def envi_read_band(data,index,interleave):
    """Read band from ENVI file.

    :param data: hyTools data method
    :type data: hyTools data method
    :param index: Zero based line index
    :type index: int
    :param interleave: Data interleave type
    :type index: str
    :return: Band array (lines,columns)
    :rtype: numpy.ndarray
    """
    if interleave == "bip":
        band =  data[:,:,index]
    elif interleave == "bil":
        band = data[:,index,:]
    elif interleave == "bsq":
        band = data[index,:,:]
    return band


def envi_read_chunk(data,col_start,col_end,line_start,line_end,interleave):
    """Read chunk from ENVI file.

    :param data: hyTools data method
    :type data: hyTools data method
    :param col_start: Zero based left column index
    :type col_start: int
    :param col_end: Zero based right column index
    :type col_end: int
    :param line_start: Zero based top line index
    :type line_start: int
    :param line_end: Zero based bottom line index
    :type line_end: int
    :param interleave: Data interleave type
    :type interleave: str
    :return: Chunk array
    :rtype: numpy.ndarray
    """
    if interleave == "bip":
        chunk = data[line_start:line_end,col_start:col_end,:]
    elif interleave == "bil":
        chunk = np.moveaxis(data[line_start:line_end,:,col_start:col_end],-1,-2)
    elif interleave == "bsq":
        chunk = np.moveaxis(data[:,line_start:line_end,col_start:col_end],0,-1)
    return chunk



def parse_envi_header(header_file):
    """Parse ENVI header file into dictionary.

    :param header_file: Header file pathname
    :type header_file: str
    :return: Populated header dictionary
    :rtype: dict
    """
    header_dict = envi_header_dict()
    header_file = open(header_file,'r')
    line = header_file.readline()

    while line :
        if "=" in line:
            key,value = line.rstrip().split("=",1)
            # Add field not in ENVI default list
            if key.strip() not in field_dict.keys():
                field_dict[key.strip()] = "str"
            val_type = field_dict[key.strip()]

            if "{" in value and not "}" in value:
                while "}" not in line:
                    line = header_file.readline()
                    value+=line

            if '{}' in value:
                value = None
            elif val_type == "list_float":
                value= np.array([float(x) for x in value.translate(str.maketrans("\n{}","   ")).split(",")])
            elif val_type == "list_int":
                value= np.array([int(x) for x in value.translate(str.maketrans("\n{}","   ")).split(",")])
            elif val_type == "list_str":
                value= [x.strip() for x in value.translate(str.maketrans("\n{}","   ")).split(",")]
            elif val_type == "int":
                value = int(value.translate(str.maketrans("\n{}","   ")))
            elif val_type == "float":
                value = float(value.translate(str.maketrans("\n{}","   ")))
            elif val_type == "str":
                value = value.translate(str.maketrans("\n{}","   ")).strip().lower()

            header_dict[key.strip()] = value
        line = header_file.readline()

    # Fill unused fields with None
    for key in field_dict:
        if key not in header_dict.keys():
            header_dict[key] = None

    header_file.close()
    return header_dict
