# -*- coding: utf-8 -*-
"""Functions for reading and writing ENVI formatted binary files

Todo:
    * Implement opening of ENVI files with different byte order

"""
import os
from collections import Counter
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


def open_envi(hy_obj,anc_path = {}):
    """Open ENVI formated image file and populate Hytools object.


    Args:
        src_file (str): Pathname of input ENVI image file, header assumed to be located in
        same directory.
        anc_path (dict): Dictionary with pathnames and band numbers of ancillary datasets.

    Returns:
        HyTools file object: Populated HyTools file object.

    """

    if not os.path.isfile(os.path.splitext(hy_obj.file_name)[0] + ".hdr"):
        print("ERROR: Header file not found.")
        return None

    header_dict = parse_envi_header(os.path.splitext(hy_obj.file_name)[0] + ".hdr")
    hy_obj.lines =  header_dict["lines"]
    hy_obj.columns =  header_dict["samples"]
    hy_obj.bands =   header_dict["bands"]
    hy_obj.bad_bands = np.array([False for band in range(hy_obj.bands)])
    hy_obj.interleave =  header_dict["interleave"]
    hy_obj.fwhm =  header_dict["fwhm"]
    hy_obj.wavelengths = header_dict["wavelength"]
    hy_obj.wavelength_units = header_dict["wavelength units"]
    hy_obj.dtype = dtype_dict[header_dict["data type"]]
    hy_obj.no_data = header_dict['data ignore value']
    hy_obj.map_info = header_dict['map info']
    hy_obj.byte_order = header_dict['byte order']
    hy_obj.anc_path = anc_path

    if isinstance(header_dict['bbl'],np.ndarray):
        hy_obj.bad_bands = np.array([x==1 for x in header_dict['bbl']])
    if header_dict["interleave"] == 'bip':
        hy_obj.shape = (hy_obj.lines, hy_obj.columns, hy_obj.bands)
    elif header_dict["interleave"] == 'bil':
        hy_obj.shape = (hy_obj.lines, hy_obj.bands, hy_obj.columns)
    elif header_dict["interleave"] == 'bsq':
        hy_obj.shape = (hy_obj.bands, hy_obj.lines, hy_obj.columns)
    else:
        print("ERROR: Unrecognized interleave type.")
        hy_obj = None

    # If no_data value is not specified guess using image corners.
    if hy_obj.no_data is None:
        hy_obj.load_data()
        up_l = hy_obj.data[0,0,0]
        up_r = hy_obj.data[0,-1,0]
        low_l = hy_obj.data[-1,0,0]
        low_r = hy_obj.data[-1,-1,0]
        counts = {v: k for k, v in Counter([up_l,up_r,low_l,low_r]).items()}
        hy_obj.no_data = counts[max(counts.keys())]
        hy_obj.close_data()

    del header_dict
    return hy_obj


class WriteENVI:
    """Iterator class for writing to an ENVI data file.

    """
    def __init__(self,output_name,header_dict):
        """
        Args:
            output_name (str): Pathname of output ENVI data file.
            header_dict (dict): Dictionary containing ENVI header information.

        Returns:
            None.

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
        """
        Args:
            line (numpy.ndarray): Line array (columns,bands).
            index (int): Zero-based line index.

        Returns:
            None.

        """

        if self.interleave == "bip":
            self.data[index,:,:] = line

        elif self.interleave == "bil":
            self.data[index,:,:] = np.moveaxis(line,0,1)

        elif self.interleave == "bsq":
            self.data[:,index,:] = np.moveaxis(line,0,1)

    def write_column(self,column,index):
        """
        Args:
            column (numpy.ndarray): Column array (lines,bands).
            index (int): Zero-based column index.

        Returns:
            None.

        """

        if self.interleave == "bip":
            self.data[:,index,:]  = column
        elif self.interleave == "bil":
            self.data[:,:,index] = column
        elif self.interleave == "bsq":
            self.data[:,:,index] =  np.moveaxis(column,0,1)

    def write_band(self,band,index):
        """
        Args:
            band (numpy.ndarray): Band array (lines,columns).
            index (int): Zero-based band index.

        Returns:
            None.

        """

        if self.interleave == "bip":
            self.data[:,:,index]  = band
        elif self.interleave == "bil":
            self.data[:,index,:] = band
        elif self.interleave == "bsq":
            self.data[index,:,:]= band

    def write_chunk(self,chunk,line_index,column_index):
        """
        Args:
            chunk (TYPE): Chunks array (chunk lines,chunk columns,bands).
            line_index (int): Zero-based upper line index.
            column_index (int): Zero-based left column index.

        Returns:
            None.

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
        """Delete numpy memmap.
        """
        del self.data

def envi_header_from_neon(hy_obj, interleave = 'bsq'):
    """Create an ENVI header dictionary from NEON metadata

    Args:
        hy_obj (Hytools object): Populated HyTools file object.
        interleave (str, optional): Date interleave type. Defaults to 'bil'.

    Returns:
        dict: Populated ENVI header dictionary.

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
    """Write ENVI header file to disk.

    Args:
        output_name (str): Header file pathname.
        header_dict (dict): Populated ENVI header dictionary..

    Returns:
        None.

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
    """
    Returns:
        dict: Empty ENVI header dictionary.

    """
    return {key:None for (key,value) in field_dict.items()}


def envi_read_line(data,index,interleave):
    """
    Args:
        data (numpy.memmap): Numpy memory-map.
        index (int): Zero-based line index.
        interleave (str): Data interleave type.

    Returns:
        numpy.ndarray: Line array (columns, bands).

    """

    if interleave == "bip":
        line = data[index,:,:]
    elif interleave == "bil":
        line = np.moveaxis(data[index,:,:],0,1)
    elif interleave == "bsq":
        line = np.moveaxis(data[index,:,:],0,1)
    return line

def envi_read_column(data,index,interleave):
    """
    Args:
        data (numpy.memmap): Numpy memory-map.
        index (int): Zero-based column index.
        interleave (str): Data interleave type.

    Returns:
        numpy.ndarray: Column array (lines,bands).

    """

    if interleave == "bip":
        column = data[:,index,:]
    elif interleave == "bil":
        column = data[:,:,index]
    elif interleave == "bsq":
        column =  np.moveaxis(data[:,:,index],0,1)
    return column

def envi_read_band(data,index,interleave):
    """
    Args:
        data (numpy.memmap): Numpy memory-map.
        index (int): Zero-based line index.
        interleave (str): Data interleave type.

    Returns:
        numpy.ndarray: Band array (lines,columns).

    """

    if interleave == "bip":
        band =  data[:,:,index]
    elif interleave == "bil":
        band = data[:,index,:]
    elif interleave == "bsq":
        band = data[index,:,:]
    return band

def envi_read_pixels(data,lines,columns,interleave):
    """
    Args:
        data (numpy.memmap): Numpy memory-map.
        lines (list): List of zero-indexed line indices.
        columns (list): List of zero-indexed column indices.
        interleave (str): Data interleave type.

    Returns:
        numpy.ndarray: Pixel array (pixels,bands).

    """
    if interleave == "bip":
        pixels =  data[lines,columns,:]
    elif interleave == "bil":
        pixels = data[lines,:,columns]
    elif interleave == "bsq":
        pixels = data[:,lines,columns]
    return pixels


def envi_read_chunk(data,col_start,col_end,line_start,line_end,interleave):
    """
    Args:
        data (numpy.memmap): Numpy memory-map.
        col_start (int):  Zero-based left column index.
        col_end (int): Non-inclusive zero-based right column index.
        line_start (int): Zero -ased top line index.
        line_end (int): Non-inclusive zero-based bottom line index.
        interleave (str): Data interleave type.

    Returns:
        numpy.ndarray: Chunk array (line_end-line_start,col_end-col_start,bands).

    """

    if interleave == "bip":
        chunk = data[line_start:line_end,col_start:col_end,:]
    elif interleave == "bil":
        chunk = np.moveaxis(data[line_start:line_end,:,col_start:col_end],-1,-2)
    elif interleave == "bsq":
        chunk = np.moveaxis(data[:,line_start:line_end,col_start:col_end],0,-1)
    return chunk



def parse_envi_header(header_file):
    """
    Args:
        header_file (str): Header file pathname.

    Returns:
        dict: Populated header dictionary.

    """

    header_dict = envi_header_dict()
    header_file = open(header_file,'r')
    line = header_file.readline()

    while line :
        if "=" in line:
            key,value = line.rstrip().split("=",1)
            # Add fields not in ENVI default list
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
