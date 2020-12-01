"""
This module contains functions to apply a topographic correction (SCS+C)
described in the following papers:

Scott A. Soenen, Derek R. Peddle,  & Craig A. Coburn (2005).
SCS+C: A Modified Sun-Canopy-Sensor Topographic Correction in Forested Terrain.
IEEE Transactions on Geoscience and Remote Sensing, 43(9), 2148-2159.

TOPO correction consists of the following steps:

    1. calculate incidence angle if it is not provided
    2. estimate C-Correction value
    3. apply C-Correction value to the image
"""

import numpy as np
from ..io import *

def gen_cosine_i(solar_zn, solar_az, aspect ,slope):
    ''' Generate cosine i image, all inputs in radians

    :param solar_zn: Solar zenith angle
    :type solar_zn: numpy.ndarray, float
    :param solar_az: Solar azimuth angle
    :type solar_az: numpy.ndarray, float
    :param aspect: Ground aspect
    :type aspect: numpy.ndarray, float
    :param slope: Ground slope
    :type slope: numpy.ndarray, float
    :return: Cosine i
    :rtype: numpy.ndarray, float
    '''
    relative_az = aspect - solar_az
    cosine_i = np.cos(solar_zn)*np.cos(slope) + np.sin(solar_zn)*np.sin(slope)*  np.cos(relative_az)

    return cosine_i

def calc_topo_c(band,cos_i):
    '''Return the topographic correction coefficients for the input data.

    :param band: Input band array
    :type band: numpy.ndarray
    :param mask: Mask
    :type mask: numpy.ndarray
    :param cos_i: Cosine i array
    :type cos_i: numpy.ndarray
    :return: Topographic correction coefficient
    :rtype: np.ndarray
    '''
    # Reshape for regression
    cos_i = np.expand_dims(cos_i,axis=1)
    X = np.concatenate([cos_i,np.ones(cos_i.shape)],axis=1)

    # Eq 7. Soenen et al., IEEE TGARS 2005
    slope, intercept = np.linalg.lstsq(X, band)[0].flatten()
    # Eq 8. Soenen et al., IEEE TGARS 2005
    C= intercept/slope

    # Set a large number if slope is zero
    if not np.isfinite(C):
        C = 100000.0
    return C

def topo_correct_img(hy_obj,output_name,cos_i = None):

    """Topographically correct an image.

    Parameters
    ----------
    hy_obj:     hyTools data object
                Data spectrum.
    output_name: str
                Path name for TOPO corrected file.
    cos_i:  np.array
                The cosine of the incidence angle (i ),
                defined as the angle between the normal to the pixel surface and the solar zenith direction

    Returns
    -------
    None

    """

    # Generate the cosine i
    # the cosine of the incidence angle (i ), defined as the angle between the normal to the pixel surface and the solar zenith direction;
    if cos_i is None:
        print("Calculating incidence angle...")
        cos_i =  calc_cosine_i(hy_obj.solar_zn, hy_obj.solar_az, hy_obj.aspect , hy_obj.slope)

    # Eq 11. Soenen et al., IEEE TGARS 2005
    # cos(alpha)* cos(theta)
    # alpha -- slope (slope), theta -- solar zenith angle (solar_zn)
    c1 = np.cos(hy_obj.solar_zn) * np.cos(hy_obj.slope)

    #Calcualate topographic correction coefficients for all bands
    topo_df = generate_topo_coeffs_img(hy_obj,cos_i = None)

    # Create writer object
    if  hy_obj.file_type == "ENVI":
        writer = WriteENVI(output_name,hy_obj.header_dict)
    elif hy_obj.file_type == "HDF":
        writer = None
    else:
        print("ERROR: File format not recognized.")

    iterator = hy_obj.iterate(by = 'chunk')

    while not iterator.complete:
        chunk = iterator.read_next()
        line_start =iterator.current_line
        line_end = iterator.current_line + chunk.shape[0]
        col_start = iterator.current_column
        col_end = iterator.current_column + chunk.shape[1]

        # Get C-Correction factor for chunks
        cos_i_chunk = cos_i[line_start:line_end,col_start:col_end]
        c1_chunk = c1[line_start:line_end,col_start:col_end]

        # Apply TOPO correction
        # Eq 11. Soenen et al., IEEE TGARS 2005
        c_factor = (c1_chunk[:,:,np.newaxis]+topo_df.c.values)/(cos_i_chunk[:,:,np.newaxis] + topo_df.c.values)
        topo_chunk = chunk*c_factor

        # Reassign no_data values
        topo_chunk[chunk == hy_obj.no_data] = hy_obj.no_data

        writer.write_chunk(topo_chunk,iterator.current_line,iterator.current_column)

    writer.close()
