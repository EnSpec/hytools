# -*- coding: utf-8 -*-
"""
This module contains functions to apply a topographic correction (SCS+C)
described in the following papers:

Scott A. Soenen, Derek R. Peddle,  & Craig A. Coburn (2005).
SCS+C: A Modified Sun-Canopy-Sensor Topographic Correction in Forested Terrain.
IEEE Transactions on Geoscience and Remote Sensing, 43(9), 2148-2159.
https://doi.org/10.1109/TGRS.2005.852480

Topographic correction consists of the following steps:

    1. calculate incidence angle if it is not provided
    2. estimate C-Correction value
    3. apply C-Correction value to the image data

"""

import numpy as np

def gen_cosine_i(solar_zn, solar_az, aspect ,slope):
    """Generate cosine i image. The cosine of the incidence angle (i) is
       defined as the angle between the normal to the pixel surface
       and the solar zenith direction.
       All input geometry units must be in radians.

    Args:
        solar_az (numpy.ndarray): Solar azimuth angle.
        solar_zn (numpy.ndarray): Solar zenith angle.
        aspect (numpy.ndarray): Ground aspect.
        slope (numpy.ndarray): Ground slope.

    Returns:
        cosine_i (numpy.ndarray): Cosine i image.

    """

    relative_az = aspect - solar_az
    cosine_i = np.cos(solar_zn)*np.cos(slope) + np.sin(solar_zn)*np.sin(slope)*  np.cos(relative_az)

    return cosine_i

def calc_scsc_c(data,cosine_i):
    """Calculate the topographic correction coefficient (c) for the input data.

    Args:
        band (numpy.ndarray): Image  array.
        cosine_i (numpy.ndarray): Cosine i array.

    Returns:
        scsc_c (numpy.ndarray): Topographic correction coefficient.

    """

    # Reshape for regression
    cosine_i = np.expand_dims(cosine_i,axis=1)
    X = np.concatenate([cosine_i,np.ones(cosine_i.shape)],axis=1)

    # Eq 7. Soenen et al. 2005
    slope, intercept = np.linalg.lstsq(X, data)[0].flatten()
    # Eq 8. Soenen et al. 2005
    scsc_c= intercept/slope

    # Set a large number if slope is zero
    if not np.isfinite(scsc_c):
        scsc_c = 100000.0
    return scsc_c

def calc_scsc_c1(solar_zn,slope):
    """ Calculate c1
        All input geometry units must be in radians.

    Args:
        solar_zn (numpy.ndarray): Solar zenith angle.
        slope (numpy.ndarray): Ground slope.

    Returns:
        scsc_c1 (numpy.ndarray): C1.

    """

    # Eq 11. Soenen et al. 2005
    scsc_c1 = np.cos(solar_zn) * np.cos(slope)
    return scsc_c1

