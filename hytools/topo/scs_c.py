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

TODO: Rationale/ examples for using different fitting algorithms

"""
import numpy as np
from scipy.optimize import nnls
from ..misc import progbar

def calc_c(data,cosine_i,fit_type = 'OLS'):
    """Calculate the topographic correction coefficient (c) for the input data.
    Used for both the cosine and SCS+S topographic corrections.

    Args:
        band (numpy.ndarray): Image  array.
        cosine_i (numpy.ndarray): Cosine i array.
        fit_type (str): Linear model fitting type.

    Returns:
        numpy.ndarray: Topographic correction coefficient.

    """

    # Reshape for regression
    cosine_i = np.expand_dims(cosine_i,axis=1)
    X = np.concatenate([cosine_i,np.ones(cosine_i.shape)],axis=1)

    # Eq 7. Soenen et al. 2005
    if fit_type == 'OLS':
        slope, intercept = np.linalg.lstsq(X, data)[0].flatten()
    elif fit_type == 'NN':
        slope, intercept = nnls(X, data)[0].flatten()

    # Eq 8. Soenen et al. 2005
    c= intercept/slope

    # Set a large number if slope is zero
    if not np.isfinite(c):
        c = 100000.0
    return c

def calc_scsc_c1(solar_zn,slope):
    """ Calculate c1
        All input geometry units must be in radians.

    Args:
        solar_zn (numpy.ndarray): Solar zenith angle.
        slope (numpy.ndarray): Ground slope.

    Returns:
        numpy.ndarray: C1.

    """

    # Eq 11. Soenen et al. 2005
    scsc_c1 = np.cos(solar_zn) * np.cos(slope)
    return scsc_c1


def calc_scsc_coeffs(hy_obj):
    '''

    Args:
        hy_obj (TYPE): DESCRIPTION.

    Returns:
        None.

    '''
    topo_coeffs = {}
    topo_coeffs['type'] = 'scs+c'
    topo_coeffs['coeffs'] = {}
    cosine_i = hy_obj.get_anc("cosine_i")

    for band_num,band in enumerate(hy_obj.bad_bands):
        if ~band:
            print(band_num)
            band = hy_obj.get_band(band_num,mask='topo')
            topo_coeffs['coeffs'][band_num] = calc_c(band,cosine_i[hy_obj.mask['topo']])

    hy_obj.topo = topo_coeffs


def apply_scsc_band(hy_obj,band,index):
    '''

    Args:
        hy_obj (TYPE): DESCRIPTION.
        band (TYPE): DESCRIPTION.
        index (TYPE): DESCRIPTION.

    Returns:
        band (TYPE): DESCRIPTION.

    '''

    c1 = np.cos(hy_obj.get_anc('slope')) * np.cos(hy_obj.get_anc('solar_zn'))
    cos_i = hy_obj.get_anc('cosine_i')

    C = hy_obj.topo['coeffs'][index]
    correction_factor = (c1 + C)/(cos_i + C)
    band = band * correction_factor
    band[~hy_obj.mask['no_data']] = hy_obj.no_data

    return band
