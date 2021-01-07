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
from .scs_c import calc_c
from ..io.envi import WriteENVI

def calc_cosine_coeffs(hy_obj):
    '''

    Args:
        hy_obj (TYPE): DESCRIPTION.

    Returns:
        None.

    '''
    topo_coeffs = {}
    topo_coeffs['type'] = 'cosine'
    topo_coeffs['coeffs'] = {}
    cosine_i = hy_obj.get_anc("cosine_i")

    for band_num,band in enumerate(hy_obj.bad_bands):
        if ~band:
            band = hy_obj.get_band(band_num,mask='topo')
            topo_coeffs['coeffs'][band_num] = calc_c(band,cosine_i[hy_obj.mask['topo']])

    hy_obj.topo = topo_coeffs


def apply_cosine_band(hy_obj,band,index):
    '''

    Args:
        hy_obj (TYPE): DESCRIPTION.
        band (TYPE): DESCRIPTION.
        index (TYPE): DESCRIPTION.

    Returns:
        band (TYPE): DESCRIPTION.

    '''
    cos_sz = np.cos(hy_obj.get_anc('solar_zn'))
    cos_i = hy_obj.get_anc('cosine_i')
    C = hy_obj.topo['coeffs'][index]
    correctionFactor = (cos_sz + C)/(cos_i + C)
    band = band * correctionFactor
    band[~hy_obj.mask['no_data']] = hy_obj.no_data

    return band





