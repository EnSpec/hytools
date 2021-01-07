# -*- coding: utf-8 -*-
"""

"""
import numpy as np
from .minneart import apply_modmin_band
from .scs_c import apply_scsc_band
from .cosine import apply_cosine_band
from ..io.envi import WriteENVI


def calc_cosine_i(solar_zn, solar_az, aspect ,slope):
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
        cnumpy.ndarray: Cosine i image.

    """

    relative_az = aspect - solar_az
    cosine_i = np.cos(solar_zn)*np.cos(slope) + np.sin(solar_zn)*np.sin(slope)*  np.cos(relative_az)

    return cosine_i


def topo_correct_band(hy_obj,band,index):
    '''

    Args:
        hy_obj (TYPE): DESCRIPTION.
        band (TYPE): DESCRIPTION.
        index (TYPE): DESCRIPTION.

    Returns:
        band (TYPE): DESCRIPTION.

    '''
    if hy_obj.topo['type'] == 'mod_minneart':
        band = apply_modmin_band(hy_obj,band,index)

    elif hy_obj.topo['type']  == 'scs+c':
        band = apply_scsc_band(hy_obj,band,index)

    elif hy_obj.topo['type']  == 'cosine':
        band = apply_cosine_band(hy_obj,band,index)

    return band


