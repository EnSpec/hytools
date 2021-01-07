# -*- coding: utf-8 -*-
"""
This module contains functions to apply the Modified topographic correction (SCS+C)
described in the following paper:

Richter, R., Kellenberger, T., & Kaufmann, H. (2009).
Comparison of topographic correction methods.
Remote Sensing, 1(3), 184-196.
https://doi.org/10.3390/rs1030184

Topographic correction consists of the following steps:


"""
import copy
import os
import numpy as np

def calc_modminn_coeffs(hy_obj):
    '''


    Args:
        hy_obj (TYPE): DESCRIPTION.

    Returns:
        None.

    '''
    hy_obj.topo = {'type': 'mod_minneart'}
    hy_obj.anc_data = {}

    cos_i = hy_obj.get_anc('cosine_i')
    i = np.rad2deg(np.arccos(cos_i))
    solar_zn = hy_obj.get_anc('solar_zn',radians=False)

    solar_zn_t = np.zeros(solar_zn.shape)
    solar_zn_t[solar_zn < 45] = solar_zn[solar_zn < 45] +20
    solar_zn_t[(solar_zn >= 45) & (solar_zn <= 55)] = solar_zn[(solar_zn >= 45) & (solar_zn <= 55)] +15
    solar_zn_t[solar_zn > 55] = solar_zn[solar_zn > 55] +10

    #Create NDVI mask to seperate vegetation
    ir = hy_obj.get_wave(850)
    red = hy_obj.get_wave(660)
    ndvi = (ir-red)/(ir+red)
    hy_obj.mask['ndvi'] = ndvi > 0.2

    c_factors = np.ones((2,hy_obj.lines,hy_obj.columns))
    c_factors[:] = cos_i/np.cos(np.radians(solar_zn_t))

    # Non vegetation correction factor
    c_factors[0][~hy_obj.mask['ndvi']] =  c_factors[0][~hy_obj.mask['ndvi']]**(1/2)
    c_factors[1][~hy_obj.mask['ndvi']] =  c_factors[1][~hy_obj.mask['ndvi']]**(1/2)

    # Vegetation correction factor
    c_factors[0][hy_obj.mask['ndvi']] =  c_factors[0][hy_obj.mask['ndvi']]**(3/4)
    c_factors[1][hy_obj.mask['ndvi']] =  c_factors[1][hy_obj.mask['ndvi']]**(1/3)

    #Adjust correction factors to prevent too strong correction
    c_factors[c_factors <.25] = .25
    c_factors[c_factors > 1] = 1

    #Correct pixels only where i > threshold
    c_factors[0][i < solar_zn_t] = 1
    c_factors[1][i < solar_zn_t] = 1

    c_factors[0][~hy_obj.mask['no_data']] = 1
    c_factors[1][~hy_obj.mask['no_data']] = 1

    hy_obj.anc_data['mm_c_factor'] = c_factors


def apply_modmin_band(hy_obj,band,index):
    '''


    Args:
        hy_obj (TYPE): DESCRIPTION.
        band (TYPE): DESCRIPTION.
        index (TYPE): DESCRIPTION.

    Returns:
        band (TYPE): DESCRIPTION.

    '''

    if hy_obj.wavelengths[index] >=720:
        cf_index = 1
    else:
        cf_index = 0
    band = band / hy_obj.anc_data['mm_c_factor'][cf_index]
    return band

