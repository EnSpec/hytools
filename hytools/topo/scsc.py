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
from .c import calc_c
from ..masks import mask_dict

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

def calc_scsc_coeffs(hy_obj,topo_dict):
    '''

    Args:
        hy_obj (TYPE): DESCRIPTION.

    Returns:
        None.

    '''

    topo_dict['coeffs'] = {}
    cosine_i = hy_obj.cosine_i()

    for band_num,band in enumerate(hy_obj.bad_bands):
        if ~band:
            band = hy_obj.get_band(band_num,mask='calc_topo')
            topo_dict['coeffs'][band_num] = calc_c(band,cosine_i[hy_obj.mask['calc_topo']])
    hy_obj.topo = topo_dict

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
    cosine_i = hy_obj.cosine_i()

    C = hy_obj.topo['coeffs'][index]
    correction_factor = (c1 + C)/(cosine_i + C)
    band[hy_obj.mask['calc_topo']] = band[hy_obj.mask['calc_topo']] * correction_factor[hy_obj.mask['calc_topo']]
    band[~hy_obj.mask['no_data']] = hy_obj.no_data

    return band

def apply_scsc(hy_obj,data,dimension,index):
    ''' Apply SCSS correction to a slice of the data

    Args:
        hy_obj (TYPE): DESCRIPTION.
        band (TYPE): DESCRIPTION.
        index (TYPE): DESCRIPTION.

    Returns:
        band (TYPE): DESCRIPTION.

    '''

    if 'c1' not in hy_obj.ancillary.keys():
        c1 = np.cos(hy_obj.get_anc('slope')) * np.cos(hy_obj.get_anc('solar_zn'))
        hy_obj.ancillary['c1'] = c1
    if 'cosine_i' not in hy_obj.ancillary.keys():
        cosine_i = hy_obj.cosine_i()
        hy_obj.ancillary['cosine_i'] = cosine_i

    C_bands = list([int(x) for x in hy_obj.topo['coeffs'].keys()])
    C = np.array(list(hy_obj.topo['coeffs'].values()))

    #Convert to float
    data = data.astype(np.float32)
    hy_obj.topo['coeffs'] =  {int(k): hy_obj.topo['coeffs'][k] for k in hy_obj.topo['coeffs']}

    if (dimension != 'band') & (dimension != 'chunk'):
        if dimension == 'line':
            #index= 3000
            #data = hy_obj.get_line(3000)
            mask = hy_obj.mask['apply_topo'][index,:]
            cosine_i = hy_obj.ancillary['cosine_i'][[index],:].T
            c1 = hy_obj.ancillary['c1'][[index],:].T

        elif dimension == 'column':
            #index= 300
            #data = hy_obj.get_column(index)
            mask = hy_obj.mask['apply_topo'][:,index]
            cosine_i = hy_obj.ancillary['cosine_i'][:,[index]]
            c1 = hy_obj.ancillary['c1'][:,[index]]

        elif dimension == 'pixels':
            #index = [[2000,2001],[200,501]]
            y,x = index
            #data = hy_obj.get_pixels(y,x)
            mask = hy_obj.mask['apply_topo'][y,x]
            cosine_i = hy_obj.ancillary['cosine_i'][[y],[x]].T
            c1 = hy_obj.ancillary['c1'][[y],[x]].T

        correction_factor = np.ones(data.shape)
        correction_factor[:,C_bands] = (c1 + C)/(cosine_i + C)
        data[mask,:] = data[mask,:]*correction_factor[mask,:]

    elif dimension  == 'chunk':
        #index = 200,501,3000,3501
        x1,x2,y1,y2 = index
        #data = hy_obj.get_chunk(x1,x2,y1,y2)
        mask = hy_obj.mask['apply_topo'][y1:y2,x1:x2]
        cosine_i = hy_obj.ancillary['cosine_i'][y1:y2,x1:x2][:,:,np.newaxis]
        c1 = hy_obj.ancillary['c1'][y1:y2,x1:x2][:,:,np.newaxis]

        correction_factor = np.ones(data.shape)
        correction_factor[:,:,C_bands] = (c1 + C)/(cosine_i + C)
        data[mask,:] = data[mask,:]*correction_factor[mask,:]

    elif (dimension  == 'band') and (index in hy_obj.topo['coeffs']):
        #index= 8
        #data = hy_obj.get_band(index)
        C = hy_obj.topo['coeffs'][index]
        correction_factor = (hy_obj.ancillary['c1'] + C)/(hy_obj.ancillary['cosine_i'] + C)
        data[hy_obj.mask['apply_topo']] = data[hy_obj.mask['apply_topo']] * correction_factor[hy_obj.mask['apply_topo']]
    return data
