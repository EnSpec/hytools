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
from ..io.envi import WriteENVI

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
        slope, intercept = np.linalg.lstsq(X, data,rcond=-1)[0].flatten()
    elif fit_type == 'NN':
        slope, intercept = nnls(X, data)[0].flatten()

    # Eq 8. Soenen et al. 2005
    c= intercept/slope

    # Set a large number if slope is zero
    if not np.isfinite(c):
        c = 100000.0
    return c

def calc_c_coeffs(hy_obj,topo_dict):
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

def apply_c(hy_obj,data,dimension,index):
    ''' Apply SCSS correction to a slice of the data

    Args:
        hy_obj (TYPE): DESCRIPTION.
        band (TYPE): DESCRIPTION.
        index (TYPE): DESCRIPTION.

    Returns:
        band (TYPE): DESCRIPTION.

    '''

    if 'cos_sz' not in hy_obj.ancillary.keys():
        cos_sz = np.cos(hy_obj.get_anc('solar_zn'))
        hy_obj.ancillary['cos_sz'] = cos_sz
    if 'cosine_i' not in hy_obj.ancillary.keys():
        cosine_i = hy_obj.cosine_i()
        hy_obj.ancillary['cosine_i'] = cosine_i

    C_bands = list(hy_obj.topo['coeffs'].keys())
    C = np.array(list(hy_obj.topo['coeffs'].values()))

    #Convert to float
    data = data.astype(np.float32)

    if dimension == 'line':
        #index= 3000
        #data = hy_obj.get_line(3000)
        data = data[:,C_bands]
        mask = hy_obj.mask['apply_topo'][index,:]
        cosine_i = hy_obj.ancillary['cosine_i'][[index],:].T
        cos_sz = hy_obj.ancillary['cos_sz'][[index],:].T
        correction_factor = (cos_sz + C)/(cosine_i + C)
        data[mask,:] = data[mask,:]*correction_factor[mask,:]

    elif dimension == 'column':
        # index= 300
        # data = hy_obj.get_column(index)
        data = data[:,C_bands]
        mask = hy_obj.mask['apply_topo'][:,index]
        cosine_i = hy_obj.ancillary['cosine_i'][:,[index]]
        cos_sz = hy_obj.ancillary['cos_sz'][:,[index]]
        correction_factor = (cos_sz + C)/(cosine_i + C)
        data[mask,:] = data[mask,:]*correction_factor[mask,:]

    elif dimension == 'band':
        #index= 8
        #data = hy_obj.get_band(index)
        C = hy_obj.topo['coeffs'][index]
        correction_factor = (hy_obj.ancillary['cos_sz'] + C)/(hy_obj.ancillary['cosine_i'] + C)
        data[hy_obj.mask['apply_topo']] = data[hy_obj.mask['apply_topo']] * correction_factor[hy_obj.mask['apply_topo']]

    elif dimension == 'chunk':
        # index = 200,501,3000,3501
        x1,x2,y1,y2 = index
        # data = hy_obj.get_chunk(x1,x2,y1,y2)
        data = data[:,:,C_bands]
        mask = hy_obj.mask['apply_topo'][y1:y2,x1:x2]
        cosine_i = hy_obj.ancillary['cosine_i'][y1:y2,x1:x2][:,:,np.newaxis]
        cos_sz = hy_obj.ancillary['cos_sz'][y1:y2,x1:x2][:,:,np.newaxis]
        correction_factor = (cos_sz + C)/(cosine_i + C)
        data[mask,:] = data[mask,:]*correction_factor[mask,:]

    elif dimension == 'pixels':
        # index = [[2000,2001],[200,501]]
        y,x = index
        # data = hy_obj.get_pixels(y,x)
        data = data[:,C_bands]
        mask = hy_obj.mask['apply_topo'][y,x]
        cosine_i = hy_obj.ancillary['cosine_i'][[y],[x]].T
        cos_sz = hy_obj.ancillary['cos_sz'][[y],[x]].T
        correction_factor = (cos_sz + C)/(cosine_i + C)
        data[mask,:] = data[mask,:]*correction_factor[mask,:]

    return data



