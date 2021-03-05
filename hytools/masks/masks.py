# -*- coding: utf-8 -*-
"""
The :mod:`hytools.masks` module include functions image correction.
"""
from .calc_apply import *
from .cloud import *

mask_dict = {'ndi' : ndi,
             'neon_edge' : neon_edge,
             'kernel_finite': kernel_finite,
             'ancillary':  ancillary,
             'cloud': cloud}

def mask_create(hy_obj,masks):
    ''' Combine a series of boolean masks using an
    and operator
    '''
    mask = np.copy(hy_obj.mask['no_data'])

    for mask_name,args in masks:
        mask &= mask_dict[mask_name](hy_obj,args)

    return mask

