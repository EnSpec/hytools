# -*- coding: utf-8 -*-
'''
HyTools:  Hyperspectral image processing library
Copyright (C) 2021 University of Wisconsin

Authors: Adam Chlus, Zhiwei Ye, Philip Townsend.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

'''
from .calc_apply import *
from .cloud import *

mask_dict = {'ndi' : ndi,
             'neon_edge' : neon_edge,
             'kernel_finite': kernel_finite,
             'ancillary':  ancillary,
             'cloud': cloud,
             'water': water,
             'external' : external}

def mask_create(hy_obj,masks):
    ''' Combine a series of boolean masks using an
    and operator
    '''
    mask = np.copy(hy_obj.mask['no_data'])

    for mask_name,args in masks:
        mask &= mask_dict[mask_name](hy_obj,args)

    return mask

