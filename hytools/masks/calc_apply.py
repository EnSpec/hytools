# -*- coding: utf-8 -*-
"""
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

This module contain functions for generating boolean masks specific to apply image corrections and
models.
"""

from scipy.ndimage.morphology import binary_erosion
import numpy as np
from .cloud import zhai_cloud


def ndi(hy_obj,args):
    mask = hy_obj.ndi(args['band_1'],args['band_2'])
    mask = (mask >= float(args['min'])) & (mask <= float(args['max']))
    return mask


def ancillary(hy_obj,args):
    ''' Mask ancillary datasets based off min and max threshold

    '''
    if args['name'] == 'cosine_i':
        mask= hy_obj.cosine_i()
    else:
        mask = hy_obj.get_anc(args['name'])
    mask = (mask >= float(args['min'])) & (mask <= float(args['max']))
    return mask


def neon_edge(hy_obj,args):
    '''
    Mask artifacts in NEON images around edges.
    '''
    radius =args['radius']
    y_grid, x_grid = np.ogrid[-radius: radius + 1, -radius: radius + 1]
    window =  (x_grid**2 + y_grid**2 <= radius**2).astype(np.float)
    buffer_edge = binary_erosion(hy_obj.mask['no_data'], window).astype(bool)
    return buffer_edge


def kernel_finite(hy_obj,args):
    '''
    Create NDVI bin class mask
    '''
    k_vol  = hy_obj.volume_kernel(hy_obj.brdf['volume'])
    k_geom = hy_obj.geom_kernel(hy_obj.brdf['geometric'],
                                b_r=hy_obj.brdf["b/r"],
                                h_b =hy_obj.brdf["h/b"])
    mask = np.isfinite(k_vol) & np.isfinite(k_geom)
    return mask


def cloud(hy_obj,args):
    if args['method'] == 'zhai_2018':
        mask = ~zhai_cloud(hy_obj,args['cloud'],args['shadow'],
                args['T1'], args['t2'], args['t3'],
                args['t4'], args['T7'], args['T8'])

    return mask


def water(hy_obj,args):
    '''
    Create water mask using NDWI threshold
    '''
    mask = hy_obj.ndi(args['band_1'],args['band_2'])
    mask = mask >= float(args['threshold'])

    mask = binary_erosion(mask)

    return mask


def external(hy_obj,args):
    '''Load a mask from an external dataset
    '''

    hy_obj.anc_path['external_mask'] = [args['files'][hy_obj.file_name], 0]
    mask = hy_obj.get_anc('external_mask') == args['class']

    return mask











