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

This module contains functions to apply the Modified topographic correction (SCS+C)
described in the following paper:

Richter, R., Kellenberger, T., & Kaufmann, H. (2009).
Comparison of topographic correction methods.
Remote Sensing, 1(3), 184-196.
https://doi.org/10.3390/rs1030184

Topographic correction consists of the following steps:


"""
import numpy as np

def calc_cosine_coeffs(hy_obj,topo_dict):
    '''

    Args:
        hy_obj (TYPE): DESCRIPTION.

    Returns:
        None.

    '''
    hy_obj.topo = topo_dict
    hy_obj.anc_data = {}

    cos_i = hy_obj.cosine_i()
    cos_solar_zn = np.cos(hy_obj.get_anc('solar_zn'))

    c_factor =  cos_solar_zn/cos_i
    c_factor[~hy_obj.mask['no_data']] = 1.
    hy_obj.ancillary['cosine_factor'] =c_factor

def apply_cosine(hy_obj,data,dimension,index):
    ''' Apply cosine correction to a slice of the data

    Args:
        hy_obj (TYPE): DESCRIPTION.
        band (TYPE): DESCRIPTION.
        index (TYPE): DESCRIPTION.

    Returns:
        band (TYPE): DESCRIPTION.

    '''

    if 'cosine_factor' not in hy_obj.ancillary.keys():
        calc_cosine_coeffs(hy_obj)

    #Convert to float
    data = data.astype(np.float32)

    if dimension == 'line':
        #index= 3000
        #data = hy_obj.get_line(3000)
        data = data*hy_obj.ancillary['cosine_factor'][np.newaxis,index,:]

    elif dimension == 'column':
        #index= 300
        #data = hy_obj.get_column(index)
        data = hy_obj.ancillary['cosine_factor'][:,index,np.newaxis]

    elif dimension == 'band':
        #index= 8
        #data = hy_obj.get_band(index)
        data = data * hy_obj.ancillary['cosine_factor']

    elif dimension == 'chunk':
        #index = 200,501,3000,3501
        x1,x2,y1,y2 = index
        #data = hy_obj.get_chunk(x1,x2,y1,y2)
        data  = data*hy_obj.ancillary['cosine_factor'][y1:y2,x1:x2][:,:,np.newaxis]

    elif dimension == 'pixels':
        #index = [[2000,2001],[200,501]]
        y,x = index
        #data = hy_obj.get_pixels(y,x)
        data = data*hy_obj.ancillary['cosine_factor'][y,x][:, np.newaxis]
    return data
