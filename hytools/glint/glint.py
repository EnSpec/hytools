# -*- coding: utf-8 -*-
"""
HyTools:  Hyperspectral image processing library
Copyright (C) 2021 University of Wisconsin

Authors: Evan Greenberg.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
import ray
from ..misc import set_glint
from .hochberg_2003 import apply_hochberg_2003_correction
from .gao_2021 import apply_gao_2021_correction
from .hedley_2005 import apply_hedley_2005_correction


def set_glint_parameters(actors, config_dict):
    # Assign glint dict
    glint_dict = config_dict['glint']

    # Set Glint dict
    _ = ray.get([
        a.do.remote(set_glint, glint_dict) for a in actors
    ])

    # Add glint correction
    _ = ray.get([
        a.do.remote(lambda x: x.corrections.append('glint')) for a in actors
    ])


def apply_glint_correct(hy_obj, data, dimension, index):
    ''' Corrects glint based on the specified algorithm in the config.
        Options include:
            Hochberg et al., 2003: hochberg
            Gao et al., 2021: gao
            Hedley et al. 2005: hedley
            ...
    '''

    # Perform one of the corrections
    if hy_obj.glint['type'] == 'hochberg':
        data = apply_hochberg_2003_correction(hy_obj, data, dimension, index)

    elif hy_obj.glint['type'] == 'gao':
        data = apply_gao_2021_correction(hy_obj, data, dimension, index)

    elif hy_obj.glint['type'] == 'hedley':
        data = apply_hedley_2005_correction(hy_obj, data, dimension, index)

    #Truncate reflectance values below 0
    if hy_obj.glint['truncate']:
        data[(data < 0) & (data != hy_obj.no_data)]= 0

    return data
