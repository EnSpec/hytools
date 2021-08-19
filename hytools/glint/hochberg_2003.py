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
import numpy as np
from ..masks import mask_create


def apply_hochberg_2003_correction(hy_obj, data, dimension, index):
    """
    Glint correction algorithm following:

    Hochberg, EJ, Andréfouët, S and Tyler, MR. 2003.
    Sea surface correction of high spatial resolution Ikonos images to
    improve bottom mapping in near‐shore environments..
    IEEE Transactions on Geoscience and Remote Sensing, 41: 1724–1729.
    """

    if 'water' not in hy_obj.mask:
        hy_obj.mask['water'] = hy_obj.get_anc('water')

    if 'hochberg_correction' not in hy_obj.ancillary:
        hy_obj.ancillary['hochberg_correction'] = (
            get_hochberg_correction(hy_obj)
        )

    if dimension == 'line':
        correction = hy_obj.ancillary['hochberg_correction'][index, :]
        bandnums = len(hy_obj.wavelengths)
        correction = np.repeat(
            correction[:, np.newaxis],
            bandnums,
            axis=1
        )

    elif dimension == 'column':
        correction = hy_obj.ancillary['hochberg_correction'][:, index]
        bandnums = len(hy_obj.wavelengths)
        correction = np.repeat(
            correction[:, np.newaxis],
            bandnums,
            axis=1
        )

    elif (dimension == 'band'):
        correction = hy_obj.ancillary['hochberg_correction']

    elif dimension == 'chunk':
        # Get Index
        x1, x2, y1, y2 = index
        correction = hy_obj.ancillary['hochberg_correction'][y1:y2, x1:x2]
        bandnums = len(hy_obj.wavelengths)
        correction = np.repeat(
            correction[:, :, np.newaxis],
            bandnums,
            axis=2
        )

    elif dimension == 'pixels':
        y, x = index
        correction = hy_obj.ancillary['hochberg_correction'][y, x]
        bandnums = len(hy_obj.wavelengths)
        correction = np.repeat(
            correction[:, np.newaxis],
            bandnums,
            axis=1
        )

    return data - correction


def get_hochberg_correction(hy_obj):
    """
    Calculates the hochberg correction across entire image.
    Uses the NIR or SWIR wavelengths to find the amount of signal
    attributed to glint. Zeros out non-water pixels
    """

    nir_swir_array = np.copy(hy_obj.get_wave(hy_obj.glint['correction_wave']))

    nir_swir_array[~hy_obj.mask['water']] = 0

    nir_swir_min = np.percentile(
        nir_swir_array[nir_swir_array > 0], .001
    )

    hochberg_correction = nir_swir_array - nir_swir_min
    hochberg_correction[~hy_obj.mask['water']] = 0

    return hochberg_correction
