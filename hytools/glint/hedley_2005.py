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
from scipy import stats
from ..masks import mask_create


def apply_hedley_2005_correction(hy_obj, data, dimension, index):
    """
    Glint correction algorithm following:

    Hedley, J. D., Harborne, A. R., & Mumby, P. J. (2005).
    Simple and robust removal of sun glint for mapping shallow‐water benthos.
    International Journal of Remote Sensing, 26(10), 2107-2112.
    """
    # Raise exception is there is no deep water sample provided
    if isinstance(hy_obj.glint.get('deep_water_sample'), type(None)):
        raise KeyError("No Deep Water Sample Provided")

    hy_obj.glint['correction_band'] = hy_obj.wave_to_band(
        hy_obj.glint['correction_wave']
    )

    if 'apply_glint' not in hy_obj.mask:
        hy_obj.gen_mask(mask_create,'apply_glint',hy_obj.glint['apply_mask'])

    if 'hedley_slopes' not in hy_obj.ancillary:
        hy_obj.ancillary['hedley_slopes'] = optimize_slopes(hy_obj)

    if 'hedley_nir_swir_diff' not in hy_obj.ancillary:
        hy_obj.ancillary['hedley_nir_swir_diff'] = nir_swir_diff(hy_obj)

    if dimension == 'line':
        correction = (
            hy_obj.ancillary['hedley_nir_swir_diff'][index, :].reshape(-1, 1)
            * hy_obj.ancillary['hedley_slopes']
        )
        correction[~hy_obj.mask['apply_glint'][index, :], :] = 0

    elif dimension == 'column':
        correction = (
            hy_obj.ancillary['hedley_nir_swir_diff'][:, index].reshape(-1, 1)
            * hy_obj.ancillary['hedley_slopes']
        )
        correction[~hy_obj.mask['apply_glint'][:, index], :] = 0

    elif (dimension == 'band'):
        correction = (
            hy_obj.ancillary['hedley_nir_swir_diff']
            * hy_obj.ancillary['hedley_slopes'][0, index]
        )
        correction[~hy_obj.mask['apply_glint']] = 0

    elif dimension == 'chunk':
        x1, x2, y1, y2 = index
        corr_diff = hy_obj.ancillary['hedley_nir_swir_diff'][y1:y2, x1:x2]
        bandnums = data.shape[2]
        corr_diff = np.repeat(
            corr_diff[:, :, np.newaxis],
            bandnums,
            axis=2
        )

        correction = corr_diff * hy_obj.ancillary['hedley_slopes']
        correction[~hy_obj.mask['apply_glint'][y1:y2, x1:x2], :] = 0

    elif dimension == 'pixels':
        y, x = index

        correction = (
            hy_obj.ancillary['hedley_nir_swir_diff'][y, x].reshape(-1, 1)
            * hy_obj.ancillary['hedley_slopes']
        )
        correction[~hy_obj.mask['apply_glint'][y, x], :] = 0

    return data - correction


def optimize_slopes(hy_obj):
    deep_water = hy_obj.get_chunk(
        *hy_obj.glint['deep_water_sample'][hy_obj.file_name]
    )

    deep_correction = (
        deep_water[:, :, hy_obj.glint['correction_band']].flatten()
    )

    # Iterate through each band to find the band-slope
    slopes = np.empty([1, len(hy_obj.wavelengths)])
    for i, band in enumerate(hy_obj.wavelengths):
        # Get flattened deep water sample of band
        wave_num = np.argmin(
            np.abs(hy_obj.wavelengths - band)
        )
        wave = deep_water[:, :, wave_num].flatten()

        # Regress
        (
            slope,
            intercept,
            r_value,
            p_value,
            std_err
        ) = stats.linregress(
            deep_correction,
            wave
        )
        slopes[0, i] = slope

    return slopes


def nir_swir_diff(hy_obj):
    nir_swir_array = np.copy(
        hy_obj.get_wave(hy_obj.glint['correction_wave'])
    )
    nir_swir_array[~hy_obj.mask['apply_glint']] = 0
    nir_swir_min = np.percentile(nir_swir_array[nir_swir_array > 0], .0001)

    return nir_swir_array - nir_swir_min
