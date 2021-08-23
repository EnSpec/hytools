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


REFRACTIVE_INDICES = np.array([
    [200, 1.396],
    [225, 1.373],
    [250, 1.362],
    [275, 1.354],
    [300, 1.349],
    [325, 1.346],
    [350, 1.343],
    [375, 1.341],
    [400, 1.339],
    [425, 1.338],
    [450, 1.337],
    [475, 1.336],
    [500, 1.335],
    [525, 1.334],
    [550, 1.333],
    [575, 1.333],
    [600, 1.332],
    [625, 1.332],
    [650, 1.331],
    [675, 1.331],
    [700, 1.331],
    [725, 1.33],
    [750, 1.33],
    [775, 1.33],
    [800, 1.329],
    [825, 1.329],
    [850, 1.329],
    [875, 1.328],
    [900, 1.328],
    [925, 1.328],
    [950, 1.327],
    [975, 1.327],
    [1000, 1.327],
    [1200, 1.324],
    [1400, 1.321],
    [1600, 1.317],
    [1800, 1.312],
    [2000, 1.306],
    [2200, 1.296],
    [2400, 1.279],
    [2600, 1.242],
    [2650, 1.219],
    [2700, 1.188],
    [2750, 1.157],
    [2800, 1.142],
    [2850, 1.149],
    [2900, 1.201],
    [2950, 1.292],
    [3000, 1.371]
])


def apply_gao_2021_correction(hy_obj, data, dimension, index):
    """
    Glint correction algorithm following:

    Gao BC, Li RR.
    Correction of Sunglint Effects in High Spatial Resolution
    Hyperspectral Imagery Using SWIR or NIR Bands and Taking Account of
    Spectral Variation of Refractive Index of Water.
    Adv Environ Eng Res 2021;2(3):16; doi:10.21926/aeer.2103017.
    """

    hy_obj.glint['correction_band'] = hy_obj.wave_to_band(
        hy_obj.glint['correction_wave']
    )

    if 'apply_glint' not in hy_obj.mask:
        hy_obj.gen_mask(mask_create,'apply_glint',hy_obj.glint['apply_mask'])

    if 'gao_b_simu' not in hy_obj.ancillary:
        hy_obj.ancillary['gao_b_simu'] = get_b_simu(hy_obj)

    if 'gao_rto' not in hy_obj.ancillary:
        hy_obj.ancillary['gao_rto'] = get_rto(hy_obj)

    if dimension == 'line':
        rto_line = hy_obj.ancillary['gao_rto'][index, :]
        rto_line = np.reshape(rto_line, (len(rto_line), 1))
        correction = rto_line * hy_obj.ancillary['gao_b_simu']

    elif dimension == 'column':
        rto_col = hy_obj.ancillary['gao_rto'][:, index]
        rto_col = np.reshape(rto_col, (len(rto_col), 1))
        correction = rto_col * hy_obj.ancillary['gao_b_simu']

    elif (dimension == 'band'):
        correction = (
            hy_obj.ancillary['gao_b_simu'][0, :][index]
            * hy_obj.ancillary['gao_rto']
        )

    elif dimension == 'chunk':
        x1, x2, y1, y2 = index
        rto_chunk = hy_obj.ancillary['gao_rto'][y1:y2, x1:x2]
        rto_chunk = np.reshape(
            rto_chunk,
            (
                rto_chunk.shape[0],
                rto_chunk.shape[1],
                1
            )
        )
        correction = rto_chunk * hy_obj.ancillary['gao_b_simu']

    elif dimension == 'pixels':
        y, x = index
        rto_pixels = hy_obj.ancillary['gao_rto'][y, x]
        rto_pixels = np.reshape(rto_pixels, (len(rto_pixels), 1))
        correction = rto_pixels * hy_obj.ancillary['gao_b_simu']

    return data - correction


def zenith_refracted(theta, n):
    """
    Find zenith of the outgoing reflected light
    n is the refractive index of water at a specific wavelength
    """
    theta_p = np.degrees(
        np.arcsin(np.sin(np.radians(theta)) / n)
    )

    return theta_p


def fresnel_reflectence(theta, theta_p):
    """
    Uses the fresnel equation to find the
    percentege of incident light reflected
    """
    theta_rad = np.radians(theta)
    theta_p_rad = np.radians(theta_p)

    return (
        (
            (np.sin(theta_rad - theta_p_rad)**2)
            / (np.sin(theta_rad + theta_p_rad)**2)
        ) + (
            (np.tan(theta_rad - theta_p_rad)**2)
            / (np.tan(theta_rad + theta_p_rad)**2)
        )
    ) / 2


def fresnel_spectra(theta, xs, ns):
    """
    Solves for the spectrum of reflected light
    according to fresnels equations
    """
    spectra = []
    for x in xs:
        n = np.interp(x, ns[:, 0], ns[:, 1])
        theta_p = zenith_refracted(theta, n)
        spectra.append(fresnel_reflectence(theta, theta_p))

    return np.array(spectra)


def get_b_simu(hy_obj):
    b_simu = fresnel_spectra(
        10**-5,
        hy_obj.wavelengths,
        REFRACTIVE_INDICES
    )

    return np.reshape(
        b_simu, (1, len(b_simu))
    )


def get_rto(hy_obj):
    b_ref = hy_obj.get_wave(hy_obj.glint['correction_wave'])

    b_ref_min = np.percentile(
        b_ref[
            (hy_obj.mask['apply_glint'])
            & (b_ref > 0)
        ],
        .0001
    )
    b_ref = b_ref - b_ref_min

    rto = (
        b_ref
        / hy_obj.ancillary['gao_b_simu'][0, :][hy_obj.glint['correction_band']]
    )
    rto[~hy_obj.mask['apply_glint']] = 0

    return rto
