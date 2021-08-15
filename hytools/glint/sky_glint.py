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
import argparse
from ctypes import c_int,c_double
import json
from multiprocessing import Pool
import os
import sys
from typing import List

import hytools as ht
from hytools.brdf.kernels import calc_volume_kernel,calc_geom_kernel
from hytools.masks import mask_create, water
from hytools.misc import progbar, pairwise, set_glint

# from matplotlib import pyplot as plt
from matplotlib.widgets import Button
from matplotlib.patches import Rectangle
import pylab as plt

import numpy as np
import numpy.ctypeslib as npct

import ray
import ray.services

from scipy.linalg import inv, svd
from scipy.interpolate import interp1d, RegularGridInterpolator

from scipy.optimize import leastsq 
from scipy.optimize import minimize
from scipy.io import loadmat

os.environ["MKL_NUM_THREADS"] = "1"


class VectorInterpolator:
    """ Linear look up table interpolator.  
        Support linear interpolation through radial space by expanding the look
        up tables with sin and cos dimensions.

        Args:
            grid_input: list of lists of floats, indicating the gridpoint 
                        elements in each grid dimension
            data_input: n dimensional array of radiative transfer engine 
                        outputs (each dimension size corresponds to the
                        given grid_input list length, with the last 
                        dimensions equal to the number of sensor channels)
            lut_interp_types: a list indicating if each dimension is in 
                              radiance (r), degrees (r), or normal (n) units.

        Notes:
            Pulled from: 
        https://github.com/isofit/isofit/blob/master/isofit/core/common.py
        """

    def __init__(self, grid_input: List[List[float]], 
                 data_input: np.array, lut_interp_types: List[str]):
        self.lut_interp_types = lut_interp_types
        self.single_point_data = None

        # Lists and arrays are mutable, so copy first
        grid = grid_input.copy()
        data = data_input.copy()

        # Check if we are using a single grid point. If so, store the grid input.
        if np.prod(list(map(len, grid))) == 1:
            self.single_point_data = data

        # expand grid dimensionality as needed
        [radian_locations] = np.where(self.lut_interp_types == 'd')
        [degree_locations] = np.where(self.lut_interp_types == 'r')

        angle_locations = np.hstack([radian_locations, degree_locations])
        angle_types = np.hstack([
            self.lut_interp_types[radian_locations],
            self.lut_interp_types[degree_locations]
        ])

        for _angle_loc in range(len(angle_locations)):
            angle_loc = angle_locations[_angle_loc]

            # get original grid at given location
            original_grid_subset = np.array(grid[angle_loc])

            # convert for angular coordinates
            if (angle_types[_angle_loc] == 'r'):
                grid_subset_cosin = np.cos(original_grid_subset)
                grid_subset_sin = np.sin(original_grid_subset)

            elif (angle_types[_angle_loc] == 'd'):
                grid_subset_cosin = np.cos(original_grid_subset / 180. * np.pi)
                grid_subset_sin = np.sin(original_grid_subset / 180. * np.pi)

            # handle the fact that the grid may no longer be in order
            grid_subset_cosin_order = np.argsort(grid_subset_cosin)
            grid_subset_sin_order = np.argsort(grid_subset_sin)

            # convert current grid location, and add a second
            grid[angle_loc] = grid_subset_cosin[grid_subset_cosin_order]
            grid.insert(angle_loc+1, grid_subset_sin[grid_subset_sin_order])

            # Copy data through the extra dimension, at the specific angle_loc axes
            data = np.swapaxes(data, -1, angle_loc)
            data_dim = list(np.shape(data))
            data_dim.append(data_dim[-1])
            data = data[..., np.newaxis] * np.ones(data_dim)

            # Copy the data between the first two axes
            for ind in range(data.shape[-1]):
                data[..., ind] = data[..., :, ind]

            # Now re-order the cosin dimension
            data = data[..., grid_subset_cosin_order, :]
            # Now re-order the sin dimension
            data = data[..., grid_subset_sin_order]

            # now re-arrange the axes so they're in the right order again,
            dst_axes = np.arange(len(data.shape)-2).tolist()
            dst_axes.insert(angle_loc, len(data.shape)-2)
            dst_axes.insert(angle_loc+1, len(data.shape)-1)
            dst_axes.remove(angle_loc)
            dst_axes.append(angle_loc)
            data = np.ascontiguousarray(np.transpose(data, axes=dst_axes))

            # update the rest of the angle locations
            angle_locations += 1

        self.n = data.shape[-1]
        grid_aug = grid + [np.arange(data.shape[-1])]
        self.itp = RegularGridInterpolator(
            grid_aug, 
            data,
            bounds_error=False, 
            fill_value=None
        )

    def __call__(self, points):

        # Can't do any interpolation with 1 point, just return original data.
        if self.single_point_data is not None:
            return self.single_point_data

        x = np.zeros((self.n, len(points) + 1 +
                      np.sum(self.lut_interp_types != 'n')))
        offset_count = 0
        for i in range(len(points)):
            if self.lut_interp_types[i] == 'n':
                x[:, i + offset_count] = points[i]
            elif self.lut_interp_types[i] == 'r':
                x[:, i + offset_count] = np.cos(points[i])
                x[:, i + 1 + offset_count] = np.sin(points[i])
                offset_count += 1
            elif self.lut_interp_types[i] == 'd':
                x[:, i + offset_count] = np.cos(points[i] / 180. * np.pi)
                x[:, i + 1 + offset_count] = np.sin(points[i] / 180. * np.pi)
                offset_count += 1

        # This last dimension is always an integer so no
        # interpolation is performed. This is done only
        # for performance reasons.
        x[:, -1] = np.arange(self.n)
        res = self.itp(x)

        return res


class LUT:

  def __init__(self, wl, ref_wave, matfile=None):

    # Load the .mat file
    self.wl = wl
    D = loadmat(matfile)

    # Data includes correction based on 1) Wind 2) Solar Zen 3) Clouds
    data = D['data']

    # Parameter grid for 1) Wind 2) Solar Zen and 3) Clouds
    self.lut_grid = [grid[0] for grid in D['grids'][0]]

    # Interpolation type
    self.interp_types = np.array(['n' for n in self.lut_grid])

    # Wavelengths included in the lookup table. Could extend to longer wl
    wl_LUT = D['wl'][0]

    # Parameter grid for 1) Wind 2) Solar Zen and 3) Clouds
    self.grids = D['grids']

    # Not exactly sure what this variable is for
    self.band_names = [str(g).strip() for g in self.lut_grid[1:]]
    
    # resample all LUT spectra to your image bands
    old_shape = data.shape
    new_shape = list(old_shape[:])
    new_shape[-1] = len(self.wl)
    nb = int(data.shape[-1])
    data = data.reshape((int(data.size / nb), nb))
    data_resamp = []
    for x in data:
        use = x > 0
        x_resamp = interp1d(
            wl_LUT[use], 
            x[use], 
            bounds_error=False, 
            fill_value='extrapolate'
        )(wl)
        x_resamp[x_resamp<0] = 0
        data_resamp.append(x_resamp)
    self.data = np.array(data_resamp).reshape(new_shape)

    # We use the 'vector interpolator' object from ISOFIT
    self.data = VectorInterpolator(
        self.lut_grid, 
        self.data, 
        self.interp_types
    ) 

    # Get band for the subtraction
    self.ref_band = np.argmin(abs(ref_wave - wl))

    # Wavelengths for the correction
    self.use = np.logical_and(self.wl > wl_LUT[0], self.wl < wl_LUT[-1])


def optimize_windspeed_correction(hy_obj):
    """
    Glint correction algorithm following:

    Extended Hochberg correction to include sky glint correction
    This is the optimization step to find the effective wind speed
    """
    corr_wave = hy_obj.get_wave(hy_obj.glint['correction_wave'])
    
    # Get solar zenith array
    solar_zn = hy_obj.get_anc('solar_zn', radians=False)

    # Get LUT
    lut = LUT(
        hy_obj.wavelengths, 
        hy_obj.glint['correction_wave'],
        hy_obj.glint['lut']
    )

    # Initial gueses for wind optimization
    x0 = [10]
    bounds = hy_obj.glint['bounds']

    # Initialize correction
    correction = np.zeros([*corr_wave.shape, len(hy_obj.wavelengths)])
    windspeeds = np.zeros([*corr_wave.shape, len(hy_obj.wavelengths)])
    # Initialize iterator
    it = np.nditer(corr_wave, flags=['multi_index'])
    for pixel in it:
        y, x = it.multi_index

        if not hy_obj.mask['water'][y, x]:
            continue

        # Get solar zenith
        sol_zen = solar_zn[y, x]

        # Optimize wind speed
        res = minimize(
             err,
             x0,
             args=(float(pixel), lut, sol_zen),
             method='tnc',
             tol=1e-9,
             bounds=bounds
        )
        glint = glint_spectrum(lut, float(res.x), sol_zen)
        correction[y, x, :] = glint
        windspeeds[y, x, :] = res.x 

     return correction 


def apply_sky_sun_glint_correction(hy_obj,data,dimension,index):
    """
    Glint correction algorithm following a Hochberg + Sky Glint correction

    This applies the correction to the sample of data
    """

    if 'sky_glint_correction' not in hy_obj.ancillary:
        hy_obj.ancillary['sky_glint_correction'] = (
            optimize_windspeed_correction(hy_obj)
        )

    if 'water' not in hy_obj.mask:
        hy_obj.gen_mask(mask_create,'water',hy_obj.glint['calc_mask'])

    if dimension == 'line':
        correction = hy_obj.ancillary['sky_glint_correciton'][index, :]

    elif dimension == 'column':
        correction = hy_obj.ancillary['sky_glint_correciton'][:, index]

    elif dimension == 'band':
        correction = hy_obj.ancillary['sky_glint_correciton'][:, :, index]

    elif dimension == 'chunk':
        correction = hy_obj.ancillary['sky_glint_correciton'][y1:y2, x1:x2, :]

    elif dimension == 'pixels':
        y, x = index
        correction = hy_obj.ancillary['sky_glint_correciton'][y, x, :]

    data = data - correction

    return data


def glint_spectrum(lut, windspeed, solzen):
    """
    Given the windspeed and solar zenith angle, produce a glint 
    spectrum in reflectance units
    """
    cloudfrac = 0 
    vector = np.array([windspeed,solzen,cloudfrac])

    return lut.data(vector) 


def err(x, corr_val, lut, sol_zen):
    """
    windspeed fit error, assuming the channel "ref_band" is pure glint.
    """
    mdl = glint_spectrum(lut, x, sol_zen)
    er = pow(mdl[lut.ref_band] - corr_val, 2)

    return er
