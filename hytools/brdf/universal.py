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

This module contains functions to calculate and apply a single ('universal')
set of multiplicative BRDF correction coefficents. Coefficients can be calculated
per flightline or across multiple flightlines.

"""
from itertools import product
from copy import deepcopy
import numpy as np
import ray
from scipy.optimize import minimize
from .kernels import calc_volume_kernel,calc_geom_kernel
from ..misc import progbar
from ..misc import update_brdf
from ..masks import mask_create
from ..plotting import universal_diagno_plot

def universal_brdf(actors,config_dict):
    brdf_dict = config_dict['brdf']

    if brdf_dict['grouped']:
        actors = calc_universal_group(actors)
    else:
        _ = ray.get([a.do.remote(calc_universal_single) for a in actors])

    if brdf_dict['diagnostic_plots']:
        print('Exporting diagnostic plots.')
        _ = ray.get([a.do.remote(universal_diagno_plot,config_dict) for a in actors])

def sample_kernels(hy_obj):
    '''Calculate and sample BRDF kernels
    '''
    #Sample kernel images
    geom_kernel = hy_obj.geom_kernel(hy_obj.brdf['geometric'],
                                     b_r=hy_obj.brdf["b/r"],
                                     h_b =hy_obj.brdf["h/b"])[hy_obj.mask['calc_brdf']]
    vol_kernel = hy_obj.volume_kernel(hy_obj.brdf['volume'])[hy_obj.mask['calc_brdf']]
    X = np.vstack([vol_kernel,geom_kernel,
                   np.ones(vol_kernel.shape)]).T
    return X

def subsample_mask(hy_obj):
    '''Subsample and update calculation mask
    '''
    if hy_obj.brdf['sample_perc'] < 1:
        idx = np.array(np.where(hy_obj.mask['calc_brdf'])).T
        idx_rand= idx[np.random.choice(range(len(idx)),
                                      int(len(idx)*(1- hy_obj.brdf['sample_perc'])),
                                      replace = False)].T
        hy_obj.mask['calc_brdf'][idx_rand[0],idx_rand[1]] = False

def calc_universal_single(hy_obj):
    '''Calculate BRDF coefficients on a per flightline basis.
    '''
    subsample_mask(hy_obj)
    X = sample_kernels(hy_obj)

    hy_obj.brdf['coeffs'] = {}
    for band_num,band in enumerate(hy_obj.bad_bands):
        if ~band:
            band = hy_obj.get_band(band_num,
                                   corrections = hy_obj.corrections, mask='calc_brdf')
            brdf_coeff = np.linalg.lstsq(X, band,rcond=None)[0].flatten().tolist()
            hy_obj.brdf['coeffs'][band_num] = brdf_coeff

def calc_universal_group(actors):
    '''Calculate BRDF coefficients using pooled data from all flightlines.
    '''
    _ = ray.get([a.do.remote(subsample_mask) for a in actors])
    X = ray.get([a.do.remote(sample_kernels) for a in actors])
    X = np.concatenate(X)

    bad_bands = ray.get(actors[0].do.remote(lambda x: x.bad_bands))
    corections = ray.get(actors[0].do.remote(lambda x: x.corrections))

    coeffs = {}

    for band_num,band in enumerate(bad_bands):
        if ~band:
            y = ray.get([a.get_band.remote(band_num,mask='calc_brdf',
                                           corrections = corections) for a in actors])
            y = np.concatenate(y)
            coeffs[band_num] = np.linalg.lstsq(X, y)[0].flatten().tolist()
            progbar(np.sum(~bad_bands[:band_num+1]),np.sum(~bad_bands))
    print('\n')

    #Update BRDF coeffs
    _ = ray.get([a.do.remote(update_brdf,{'key':'coeffs',
                                          'value': coeffs}) for a in actors])

    return actors


def apply_universal(hy_obj,data,dimension,index):
    ''' Apply universal BRDF correction to a slice of the data

    Args:
        hy_obj : Hytools class object.
        data (np.ndarray): Data slice.
        index (int,list): Data index(ices).

    Returns:
        data (np.ndarray): BRDF correct data slice.
    '''

    if 'k_vol' not in hy_obj.ancillary:
        hy_obj.ancillary['k_vol'] = hy_obj.volume_kernel(hy_obj.brdf['volume'])
    if 'k_geom' not in hy_obj.ancillary:
        hy_obj.ancillary['k_geom'] = hy_obj.geom_kernel(hy_obj.brdf['geometric'],
                                     b_r=hy_obj.brdf["b/r"],
                                     h_b =hy_obj.brdf["h/b"])
    if ('k_vol_nadir' not in hy_obj.ancillary) or ('k_geom_nadir' not in hy_obj.ancillary):
        solar_zn = hy_obj.brdf['solar_zn_norm_radians']  * np.ones((hy_obj.lines,hy_obj.columns))
        hy_obj.ancillary['k_vol_nadir']  = calc_volume_kernel(0,solar_zn,
                                                               0,0,hy_obj.brdf['volume'])
        hy_obj.ancillary['k_geom_nadir']  = calc_geom_kernel(0,solar_zn,
                                                             0,0,hy_obj.brdf['geometric'],
                                                             b_r=hy_obj.brdf["b/r"],
                                                             h_b =hy_obj.brdf["h/b"])
    if 'apply_brdf' not in hy_obj.mask:
        hy_obj.gen_mask(mask_create,'apply_brdf',hy_obj.brdf['apply_mask'])

    brdf_bands = [int(x) for x in hy_obj.brdf['coeffs'].keys()]
    fvol, fgeo, fiso  = np.array([hy_obj.brdf['coeffs'][band] for band in hy_obj.brdf['coeffs'].keys()]).T

    #Convert to float
    data = data.astype(np.float32)

    if dimension == 'line':

        brdf = fvol[:,np.newaxis]*hy_obj.ancillary['k_vol'][[index],:]
        brdf+= fgeo[:,np.newaxis]*hy_obj.ancillary['k_geom'][[index],:]
        brdf+= fiso[:,np.newaxis]

        brdf_nadir = fvol[:,np.newaxis]*hy_obj.ancillary['k_vol_nadir'][[index],:]
        brdf_nadir+= fgeo[:,np.newaxis]*hy_obj.ancillary['k_geom_nadir'][[index],:]
        brdf_nadir+= fiso[:,np.newaxis]

        correction_factor = brdf_nadir/brdf
        correction_factor[:,~hy_obj.mask['apply_brdf'][index,:]] = 1
        data[:,brdf_bands] = data[:,brdf_bands]*correction_factor.T

    elif dimension == 'column':

        brdf = fvol[np.newaxis,:]*hy_obj.ancillary['k_vol'][:,[index]]
        brdf+= fgeo[np.newaxis,:]*hy_obj.ancillary['k_geom'][:,[index]]
        brdf+= fiso[np.newaxis,:]

        brdf_nadir = fvol[np.newaxis,:]*hy_obj.ancillary['k_vol_nadir'][:,[index]]
        brdf_nadir+= fgeo[np.newaxis,:]*hy_obj.ancillary['k_geom_nadir'][:,[index]]
        brdf_nadir+= fiso[np.newaxis,:]

        correction_factor = brdf_nadir/brdf
        correction_factor[~hy_obj.mask['apply_brdf'][:,index],:] = 1

        data[:,brdf_bands] = data[:,brdf_bands]*correction_factor.T

    elif dimension == 'band':
        fvol, fgeo, fiso  = hy_obj.brdf['coeffs'][index]
        brdf = fvol*hy_obj.ancillary['k_vol']
        brdf += fgeo*hy_obj.ancillary['k_geom']
        brdf+=fiso

        brdf_nadir = fvol*hy_obj.ancillary['k_vol_nadir']
        brdf_nadir+= fgeo*hy_obj.ancillary['k_geom_nadir']
        brdf_nadir+= fiso

        correction_factor = brdf_nadir/brdf
        correction_factor[~hy_obj.mask['apply_brdf']] = 1
        data= data* correction_factor

    elif dimension == 'chunk':
        x1,x2,y1,y2 = index

        brdf = fvol[np.newaxis,np.newaxis,:]*hy_obj.ancillary['k_vol'][y1:y2,x1:x2,np.newaxis]
        brdf+= fgeo[np.newaxis,np.newaxis,:]*hy_obj.ancillary['k_geom'][y1:y2,x1:x2,np.newaxis]
        brdf+= fiso[np.newaxis,np.newaxis,:]

        brdf_nadir = fvol[np.newaxis,np.newaxis,:]*hy_obj.ancillary['k_vol_nadir'][y1:y2,x1:x2,np.newaxis]
        brdf_nadir+= fgeo[np.newaxis,np.newaxis,:]*hy_obj.ancillary['k_geom_nadir'][y1:y2,x1:x2,np.newaxis]
        brdf_nadir+= fiso[np.newaxis,np.newaxis,:]

        correction_factor = brdf_nadir/brdf
        correction_factor[~hy_obj.mask['apply_brdf'][y1:y2,x1:x2]] = 1

        data[:,:,brdf_bands] = data[:,:,brdf_bands]*correction_factor

    elif dimension == 'pixels':
        y,x = index

        brdf = fvol[np.newaxis,:]*hy_obj.ancillary['k_vol'][y,x,np.newaxis]
        brdf+= fgeo[np.newaxis,:]*hy_obj.ancillary['k_geom'][y,x,np.newaxis]
        brdf+= fiso[np.newaxis,:]

        brdf_nadir = fvol[np.newaxis,:]*hy_obj.ancillary['k_vol_nadir'][y,x,np.newaxis]
        brdf_nadir+= fgeo[np.newaxis,:]*hy_obj.ancillary['k_geom_nadir'][y,x,np.newaxis]
        brdf_nadir+= fiso[np.newaxis,:]

        correction_factor = brdf_nadir/brdf
        correction_factor[~hy_obj.mask['apply_brdf'][y,x]] = 1

        data[:,brdf_bands] = data[:,brdf_bands]*correction_factor

    return data
