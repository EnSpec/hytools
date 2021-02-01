# -*- coding: utf-8 -*-
"""
This module contains functions apply an empirical BRDF correction described
the following paper:

Equations and constants can be found in the following papers:

"""

from itertools import tee
import numpy as np
import ray
from scipy.interpolate import interp1d
from .kernels import calc_volume_kernel,calc_geom_kernel
from ..masks import neon_edge
from ..misc import progbar

def flex_brdf(actors,config_dict):
    brdf_dict= config_dict['brdf']
    if brdf_dict['grouped']:
        actors = calc_flex_group(actors,brdf_dict)
    else:
        _ = ray.get([a.do.remote(calc_flex_single,brdf_dict) for a in actors])

def ndvi_class(hy_obj):
    '''
    Create NDVI bin class mask
    '''
    ndvi = hy_obj.ndi()
    k_vol  = hy_obj.volume_kernel(hy_obj.brdf['volume'])
    k_geom = hy_obj.geom_kernel(hy_obj.brdf['geometric'],
                                b_r=hy_obj.brdf["b/r"],
                                h_b =hy_obj.brdf["h/b"])
    k_finite = np.isfinite(k_vol) & np.isfinite(k_geom)

    class_mask = np.zeros((hy_obj.lines, hy_obj.columns))

    for bin_num in hy_obj.brdf['bins']:
        start,end =  hy_obj.brdf['bins'][bin_num]
        class_mask[(ndvi > start) & (ndvi <= end)] = bin_num

    #Exclude from analysis
    class_mask[hy_obj.get_anc('sensor_zn') < hy_obj.brdf['sensor_zn_min']] = 0
    class_mask[~hy_obj.mask['no_data']] = 0
    class_mask[~k_finite] = 0

    #Exclude edge pixels in NEON imagery
    if hy_obj.brdf['neon_buffer']:
        class_mask[neon_edge(hy_obj)] = 0

    #Subsample data
    idx = np.array(np.where(class_mask!=0)).T
    idxRand= idx[np.random.choice(range(len(idx)),int(len(idx)*(1-hy_obj.brdf['sample_perc'])), replace = False)].T
    class_mask[idxRand[0],idxRand[1]] = 0
    class_mask = class_mask.astype(np.int8)
    hy_obj.ancillary['ndvi_classes'] = class_mask

def ndvi_bins(ndvi,brdf_dict):
    '''Calculate NDVI bin extents
    '''
    perc_range = brdf_dict['ndvi_perc_max'] - brdf_dict['ndvi_perc_min'] + 1

    ndvi_break_dyn_bin = np.percentile(ndvi[ndvi > 0],
                                       np.arange(brdf_dict['ndvi_perc_min'],
                                                brdf_dict['ndvi_perc_max'] + 1,
                                                perc_range / (brdf_dict['num_bins'] - 1)))
    ndvi_thres = [brdf_dict['ndvi_bin_min']]
    ndvi_thres += ndvi_break_dyn_bin.tolist()
    ndvi_thres += [brdf_dict['ndvi_bin_max']]
    ndvi_thres = sorted(list(set(ndvi_thres)))
    bins = [[x,y] for x,y in pairwise(ndvi_thres)]
    return bins

def pairwise(iterable):
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def set_brdf_coeffs(hy_obj,brdf_coeffs):
    hy_obj.brdf  = brdf_coeffs

def get_kernel_samples(hy_obj):
    '''Calculate and sample BRDF kernels
    '''
    geom_kernel = hy_obj.geom_kernel(hy_obj.brdf['geometric'],
                                     b_r=hy_obj.brdf["b/r"] ,
                                     h_b =hy_obj.brdf["h/b"])
    geom_kernel = geom_kernel[hy_obj.ancillary['ndvi_classes'] !=0]

    vol_kernel = hy_obj.volume_kernel(hy_obj.brdf['volume'])
    vol_kernel = vol_kernel[hy_obj.ancillary['ndvi_classes'] !=0]

    classes = hy_obj.ancillary['ndvi_classes'][hy_obj.ancillary['ndvi_classes'] !=0]
    X = np.vstack([vol_kernel,geom_kernel,
                   np.ones(vol_kernel.shape),classes]).T
    return X

def get_band_samples(hy_obj,args):
    band = hy_obj.get_band(args['band_num'],
                           corrections = hy_obj.corrections)
    return band[hy_obj.ancillary['ndvi_classes'] !=0]

def calc_flex_single(hy_obj,brdf_dict):
    ''' Calculate BRDF coefficents for a single image
    '''
    hy_obj.brdf = brdf_dict
    hy_obj.brdf['coeffs'] ={}

    # Determine bin dimensions and create class mask
    if hy_obj.brdf['bin_type'] == 'dynamic':
        bins = ndvi_bins(hy_obj.ndi()[~hy_obj.mask['no_data']],brdf_dict)
    else:
        bins = brdf_dict['bins']

    hy_obj.brdf['bins'] = {k:v for (k,v) in enumerate(bins,start=1)}

    ndvi_class(hy_obj)
    kernel_samples= get_kernel_samples(hy_obj)

    # Calculate coefficients for each band and class
    for band_num,band in enumerate(hy_obj.bad_bands):
        if ~band:
            hy_obj.brdf['coeffs'][band_num] = {}
            band_samples = hy_obj.do(get_band_samples, {'band_num':band_num})
            coeffs= []

            for bin_num in hy_obj.brdf['bins']:
                bin_mask = [kernel_samples[:,3] == bin_num]
                X = kernel_samples[:,:3][bin_mask]
                y = band_samples[bin_mask]
                coeffs.append(np.linalg.lstsq(X, y,rcond=-1)[0].flatten().tolist())
            hy_obj.brdf['coeffs'][band_num]  = coeffs

def calc_flex_group(actors,brdf_dict):
    ''' Calculate BRDF coefficents for a group of images
    '''
    # Aggregate NDVI values from images
    ndvi = ray.get([a.ndi.remote(mask = 'no_data') for a in actors])
    ndvi = np.concatenate([n.flatten() for n in ndvi])

    # Determine bin dimensions
    if  brdf_dict['bin_type'] == 'dynamic':
        bins = ndvi_bins(ndvi,brdf_dict)
    else:
        bins = brdf_dict['bins']
    brdf_dict['bins']  = {k:v for (k,v) in enumerate(bins,start=1)}

    #Update BRDF coeffs
    _ = ray.get([a.do.remote(set_brdf_coeffs,brdf_dict) for a in actors])

    #Create NDVI class mask and sample kernels
    _ = ray.get([a.do.remote(ndvi_class) for a in actors])
    kernel_samples = ray.get([a.do.remote(get_kernel_samples) for a in actors])
    kernel_samples = np.concatenate(kernel_samples)

    bad_bands = ray.get(actors[0].do.remote(lambda x: x.bad_bands))
    brdf_dict['coeffs'] = {}

    for band_num,band in enumerate(bad_bands):
        if ~band:
            brdf_dict['coeffs'][band_num] = {}
            band_samples = ray.get([a.do.remote(get_band_samples,
                                     {'band_num':band_num}) for a in actors])
            band_samples = np.concatenate(band_samples)
            coeffs= []
            for bin_num in brdf_dict['bins']:
                bin_mask = [kernel_samples[:,3] == bin_num]
                X = kernel_samples[:,:3][bin_mask]
                y = band_samples[bin_mask]
                coeffs.append(np.linalg.lstsq(X, y,rcond=-1)[0].flatten().tolist())
            brdf_dict['coeffs'][band_num]  = coeffs
            progbar(np.sum(~bad_bands[:band_num+1]),np.sum(~bad_bands))

    print('\n')
    #Update BRDF coeffs
    _ = ray.get([a.do.remote(set_brdf_coeffs,brdf_dict) for a in actors])


def apply_flex(hy_obj,data,dimension,index):
    ''' Apply flex BRDF correction to a slice of the data

    Args:
        hy_obj : Hytools class object.
        data (np.ndarray): Data slice.
        index (int,list): Data index.

    Returns:
        data (np.ndarray): BRDF corrected data slice.
    '''

    # Load kernels and NDVI to memory if not already there
    if 'k_vol' not in hy_obj.ancillary:
        hy_obj.ancillary['k_vol'] = hy_obj.volume_kernel(hy_obj.brdf['volume'])

    if 'k_geom' not in hy_obj.ancillary:
        hy_obj.ancillary['k_geom'] = hy_obj.geom_kernel(hy_obj.brdf['geometric'],
                                     b_r=hy_obj.brdf["b/r"],
                                     h_b =hy_obj.brdf["h/b"])

    if hy_obj.brdf['solar_zn_norm']:
        solar_zn = hy_obj.brdf['mean_solar_zenith']  * np.ones((hy_obj.lines,hy_obj.columns))
    else:
        solar_zn = hy_obj.get_anc('solar_zn')

    if 'k_vol_nadir' not in hy_obj.ancillary:
        hy_obj.ancillary['k_vol_nadir']  = calc_volume_kernel(0,solar_zn,
                                                               0,0,hy_obj.brdf['volume'])
    if 'k_geom_nadir' not in hy_obj.ancillary:
        hy_obj.ancillary['k_geom_nadir']  = calc_geom_kernel(0,solar_zn,
                                                             0,0,hy_obj.brdf['geometric'],
                                                             b_r=hy_obj.brdf["b/r"],
                                                             h_b =hy_obj.brdf["h/b"])
    if 'ndvi' not in hy_obj.ancillary:
        hy_obj.ancillary['ndvi'] =  hy_obj.ndi()

    if 'brdf' not in hy_obj.mask:
        hy_obj.mask['brdf']  = (hy_obj.ancillary['ndvi'] > hy_obj.brdf['ndvi_min']) & (hy_obj.ancillary['ndvi'] < hy_obj.brdf['ndvi_max'])

    if 'interpolators' not in hy_obj.brdf:
        bin_centers = np.mean(list(hy_obj.brdf['bins'].values()),axis=1)
        hy_obj.ancillary['interpolators'] ={}

        #Generate interpolators
        for i in hy_obj.brdf['coeffs']:
            coeffs= np.array(hy_obj.brdf['coeffs'][i])
            interpolator = interp1d(bin_centers, coeffs, kind = hy_obj.brdf['interp_kind'],
                                    axis=0,fill_value="extrapolate")
            hy_obj.ancillary['interpolators'][int(i)] = interpolator

    #Convert to float
    data = data.astype(np.float32)
    brdf_bands = [int(x) for x in hy_obj.ancillary['interpolators']]
    if dimension == 'line':
        # index= 3000
        # data = hy_obj.get_line(3000)

        interpolated_f = [hy_obj.ancillary['interpolators'][band](hy_obj.ancillary['ndvi'][index,:]) for band in brdf_bands]
        interpolated_f = np.array(interpolated_f)
        fvol, fgeo, fiso  = interpolated_f[:,:,0], interpolated_f[:,:,1], interpolated_f[:,:,2]

        brdf = fvol*hy_obj.ancillary['k_vol'][index,:]
        brdf+= fgeo*hy_obj.ancillary['k_geom'][index,:]
        brdf+= fiso

        brdf_nadir = fvol*hy_obj.ancillary['k_vol_nadir'][index,:]
        brdf_nadir+= fgeo*hy_obj.ancillary['k_geom_nadir'][index,:]
        brdf_nadir+= fiso

        correction_factor = brdf_nadir/brdf
        correction_factor[:,~hy_obj.mask['brdf'][index] == 0] = 1
        correction_factor = np.moveaxis(correction_factor,0,1)

        data[:,brdf_bands] = data[:,brdf_bands]*correction_factor

    elif dimension == 'column':
        #index= 300
        #data = hy_obj.get_column(index)

        interpolated_f = [hy_obj.ancillary['interpolators'][band](hy_obj.ancillary['ndvi'][:,index]) for band in brdf_bands]
        interpolated_f = np.array(interpolated_f)
        fvol, fgeo, fiso  = interpolated_f[:,:,0], interpolated_f[:,:,1], interpolated_f[:,:,2]

        brdf = fvol*hy_obj.ancillary['k_vol'][:,index]
        brdf+= fgeo*hy_obj.ancillary['k_geom'][:,index]
        brdf+= fiso

        brdf_nadir = fvol*hy_obj.ancillary['k_vol_nadir'][:,index]
        brdf_nadir+= fgeo*hy_obj.ancillary['k_geom_nadir'][:,index]
        brdf_nadir+= fiso

        correction_factor = brdf_nadir/brdf
        correction_factor = np.moveaxis(correction_factor,0,1)

        data[:,brdf_bands] = data[:,brdf_bands]

    elif (dimension == 'band') & (index in brdf_bands):
        # index= 8
        # data = hy_obj.get_band(index)

        interpolated_f = hy_obj.ancillary['interpolators'][index](hy_obj.ancillary['ndvi'])
        fvol, fgeo, fiso  = interpolated_f[:,:,0], interpolated_f[:,:,1], interpolated_f[:,:,2]

        brdf = fvol*hy_obj.ancillary['k_vol']
        brdf += fgeo*hy_obj.ancillary['k_geom']
        brdf += fiso

        brdf_nadir = fvol*hy_obj.ancillary['k_vol_nadir']
        brdf_nadir += fgeo*hy_obj.ancillary['k_geom_nadir']
        brdf_nadir += fiso

        correction_factor = brdf_nadir/brdf
        correction_factor[~hy_obj.mask['brdf']] = 1
        data= data* correction_factor

    elif dimension == 'chunk':
        # index = 200,501,3000,3501
        x1,x2,y1,y2 = index
        # data = hy_obj.get_chunk(x1,x2,y1,y2)

        interpolated_f = [hy_obj.ancillary['interpolators'][band](hy_obj.ancillary['ndvi'][y1:y2,x1:x2]) for band in brdf_bands]
        interpolated_f = np.array(interpolated_f)
        interpolated_f = np.swapaxes(interpolated_f,0,-1)
        fvol, fgeo, fiso  = interpolated_f[0,:,:,:], interpolated_f[1,:,:,:], interpolated_f[2,:,:,:]

        brdf = fvol*hy_obj.ancillary['k_vol'][y1:y2,x1:x2,np.newaxis]
        brdf+= fgeo*hy_obj.ancillary['k_geom'][y1:y2,x1:x2,np.newaxis]
        brdf+= fiso

        brdf_nadir = fvol*hy_obj.ancillary['k_vol_nadir'][y1:y2,x1:x2,np.newaxis]
        brdf_nadir+= fgeo*hy_obj.ancillary['k_geom_nadir'][y1:y2,x1:x2,np.newaxis]
        brdf_nadir+= fiso

        correction_factor = brdf_nadir/brdf
        correction_factor[~hy_obj.mask['brdf'][y1:y2,x1:x2] == 0] = 1
        data[:,:,brdf_bands] = data[:,:,brdf_bands]*correction_factor

    elif dimension == 'pixels':
        # index = [[2000,2001],[200,501]]
        y,x = index
        # data = hy_obj.get_pixels(y,x)

        interpolated_f = [hy_obj.ancillary['interpolators'][band](hy_obj.ancillary['ndvi'][y,x]) for band in brdf_bands]
        interpolated_f = np.array(interpolated_f)
        interpolated_f = np.swapaxes(interpolated_f,0,1)
        fvol, fgeo, fiso  = interpolated_f[:,:,0], interpolated_f[:,:,1], interpolated_f[:,:,2]

        brdf = fvol*hy_obj.ancillary['k_vol'][y,x,np.newaxis]
        brdf+= fgeo*hy_obj.ancillary['k_geom'][y,x,np.newaxis]
        brdf+= fiso

        brdf_nadir = fvol*hy_obj.ancillary['k_vol_nadir'][y,x,np.newaxis]
        brdf_nadir+= fgeo*hy_obj.ancillary['k_geom_nadir'][y,x,np.newaxis]
        brdf_nadir+= fiso

        correction_factor = brdf_nadir/brdf
        correction_factor[~hy_obj.mask['brdf'][y,x] == 0] = 1
        data[:,brdf_bands] = data[:,brdf_bands]*correction_factor

    return data
