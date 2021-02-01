# -*- coding: utf-8 -*-
"""
This module contains functions to calculate and apply a single ('universal')
set of multiplicative BRDF correction coefficents. Coefficients can be calculated
per flightline or across multiple flightlines.
"""
from itertools import product
from copy import deepcopy
import numpy as np
import ray
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import minimize
from .kernels import calc_volume_kernel,calc_geom_kernel
from ..misc import progbar
from ..masks import mask_dict

diagnostic_waves = [450,850,1660,2200]

def universal_brdf(actors,config_dict):
    brdf_dict= config_dict['brdf']

    #Generate binary masks for coefficient calulation
    brdf_masker = mask_dict[brdf_dict['calc_mask'][0]]
    _ = ray.get([a.gen_mask.remote(brdf_masker,
                                   'calc_brdf',
                                   brdf_dict['calc_mask'][1]) for a in actors])

    #Set BRDF correction params
    _ = ray.get([a.do.remote(set_brdf_coeffs,brdf_dict) for a in actors])

    #Automatically determine optimal kernel
    if brdf_dict['auto_kernel'] == True:
        auto_kernel(actors,brdf_dict)
        #Rerun mask
        _ = ray.get([a.gen_mask.remote(brdf_masker,
                                       'calc_brdf',
                                       brdf_dict['calc_mask'][1]) for a in actors])

    if brdf_dict['grouped']:
        actors = calc_universal_group(actors,brdf_dict)
    else:
        _ = ray.get([a.do.remote(calc_universal_single,brdf_dict) for a in actors])

    if brdf_dict['diagnostic_plots'] == True:
        print('Exporting diagnostic plots.')
        _ = ray.get([a.do.remote(brdf_diagnostic_plot,config_dict) for a in actors])


def set_brdf_coeffs(hy_obj,brdf_dict):
    hy_obj.brdf = brdf_dict

def sample_kernels(hy_obj):
    '''Calculate and sample BRDF kernels
    '''
    hy_obj.brdf['coeffs'] = {}

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

def auto_kernel(actors,brdf_dict):
    '''Automatically determine the optimal kernel combination and parameters for
    BRDF correction.
    '''
    print('EXPERIMENTAL FEATURE!')
    print('Determining optimal kernel combination and parameters.')

    geometric = ['li_sparse_r','li_dense_r','roujean']
    volume = ['ross_thin','ross_thick','hotspot']

    orig_sample = deepcopy(brdf_dict['sample_perc'])

    #Update BRDF correction params
    brdf_dict['sample_perc'] = 0.01

    _ = ray.get([a.do.remote(set_brdf_coeffs,brdf_dict) for a in actors])

    #Create subsampling mask
    _ = ray.get([a.do.remote(subsample_mask) for a in actors])

    corections = ray.get(actors[0].do.remote(lambda x: x.corrections))

    # Load viewing geometry to memory
    solar_zn = ray.get([a.get_anc.remote('solar_zn',mask='calc_brdf') for a in actors])
    solar_zn = np.concatenate(solar_zn)
    solar_az = ray.get([a.get_anc.remote('solar_az',mask='calc_brdf') for a in actors])
    solar_az = np.concatenate(solar_az)
    sensor_zn = ray.get([a.get_anc.remote('sensor_zn',mask='calc_brdf') for a in actors])
    sensor_zn = np.concatenate(sensor_zn)
    sensor_az = ray.get([a.get_anc.remote('sensor_az',mask='calc_brdf') for a in actors])
    sensor_az = np.concatenate(sensor_az)

    band_samples = []
    diagnostic_bands = [ray.get(actors[0].wave_to_band.remote(wave)) for wave in diagnostic_waves]

    for band in diagnostic_bands:
        y = ray.get([a.get_band.remote(band,mask='calc_brdf',
                               corrections = corections) for a in actors])
        band_samples.append(np.concatenate(y))

    #Use only 500,000 samples
    if len(sensor_az) > 500000:
        regress_mask = np.zeros(len(sensor_az)).astype(bool)
        subset_indices = np.random.choice(np.arange(len(sensor_az)),
                                        size=500000,replace= False)
        regress_mask[subset_indices] = True
    else:
        regress_mask = np.ones(len(sensor_az)).astype(bool)

    def kernel_opt(args):
        '''Kernel minimization function. Returns the average RMSE across
        five diagnostic bands.

        '''
        b_r,h_b = args

        if geom:
            geom_kernel= calc_geom_kernel(solar_az,solar_zn,sensor_az,
                                          sensor_zn,geom,b_r=b_r,h_b =h_b)
            X = np.vstack([geom_kernel,np.ones(geom_kernel.shape)]).T
        if vol:
            vol_kernel = calc_volume_kernel(solar_az,solar_zn,
                                            sensor_az,sensor_zn,vol)
            X = np.vstack([vol_kernel,np.ones(vol_kernel.shape)]).T
        if vol and geom:
            X = np.vstack([vol_kernel,geom_kernel,
                       np.ones(vol_kernel.shape)]).T
        elif (not isinstance(vol,str)) & (not isinstance(geom,str)):
            return 1E11

        band_rmse = []
        for y in band_samples:
            coeffs = np.linalg.lstsq(X[regress_mask],
                                     y[regress_mask])[0].flatten().tolist()
            pred = (X[regress_mask]*coeffs).sum(axis=1)
            obs = y[regress_mask]
            band_rmse.append(np.sqrt(np.mean((obs-pred)**2)))
        return np.mean(band_rmse)

    min_rmse = 1E10

    #Cycle through all combinations of kernels
    for vol, geom in product(volume, geometric):

        result = minimize(kernel_opt, (1,2),
                       method='Nelder-Mead',
                       tol=1e-6,
                       bounds = ((0.25,20),(0.25,20)))

        if isinstance(vol,str) | isinstance(geom,str):
            print("Kernel combination: %s, %s" % (vol,geom))
            print("Object shapes: %s, %s" % (round(result.x[0],2),round(result.x[1],2)))
            rmse = kernel_opt(result.x)
            print("RMSE: %s" % round(rmse,8))
            print("\n##############################")

            if (rmse < min_rmse) & ((result.x > 0).sum() == 2):
                min_rmse= rmse
                opt_geom = geom
                opt_vol = vol
                opt_b_r = result.x[0]
                opt_h_b = result.x[1]

    print("Optimal kernel combination: %s, %s" % (opt_vol,opt_geom))
    #Update BRDF correction params
    brdf_dict['sample_perc'] = orig_sample
    brdf_dict['volume']  = opt_vol
    brdf_dict['geometric'] = opt_geom
    brdf_dict['b/r'] =opt_b_r
    brdf_dict['h/b'] = opt_h_b
    _ = ray.get([a.do.remote(set_brdf_coeffs,brdf_dict) for a in actors])

def calc_universal_single(hy_obj):
    '''Calculate BRDF coefficients on a per flightline basis.
    '''
    subsample_mask(hy_obj)
    X = sample_kernels(hy_obj)

    for band_num,band in enumerate(hy_obj.bad_bands):
        if ~band:
            band = hy_obj.get_band(band_num,
                                   corrections = hy_obj.corrections, mask='calc_brdf')
            brdf_coeff = np.linalg.lstsq(X, band,rcond=None)[0].flatten().tolist()
            hy_obj.brdf['coeffs'][band_num] = brdf_coeff

def calc_universal_group(actors,brdf_dict):
    '''Calculate BRDF coefficients using pooled data from all flightlines.
    '''
    _ = ray.get([a.do.remote(subsample_mask) for a in actors])
    X = ray.get([a.do.remote(sample_kernels) for a in actors])
    X = np.concatenate(X)

    bad_bands = ray.get(actors[0].do.remote(lambda x: x.bad_bands))
    corections = ray.get(actors[0].do.remote(lambda x: x.corrections))

    brdf_dict['coeffs'] = {}

    for band_num,band in enumerate(bad_bands):
        if ~band:
            y = ray.get([a.get_band.remote(band_num,mask='calc_brdf',
                                           corrections = corections) for a in actors])
            y = np.concatenate(y)
            coeffs = np.linalg.lstsq(X, y)[0].flatten().tolist()
            brdf_dict['coeffs'][band_num] = coeffs
            progbar(np.sum(~bad_bands[:band_num+1]),np.sum(~bad_bands))
    print('\n')

    #Set BRDF coefficients
    _ = ray.get([a.do.remote(set_brdf_coeffs,brdf_dict) for a in actors])

    return actors


def brdf_diagnostic_plot(hy_obj,config_dict):
    ''' Generate a diagnostic plot of BRDF correction results.
    '''
    #Flip sign of zenith angle
    sensor_zn = hy_obj.get_anc('sensor_zn',radians =False)
    sensor_zn[~hy_obj.mask['no_data']] = np.nan
    for i,line in enumerate(sensor_zn):
        line[:np.nanargmin(line)] *= -1
        sensor_zn[i] = line
    sensor_zn = (sensor_zn[hy_obj.mask['calc_brdf']]//2)*2

    diagno_df = pd.DataFrame()
    diagno_df['sensor_zn'] =sensor_zn

    diagnostic_bands = [hy_obj.wave_to_band(wave) for wave in diagnostic_waves]
    for band_num in diagnostic_bands:
        band = hy_obj.get_band(band_num,mask='calc_brdf')
        diagno_df['uncorr_%s' % band_num] =  band
        band = hy_obj.get_band(band_num,
                               corrections = hy_obj.corrections,
                               mask='calc_brdf')
        diagno_df['corr_%s' % band_num] =  band

        fvol, fgeo, fiso  = hy_obj.brdf['coeffs'][band_num]
        brdf = fvol*hy_obj.ancillary['k_vol']
        brdf += fgeo*hy_obj.ancillary['k_geom']
        brdf+=fiso
        brdf = brdf[hy_obj.mask['calc_brdf']]
        diagno_df['brdf_%s' % band_num] =  brdf

    #Average every 2 degrees of zenith angle
    diagno_df=  diagno_df.groupby(by= 'sensor_zn').mean()

    fig = plt.figure(figsize= (8,6))
    for a,band_num in enumerate(diagnostic_bands,start=1):
        ax = fig.add_subplot(2,2,a)
        ax.plot(diagno_df.index,diagno_df['brdf_%s' % band_num],c='k',ls ='--')
        ax.scatter(diagno_df.index,diagno_df['uncorr_%s' % band_num],marker ='o',fc='w',ec='k')
        ax.scatter(diagno_df.index,diagno_df['corr_%s' % band_num],marker ='o',fc='k',ec='k')
        ax.text(.85,.9, "%s nm" % int(hy_obj.wavelengths[band_num]), transform=ax.transAxes,
                ha = 'center', fontsize = 12)
        if a > 2:
            ax.set_xlabel('View zenith angle')
        if a in [1,3]:
            ax.set_ylabel('Reflectance')

    #Create legend
    custom_points = []
    custom_points.append(Line2D([0],[0], marker = 'o',label='Uncorrected',
                          markerfacecolor='w', markersize=10,lw=0,markeredgecolor='k'))
    custom_points.append(Line2D([0],[0], marker = 'o',label='Corrected',
                          markerfacecolor='k', markersize=10,lw=0,markeredgecolor='k'))
    custom_points.append(Line2D([0],[1],label='Modeled BRDF',c='k', ls ='--'))
    ax.legend(handles=custom_points, loc='center',frameon=False,
              bbox_to_anchor=(-.15, -.3), ncol =3,columnspacing = 1.5,labelspacing=.25)

    plt.savefig("%s%s_brdf_plot.png" % (config_dict['export']['output_dir'],hy_obj.base_name))
    plt.close()

def apply_universal(hy_obj,data,dimension,index):
    ''' Apply SCSS correction to a slice of the data

    Args:
        hy_obj : Hytools class object.
        data (np.ndarray): Data slice.
        index (int,list): Data index(ices).

    Returns:
        data (np.ndarray): BRDF correct data slice.
    '''

    # Load kernels and mask to memory if not already there
    if 'k_vol' not in hy_obj.ancillary.keys():
        hy_obj.ancillary['k_vol'] = hy_obj.volume_kernel(hy_obj.brdf['volume'])

    if 'k_geom' not in hy_obj.ancillary.keys():
        hy_obj.ancillary['k_geom'] = hy_obj.geom_kernel(hy_obj.brdf['geometric'],
                                    b_r=hy_obj.brdf["b/r"],
                                    h_b =hy_obj.brdf["h/b"])

    if hy_obj.brdf['solar_zn_norm']:
        solar_zn = hy_obj.brdf['mean_solar_zenith']  * np.ones((hy_obj.lines,hy_obj.columns))
    else:
        solar_zn = hy_obj.get_anc('solar_zn')

    if 'k_vol_nadir' not in hy_obj.ancillary.keys():
        hy_obj.ancillary['k_vol_nadir'] = calc_volume_kernel(0,solar_zn,
                                     0,0,hy_obj.brdf['volume'])

    if 'k_geom_nadir' not in hy_obj.ancillary.keys():
        hy_obj.ancillary['k_geom_nadir'] = calc_geom_kernel(0,solar_zn,
                                        0,0,hy_obj.brdf['geometric'],
                                        b_r=hy_obj.brdf["b/r"],
                                        h_b =hy_obj.brdf["h/b"])

    if 'apply_brdf' not in hy_obj.mask:
        brdf_masker = mask_dict[hy_obj.brdf['apply_mask'][0]]
        hy_obj.gen_mask(brdf_masker,
                        'apply_brdf',
                        hy_obj.brdf['apply_mask'][1])

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

        data[brdf_bands,:] = data[brdf_bands,:]*correction_factor

    elif dimension == 'column':

        brdf = fvol[np.newaxis,:]*hy_obj.ancillary['k_vol'][:,[index]]
        brdf+= fgeo[np.newaxis,:]*hy_obj.ancillary['k_geom'][:,[index]]
        brdf+= fiso[np.newaxis,:]

        brdf_nadir = fvol[np.newaxis,:]*hy_obj.ancillary['k_vol_nadir'][:,[index]]
        brdf_nadir+= fgeo[np.newaxis,:]*hy_obj.ancillary['k_geom_nadir'][:,[index]]
        brdf_nadir+= fiso[np.newaxis,:]

        correction_factor = brdf_nadir/brdf
        correction_factor[~hy_obj.mask['apply_brdf'][:,index],:] = 1

        data[:,brdf_bands] = data[:,brdf_bands]*correction_factor

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
