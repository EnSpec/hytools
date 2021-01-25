# -*- coding: utf-8 -*-
"""
This module contains functions apply multiplicative BRDF correction:
"""
import numpy as np
import ray
from .kernels import calc_volume_kernel,calc_geom_kernel
from ..misc import progbar
from ..masks import mask_dict

def global_brdf(actors,brdf_dict):
    
    brdf_masker = mask_dict[brdf_dict['mask']]
    _ = ray.get([a.gen_mask.remote(brdf_masker,'brdf') for a in actors])
    
    if brdf_dict['grouped']:
        actors = calc_global_group(actors,brdf_dict)
    else:
        _ = ray.get([a.do.remote(calc_global_single,brdf_dict) for a in actors])


def calc_global_single(hy_obj,brdf_dict):
    '''
        Calculate BRDF coefficients on a per flightline basis.
    '''

    volume = brdf_dict['volume'] 
    geometric = brdf_dict['geometric']
    brdf_dict['coeffs'] = {}

    geom_kernel = hy_obj.geom_kernel(geometric,b_r=brdf_dict["b/r"] ,
                                     h_b =brdf_dict["h/b"])[hy_obj.mask['brdf']]
    vol_kernel = hy_obj.volume_kernel(volume)[hy_obj.mask['brdf']]
    X = np.vstack([vol_kernel,geom_kernel,np.ones(geom_kernel.shape)]).T

    for band_num,band in enumerate(hy_obj.bad_bands):
        if ~band:
            band = hy_obj.get_band(band_num,
                                   corrections = hy_obj.corrections, mask='brdf')
            brdf_coeff = np.linalg.lstsq(X, band,rcond=None)[0].flatten().tolist()
            brdf_dict['coeffs'][band_num] = brdf_coeff            
            
    hy_obj.brdf = brdf_dict
    hy_obj.corrections.append('brdf')
    
def calc_global_group(actors,brdf_dict):
    '''
        Calculate BRDF coefficients using pooled data from all flightlines.
    '''
    
    volume = brdf_dict['volume'] 
    geometric = brdf_dict['geometric']
    brdf_dict['coeffs'] = {}
    
    def get_kernel_samples(hy_obj):
        geom_kernel = hy_obj.geom_kernel(geometric,b_r=brdf_dict["b/r"] ,
                                     h_b =brdf_dict["h/b"])[hy_obj.mask['brdf']]
        vol_kernel = hy_obj.volume_kernel(volume)[hy_obj.mask['brdf']]
        X = np.vstack([vol_kernel,geom_kernel,
                       np.ones(vol_kernel.shape)]).T
        return X
    
    def set_brdf_coeffs(hy_obj,brdf_coeffs):
        hy_obj.brdf  = brdf_coeffs
        hy_obj.corrections.append('brdf')

    def sample_mask(hy_obj):    
        sub_samples = np.zeros((hy_obj.lines,hy_obj.columns)).astype(bool)
        idx = np.array(np.where(hy_obj.mask['brdf'])).T
        idxRand= idx[np.random.choice(range(len(idx)),int(len(idx)*brdf_dict['sample_perc']), replace = False)].T
        sub_samples[idxRand[0],idxRand[1]] = True
        return sub_samples    

    _ = ray.get([a.gen_mask.remote(sample_mask,'brdf_samples') for a in actors])
    X = ray.get([a.do.remote(get_kernel_samples) for a in actors])
    X = np.concatenate(X)
    
    bad_bands = ray.get(actors[0].do.remote(lambda x: x.bad_bands))
    corections = ray.get(actors[0].do.remote(lambda x: x.corrections))

    for band_num,band in enumerate(bad_bands):
        if ~band:
            y = ray.get([a.get_band.remote(band_num,mask='brdf_samples',
                                           corrections = corections) for a in actors])
            y = np.concatenate(y)
            coeffs = np.linalg.lstsq(X, y)[0].flatten().tolist()
            brdf_dict['coeffs'] [band_num] = coeffs
            progbar(np.sum(~bad_bands[:band_num+1]),np.sum(~bad_bands))

    print('\n')
    _ = ray.get([a.do.remote(set_brdf_coeffs,brdf_dict) for a in actors])
    
    return actors

def apply_standard(hy_obj,data,dimension,index):
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
        k_vol = hy_obj.volume_kernel(hy_obj.brdf['volume'])
        hy_obj.ancillary['k_vol'] = k_vol

    if 'k_geom' not in hy_obj.ancillary.keys():
        k_geom = hy_obj.geom_kernel(hy_obj.brdf['geometric'],
                                    b_r=hy_obj.brdf["b/r"],
                                    h_b =hy_obj.brdf["h/b"])
        hy_obj.ancillary['k_geom'] = k_geom

    if 'k_vol_nadir' not in hy_obj.ancillary.keys():
        k_vol_nadir = calc_volume_kernel(0,hy_obj.get_anc('solar_zn'),
                                     0,0,hy_obj.brdf['volume'])
        hy_obj.ancillary['k_vol_nadir'] = k_vol_nadir

    if 'k_geom_nadir' not in hy_obj.ancillary.keys():
        k_geom_nadir = calc_geom_kernel(0,hy_obj.get_anc('solar_zn'),
                                        0,0,hy_obj.brdf['geometric'],
                                        b_r=hy_obj.brdf["b/r"],
                                        h_b =hy_obj.brdf["h/b"])
        hy_obj.ancillary['k_geom_nadir'] = k_geom_nadir

    brdf_bands = [int(x) for x in hy_obj.brdf['coeffs'].keys()]
    fvol, fgeo, fiso  = np.array([hy_obj.brdf['coeffs'][band] for band in hy_obj.brdf['coeffs'].keys()]).T

    #Convert to float
    data = data.astype(np.float32)

    if dimension == 'line':
        #index= 3000
        #data = hy_obj.get_line(3000)

        brdf = fvol[:,np.newaxis]*hy_obj.ancillary['k_vol'][[index],:]
        brdf+= fgeo[:,np.newaxis]*hy_obj.ancillary['k_geom'][[index],:]
        brdf+= fiso[:,np.newaxis]

        brdf_nadir = fvol[:,np.newaxis]*hy_obj.ancillary['k_vol_nadir'][[index],:]
        brdf_nadir+= fgeo[:,np.newaxis]*hy_obj.ancillary['k_geom_nadir'][[index],:]
        brdf_nadir+= fiso[:,np.newaxis]

        correction_factor = brdf_nadir/brdf
        correction_factor[:,hy_obj.mask['brdf'][index,:]] = 1

        data[brdf_bands,:] = data[brdf_bands,:]*correction_factor

    elif dimension == 'column':
        # index= 300
        # data = hy_obj.get_column(index)

        brdf = fvol[np.newaxis,:]*hy_obj.ancillary['k_vol'][:,[index]]
        brdf+= fgeo[np.newaxis,:]*hy_obj.ancillary['k_geom'][:,[index]]
        brdf+= fiso[np.newaxis,:]

        brdf_nadir = fvol[np.newaxis,:]*hy_obj.ancillary['k_vol_nadir'][:,[index]]
        brdf_nadir+= fgeo[np.newaxis,:]*hy_obj.ancillary['k_geom_nadir'][:,[index]]
        brdf_nadir+= fiso[np.newaxis,:]

        correction_factor = brdf_nadir/brdf
        correction_factor[hy_obj.mask['brdf'][:,index],:] = 1

        data[:,brdf_bands] = data[:,brdf_bands]*correction_factor

    elif dimension == 'band':
        #index= 8
        #data = hy_obj.get_band(index)
        fvol, fgeo, fiso  = hy_obj.brdf['coeffs'][index]
        brdf = fvol*hy_obj.ancillary['k_vol']
        brdf += fgeo*hy_obj.ancillary['k_geom']
        brdf+=fiso

        brdf_nadir = fvol*hy_obj.ancillary['k_vol_nadir']
        brdf_nadir+= fgeo*hy_obj.ancillary['k_geom_nadir']
        brdf_nadir+= fiso

        correction_factor = brdf_nadir/brdf
        correction_factor[~hy_obj.mask['brdf']] = 1
        data= data* correction_factor

    elif dimension == 'chunk':
        #index = 200,501,3000,3501
        x1,x2,y1,y2 = index
        #data = hy_obj.get_chunk(x1,x2,y1,y2)

        brdf = fvol[np.newaxis,np.newaxis,:]*hy_obj.ancillary['k_vol'][y1:y2,x1:x2,np.newaxis]
        brdf+= fgeo[np.newaxis,np.newaxis,:]*hy_obj.ancillary['k_geom'][y1:y2,x1:x2,np.newaxis]
        brdf+= fiso[np.newaxis,np.newaxis,:]

        brdf_nadir = fvol[np.newaxis,np.newaxis,:]*hy_obj.ancillary['k_vol_nadir'][y1:y2,x1:x2,np.newaxis]
        brdf_nadir+= fgeo[np.newaxis,np.newaxis,:]*hy_obj.ancillary['k_geom_nadir'][y1:y2,x1:x2,np.newaxis]
        brdf_nadir+= fiso[np.newaxis,np.newaxis,:]

        correction_factor = brdf_nadir/brdf
        correction_factor[hy_obj.mask['brdf'][y1:y2,x1:x2]] = 1

        data[:,:,brdf_bands] = data[:,:,brdf_bands]*correction_factor

    elif dimension == 'pixels':
        #index = [[2000,2001],[200,501]]
        y,x = index
        #data = hy_obj.get_pixels(y,x)

        brdf = fvol[np.newaxis,:]*hy_obj.ancillary['k_vol'][y,x,np.newaxis]
        brdf+= fgeo[np.newaxis,:]*hy_obj.ancillary['k_geom'][y,x,np.newaxis]
        brdf+= fiso[np.newaxis,:]

        brdf_nadir = fvol[np.newaxis,:]*hy_obj.ancillary['k_vol_nadir'][y,x,np.newaxis]
        brdf_nadir+= fgeo[np.newaxis,:]*hy_obj.ancillary['k_geom_nadir'][y,x,np.newaxis]
        brdf_nadir+= fiso[np.newaxis,:]

        correction_factor = brdf_nadir/brdf
        correction_factor[hy_obj.mask['brdf'][y,x]] = 1

        data[:,brdf_bands] = data[:,brdf_bands]*correction_factor

    return data
