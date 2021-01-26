import hytools as ht
import glob
import ray
import warnings
import numpy as np
import matplotlib.pyplot as plt
import json
from hytools.io import *
from hytools.correction.brdf import *
import os
import copy
from scipy.ndimage.morphology import binary_erosion
import pandas as pd

files = glob.glob('/home/chlus/dev_hytools/data/neon/*h5')

if ray.is_initialized():
    ray.shutdown()
ray.init(num_cpus = len(files))

# Create actors and load files
HyTools = ray.remote(ht.HyTools)
actors = [HyTools.remote() for file in files]
_ = ray.get([a.read_file.remote(file,'neon') for a,file in zip(actors,files)])

#test object
hy_obj = ht.HyTools()
hy_obj.read_file(files[0],'neon')

args = {"samp_perc": 0.01,
            "buffer": 30,
            "ndvi_ir": 850,
            "ndvi_red": 660,
            "cos_i_min": 0.12,
            "slope_min":0.87}

def neon_masks(hy_obj):
    '''Generate 3 masks for the NEON data:
        1. NDVI
        2. X % sampling mask, with edges excluded
        3. Topographic correction mask, 1 = correct
    '''

    ir = hy_obj.get_wave(args['ndvi_ir'])
    red = hy_obj.get_wave(args['ndvi_red'])
    ndvi = (ir-red)/(ir+red)
    hy_obj.mask['ndvi'] = ndvi > 0.05

    buffer_mask = binary_erosion(ir != hy_obj.no_data,
                                np.ones((args['buffer']*2,args['buffer']*2)),
                                border_value=0)

    # Select 'sample_perc' % of pixels for modeling
    # This can probably be written more concisely
    sub_samples = np.zeros(buffer_mask.shape).astype(bool)
    idx = np.array(np.where(buffer_mask== True)).T
    idxRand= idx[np.random.choice(range(len(idx)),int(len(idx)*args['samp_perc']), replace = False)].T
    sub_samples[idxRand[0],idxRand[1]] = True
    hy_obj.mask['sample']  = sub_samples

    cosine_i = hy_obj.cosine_i()
    slope = hy_obj.get_anc('slope',radians=False)
    hy_obj.mask['topo'] =  (slope>= args['slope_min']) & (cosine_i >= args['cos_i_min'])


#Generate masks
_ = ray.get([a.do.remote(neon_masks) for a in actors])

waves=  [480, 560, 660, 850, 950, 1050, 1240, 1650, 2217]

def topo_correction(hy_obj):
    topo_coeffs = {}
    topo_coeffs['C'] = []
    topo_coeffs['slope'] = []

    topo_mask = hy_obj.mask['topo'] & hy_obj.mask['sample']


    ir = hy_obj.get_wave(args['ndvi_ir'])
    red = hy_obj.get_wave(args['ndvi_red'])
    ndvi = (ir-red)/(ir+red)

    cos_i = hy_obj.cosine_i()
    cos_i = cos_i[topo_mask]
    X = np.vstack([cos_i, np.ones(cos_i.shape)]).T

    c1 = np.cos(hy_obj.get_anc('solar_zn')) * np.cos(hy_obj.get_anc('solar_zn'))
    c1 = c1[topo_mask]

    k_vol = hy_obj.volume_kernel('ross_thick')[hy_obj.mask['sample']]
    k_geom = hy_obj.geom_kernel('li_sparse')[hy_obj.mask['sample']]

    outlier_data= pd.DataFrame()
    outlier_data['vol'] = k_vol
    outlier_data['geom'] = k_geom
    outlier_data['ndvi'] = ndvi[hy_obj.mask['sample']]
    outlier_data['image']  = os.path.basename(hy_obj.file_name)

    for wave in waves:
        band = hy_obj.get_wave(wave).astype(float)
        y = band[topo_mask]

        # Eq 7. Soenen et al., IEEE TGARS 2005
        slope, intercept = np.linalg.lstsq(X, y)[0].flatten()
        # Eq 8. Soenen et al., IEEE TGARS 2005
        C= intercept/slope

        # Set a large number if slope is zero
        if not np.isfinite(C):
          C = 100000.0

        topo_coeffs['C'].append(C)
        topo_coeffs['slope'].append(slope)

        #Apply correction to subset of samples
        correctionFactor = (c1 * C)/(cos_i * C)
        band[topo_mask] = band[topo_mask] * correctionFactor
        outlier_data[wave] =   band[hy_obj.mask['sample']]

    hy_obj.topo_coeffs = topo_coeffs
    return outlier_data


    outlier_data = ray.get([a.do.remote(topo_correction) for a in actors])
    outlier_data = pd.concat(outlier_data)
    outlier_data.replace([np.inf, -np.inf], np.nan,inplace =True)
    outlier_data.dropna(axis = 0,how ='any',inplace= True)

    quantiles = outlier_data.quantile([.05,.33,.67,.95])

    outlier_data['bins_geom'] = np.digitize(outlier_data['geom'],quantiles['geom'])
    outlier_data['bins_vol'] = np.digitize(outlier_data['vol'],quantiles['vol'])
    outlier_data['bins_ndvi'] = np.digitize(outlier_data['ndvi'],[.15,.4,.7,.95])

    median = outlier_data.groupby(by = ['image','bins_geom','bins_vol','bins_ndvi']).median()



    pca = PCA(n_components=2)
    pca.fit(outlier_data.drop('image',axis=1))
    trans  =pd.DataFrame(pca.transform(outlier_data.drop('image',axis=1)))


