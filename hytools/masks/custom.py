''' Custom masks

'''
from scipy.ndimage.morphology import binary_erosion
import numpy as np

def ndi(hy_obj,args):
    mask = hy_obj.ndi(args['band_1'],args['band_2'])
    mask = (mask >= float(args['min'])) & (mask <= float(args['max']))
    return mask

def ancillary(hy_obj,args):
    ''' Mask ancillary datasets based off min and max threshold

    '''
    if args['name'] == 'cosine_i':
        mask= hy_obj.cosine_i()
    else:
        mask = hy_obj.get_anc(args['name'])

    mask = (mask >= float(args['min'])) & (mask <= float(args['max']))
    return mask

def neon_edge(hy_obj,args):
    '''
    Mask artifacts in NEON images around edges.
    '''
    radius =args['radius']
    y_grid, x_grid = np.ogrid[-radius: radius + 1, -radius: radius + 1]
    window =  (x_grid**2 + y_grid**2 <= radius**2).astype(np.float)
    buffer_edge = binary_erosion(hy_obj.mask['no_data'], window).astype(bool)
    return ~buffer_edge

def kernel_finite(hy_obj,args):
    '''
    Create NDVI bin class mask
    '''
    k_vol  = hy_obj.volume_kernel(hy_obj.brdf['volume'])
    k_geom = hy_obj.geom_kernel(hy_obj.brdf['geometric'],
                                b_r=hy_obj.brdf["b/r"],
                                h_b =hy_obj.brdf["h/b"])
    mask = np.isfinite(k_vol) & np.isfinite(k_geom)
    return mask

mask_dict = {'ndi' : ndi,
             'neon_edge' : neon_edge,
             'kernel_finite': kernel_finite,
             'ancillary':  ancillary }

def mask_create(hy_obj,masks):
    ''' Combine a series of boolean masks using an
    and operator
    '''
    mask = np.copy(hy_obj.mask['no_data'])

    for mask_name,args in masks:
        mask &= mask_dict[mask_name](hy_obj,args)

    return mask