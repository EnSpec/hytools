''' Custom masks

'''
from scipy.ndimage.morphology import binary_erosion
import numpy as np

def topo1(hy_obj):
    ndvi = hy_obj.ndi()
    cosine_i = hy_obj.cosine_i()
    slope = hy_obj.get_anc('slope')
    mask = (ndvi > 0.2) & hy_obj.mask['no_data'] & (cosine_i > .12) & (slope > .03)
    return mask

def ndvi_threshold(hy_obj,args):
    ndvi = hy_obj.ndi()
    return (ndvi > args['ndvi_min']) & (ndvi < args['ndvi_max'])

def neon_edge(hy_obj,radius = 30):
    '''
    Mask artifacts in NEON images around edges.
    '''
    y_grid, x_grid = np.ogrid[-radius: radius + 1, -radius: radius + 1]
    window =  (x_grid**2 + y_grid**2 <= radius**2).astype(np.float)
    buffer_edge = binary_erosion(hy_obj.mask['no_data'], window)

    return ~buffer_edge

mask_dict = {'topo1' : topo1,
             'ndvi' : ndvi_threshold,
             'neon_edge' : neon_edge}