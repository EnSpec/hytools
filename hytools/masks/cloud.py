''' Cloud masks

'''
from scipy.ndimage import median_filter
import numpy as np


def zhai_cloud(hy_obj,cloud,shadow,T1=0.01,t2=.1,t3=.25,t4=.5,T7= 9,T8= 9):
    '''This function replicates the method of Zhai et al. (2018) for detecting clouds and shadows in
    multispectral and hyperspectral imagery but does not apply shadow spatial refinement.

    Suggested values for coefficients and params:
        T1 : 0.01, 0.1, 1, 10, 100
        t2 : 1/10, 1/9, 1/8, 1/7, 1/6, 1/5, 1/4, 1/3, 1/2
        t3 : 1/4, 1/3, 1/2, 2/3, 3/4
        t4 : 1/2, 2/3, 3/4, 4/5, 5/6
        T7 : 3, 5, 7, 9, 11
        T8 : 3, 5, 7, 9, 11

    Zhai, H., Zhang, H., Zhang, L., & Li, P. (2018).
    Cloud/shadow detection based on spectral indices for multi/hyperspectral optical remote sensing imagery.
    ISPRS journal of photogrammetry and remote sensing, 144, 235-253.
    https://doi.org/10.1016/j.isprsjprs.2018.07.006

    Args:
        hy_obj : HyTools data container object:
        cloud (bool): Detect clouds.
        shadow (bool): Detect clouds.
        T1 (float): Threshold T1.
        t2 (float): Adjusting coefficient t2.
        t3 (float): Adjusting coefficient t3.
        t4 (float): Adjusting coefficient t4.
        T7 (float): Parameter T7.
        T8 (float): Parameter T8.

    Returns:
        mask (nd.array): Boolean array where detected clouds and/or shadows = True.

    '''

    blue= hy_obj.get_wave(440)
    green= hy_obj.get_wave(550)
    red= hy_obj.get_wave(660)
    nir = hy_obj.get_wave(850)

    #If SWIR not available
    if hy_obj.wavelengths.max() < 1570:
        # Zhai et al. 2018 Eq. 1a,b
        CI_1 = (3*nir)/(blue+green+red)
        CI_2 = (blue+green+red+nir)/4
        # Zhai et al. 2018 Eq. 3
        CSI = nir

    else:
        swir1 = hy_obj.get_wave(1570)
        swir2= hy_obj.get_wave(2110)
        # Zhai et al. 2018 Eq. 1a,b
        CI_1 = (nir+ 2*swir1)/(blue+green+red)
        CI_2 = (blue+green+red+nir+swir1+swir2)/6
        # Zhai et al. 2018 Eq. 3
        CSI = (nir + swir1)/2

    # Zhai et al. 2018 Eq.5
    T2 = np.mean(CI_2[hy_obj.mask['no_data']]) + t2*(np.max(CI_2[hy_obj.mask['no_data']])-np.mean(CI_2[hy_obj.mask['no_data']]))
    # Zhai et al. 2018 Eq.6
    T3 = np.min(CSI[hy_obj.mask['no_data']]) + t3*(np.mean(CSI[hy_obj.mask['no_data']])-np.min(CSI[hy_obj.mask['no_data']]))
    # Zhai et al. 2018 Eq.7
    T4 = np.min(blue[hy_obj.mask['no_data']]) + t4*(np.mean(blue[hy_obj.mask['no_data']])-np.min(blue[hy_obj.mask['no_data']]))

    mask = np.zeros((hy_obj.lines,hy_obj.columns)).astype(bool)

    if cloud:
        clouds = (np.abs(CI_1) < T1) | (CI_2 >  T2)
        clouds = median_filter(clouds, T7)
        mask[clouds] = True

    if shadow:
        shadows = (CSI<T3) & (blue<T4)
        shadows = median_filter(shadows,T8)
        mask[shadows] = True

    return mask


















