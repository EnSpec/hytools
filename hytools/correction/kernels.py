"""
This module contains functions to calculate BRDF scattering kernels.
Equations and constants can be found in the following papers:

Colgan, M. S., Baldeck, C. A., Feret, J. B., & Asner, G. P. (2012).
Mapping savanna tree species at ecosystem scales using support vector machine classification
and BRDF correction on airborne hyperspectral and LiDAR data.
Remote Sensing, 4(11), 3462-3480.

Schlapfer, D., Richter, R., & Feingersh, T. (2015).
Operational BRDF effects correction for wide-field-of-view optical scanners (BREFCOR).
IEEE Transactions on Geoscience and Remote Sensing, 53(4), 1855-1864.

Wanner, W., Li, X., & Strahler, A. H. (1995).
On the derivation of kernels for kernel-driven models of bidirectional reflectance.
Journal of Geophysical Research: Atmospheres, 100(D10), 21077-21089.
"""
import numpy as np

def generate_geom_kernel(solar_az,solar_zn,sensor_az,sensor_zn,kernel,b_r=10.,h_b =2.):
    '''Calculate the geometric scattering kernel.
       Constants b_r (b/r) and h_b (h/b) from Colgan et al. RS 2012
       All input geometry in radians.

    :param solar_az: Solar azimuth angle
    :type solar_az: np.ndarray
    :param solar_zn: Solar zenith angle
    :type solar_zn: np.ndarray, float
    :param sensor_az: Sensor view azimuth angle
    :type sensor_az: np.ndarray, float
    :param sensor_zn: Sensor view zenith angle
    :type sensor_zn: np.ndarray, float
    :param kernel: Li geometric scattering kernel type [li_dense,li_parse]
    :type kernel: str
    :param b_r: Object height, defaults to 10.
    :type b_r: float, optional
    :param h_b: Object shape, defaults to 2
    :type h_b: float, optional
    :return: Geometric scattering kernel
    :rtype: nd.array, float
    '''
    relative_az = sensor_az - solar_az

    # Eq. 37,52. Wanner et al. JGRA 1995
    solar_zn_ = np.arctan(b_r * np.tan(solar_zn))
    sensor_zn_ = np.arctan(b_r * np.tan(sensor_zn))
    # Eq 50. Wanner et al. JGRA 1995
    D = np.sqrt((np.tan(solar_zn_)**2) + (np.tan(sensor_zn_)**2) - 2*np.tan(solar_zn_)*np.tan(sensor_zn_)*np.cos(relative_az))
    # Eq 49. Wanner et al. JGRA 1995
    t_num = h_b * np.sqrt(D**2 + (np.tan(solar_zn_)*np.tan(sensor_zn_)*np.sin(relative_az))**2)
    t_denom = (1/np.cos(solar_zn_))  + (1/np.cos(sensor_zn_))
    t = np.arccos(np.clip(t_num/t_denom,-1,1))
    # Eq 33,48. Wanner et al. JGRA 1995
    O = (1/np.pi) * (t - np.sin(t)*np.cos(t)) * t_denom
    # Eq 51. Wanner et al. JGRA 1995
    cos_phase_ =  np.cos(solar_zn_)*np.cos(sensor_zn_) + np.sin(solar_zn_)*np.sin(sensor_zn_)*np.cos(relative_az)

    if kernel == 'li_sparse':
        # Eq 32. Wanner et al. JGRA 1995
        k_geom = O - (1/np.cos(solar_zn_)) - (1/np.cos(sensor_zn_)) + .5*(1+ cos_phase_) * (1/np.cos(sensor_zn_))
    elif kernel == 'li_dense':
        # Eq 47. Wanner et al. JGRA 1995
        k_geom = (((1+cos_phase_) * (1/np.cos(sensor_zn_)))/ (t_denom - O)) - 2
    return k_geom


def generate_volume_kernel(solar_az,solar_zn,sensor_az,sensor_zn,kernel):
    '''Calculate the volumetric scattering kernel.
       Input geometry in radians.

    :param solar_az: Solar azimuth angle
    :type solar_az: np.ndarray
    :param solar_zn: Solar zenith angle
    :type solar_zn: np.ndarray, float
    :param sensor_az: Sensor view azimuth angle
    :type sensor_az: np.ndarray, float
    :param sensor_zn: Sensor view zenith angle
    :type sensor_zn: np.ndarray, float
    :param kernel: Volume scattering kernel type [ross_thick,ross_thin]
    :type kernel: str
    :return: Volumetric scattering kernel
    :rtype: nd.array, float
    '''
    relative_az = sensor_az - solar_az

    # Eq 2. Schlapfer et al. IEEE-TGARS 2015
    phase = np.arccos(np.cos(solar_zn)*np.cos(sensor_zn) + np.sin(solar_zn)*np.sin(sensor_zn)*  np.cos(relative_az))

    if kernel == 'ross_thin':
        # Eq 13. Wanner et al. JGRA 1995
        k_vol = ((np.pi/2 - phase)*np.cos(phase) + np.sin(phase))/ (np.cos(sensor_zn)*np.cos(solar_zn)) - np.pi/2
    elif kernel == 'ross_thick':
        # Eq 7. Wanner et al. JGRA 1995
        k_vol = ((np.pi/2 - phase)*np.cos(phase) + np.sin(phase))/ (np.cos(sensor_zn)*np.cos(solar_zn)) - np.pi/4
    return k_vol
