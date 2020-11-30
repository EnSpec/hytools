"""
This module contains functions to calculate BRDF scattering kernels.
Equations and constants can be found in the following papers:

Colgan, M. S., Baldeck, C. A., Feret, J. B., & Asner, G. P. (2012).
Mapnp.ping savanna tree species at ecosystem scales using support vector machine classification
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

def generate_geom_kernel(solar_az,solar_zn,sensor_az,sensor_zn,kernel):
    '''Calculate the Li geometric scattering kernel.
       Constants from Colgan et al. RS 2012

    :param solar_az: Solar azimuth angle in radians
    :type solar_az: np.ndarray
    :param solar_zn: Solar zenith angle in radians
    :type solar_zn: np.ndarray, float
    :param sensor_az: Sensor view azimuth angle in radians
    :type sensor_az: np.ndarray, float
    :param sensor_zn: Sensor view zenith angle in radians
    :type sensor_zn: np.ndarray, float
    :param kernel: Li geometric scattering kernel type [dense,sparse]
    :type kernel: str
    :return: Geometric scattering kernel
    :rtype: nd.array, float
    '''

    relative_az = sensor_az - solar_az

    # Eq. 37,52. Wanner et al. JGRA 1995
    solar_zn_ = np.arctan(10* np.tan(solar_zn))
    sensor_zn_ = np.arctan(10* np.tan(sensor_zn))
    # Eq 50. Wanner et al. JGRA 1995
    D = np.sqrt((np.tan(solar_zn_)**2) + (np.tan(sensor_zn_)**2) - 2*np.tan(solar_zn_)*np.tan(sensor_zn_)*np.cos(relative_az))
    # Eq 49. Wanner et al. JGRA 1995
    t_num = 2. * np.sqrt(D**2 + (np.tan(solar_zn_)*np.tan(sensor_zn_)*np.sin(relative_az))**2)
    t_denom = (1/np.cos(solar_zn_))  + (1/np.cos(sensor_zn_))
    t = np.arccos(np.clip(t_num/t_denom,-1,1))
    # Eq 33,48. Wanner et al. JGRA 1995
    O = (1/np.pi) * (t - np.sin(t)*np.cos(t)) * t_denom
    # Eq 51. Wanner et al. JGRA 1995
    cos_phase_ =  np.cos(solar_zn_)*np.cos(sensor_zn_) + np.sin(solar_zn_)*np.sin(sensor_zn_)*np.cos(relative_az)

    if kernel == 'sparse':
        # Eq 32. Wanner et al. JGRA 1995
        k_geom = O - (1/np.cos(solar_zn_)) - (1/np.cos(sensor_zn_)) + .5*(1+ cos_phase_) * (1/np.cos(sensor_zn_))
    elif kernel == 'dense':
        # Eq 47. Wanner et al. JGRA 1995
        k_geom = (((1+cos_phase_) * (1/np.cos(sensor_zn_)))/ (t_denom - O)) - 2
    return k_geom


def generate_volume_kernel(solar_az,solar_zn,sensor_az,sensor_zn,kernel,hotspot = False):
    '''Calculate the Ross volumetric scattering kernel.

    :param solar_az: Solar azimuth angle in radians
    :type solar_az: np.ndarray
    :param solar_zn: Solar zenith angle in radians
    :type solar_zn: np.ndarray, float
    :param sensor_az: Sensor view azimuth angle in radians
    :type sensor_az: np.ndarray, float
    :param sensor_zn: Sensor view zenith angle in radians
    :type sensor_zn: np.ndarray, float
    :param kernel: Volume scattering kernel type [thick,thin]
    :type kernel: str
    :param hotspot: Apply hotspot correction, defaults to False
    :type hotspot: bool, optional
    :return: Volumetric scattering kernel
    :rtype: nd.array, float
    '''
    relative_az = sensor_az - solar_az

    # Eq 2. Schlapfer et al. IEEE-TGARS 2015
    phase = np.arccos(np.cos(solar_zn)*np.cos(sensor_zn) + np.sin(solar_zn)*np.sin(sensor_zn)*  np.cos(relative_az))

    if kernel == 'thick':
        # Eq 13. Wanner et al. JGRA 1995
        k_vol = ((np.pi/2 - phase)*np.cos(phase) + np.sin(phase))/ (np.cos(sensor_zn)*np.cos(solar_zn)) - np.pi/4
    elif kernel == 'thin':
        # Eq 13. Wanner et al. JGRA 1995
        k_vol = ((np.pi/2 - phase)*np.cos(phase) + np.sin(phase))/ (np.cos(sensor_zn)*np.cos(solar_zn)) - np.pi/2

    if hotspot:
        # TODO: Add citation
        k_vol = (k_vol + (1/3))* (1 + 1/(1 + (phase/np.radians(1.3))))
    return k_vol
