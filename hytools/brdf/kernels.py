# -*- coding: utf-8 -*-
"""
This module contains functions to calculate BRDF scattering kernels.

Equations and constants can be found in the following papers:

Colgan, M. S., Baldeck, C. A., Feret, J. B., & Asner, G. P. (2012).
Mapping savanna tree species at ecosystem scales using support vector machine classification
and BRDF correction on airborne hyperspectral and LiDAR data.
Remote Sensing, 4(11), 3462-3480.
https://doi.org/10.3390/rs4113462

Lucht, W., Schaaf, C. B., & Strahler, A. H. (2000).
An algorithm for the retrieval of albedo from space using semiempirical BRDF models.
IEEE Transactions on Geoscience and Remote sensing, 38(2), 977-998.
https://doi.org/10.1109/36.841980

Maignan, F., Bréon, F. M., & Lacaze, R. (2004).
Bidirectional reflectance of Earth targets: Evaluation of analytical
models using a large set of spaceborne measurements with emphasis on the Hot Spot.
Remote Sensing of Environment, 90(2), 210-220.
https://doi.org/10.1016/j.rse.2003.12.006

Roujean, J. L., Leroy, M., & Deschamps, P. Y. (1992).
A bidirectional reflectance model of the Earth's surface for the correction
of remote sensing data.
Journal of Geophysical Research: Atmospheres, 97(D18), 20455-20468.
https://doi.org/10.1029/92JD01411

Schlapfer, D., Richter, R., & Feingersh, T. (2015).
Operational BRDF effects correction for wide-field-of-view optical scanners (BREFCOR).
IEEE Transactions on Geoscience and Remote Sensing, 53(4), 1855-1864.
https://doi.org/10.1109/TGRS.2014.2349946

Wanner, W., Li, X., & Strahler, A. H. (1995).
On the derivation of kernels for kernel-driven models of bidirectional reflectance.
Journal of Geophysical Research: Atmospheres, 100(D10), 21077-21089.
https://doi.org/10.1029/95JD02371

Zhang, X., Jiao, Z., Dong, Y., Zhang, H., Li, Y., He, D., ... & Chang, Y. (2018).
Potential investigation of linking PROSAIL with the ross-li BRDF model for
vegetation characterization.
Remote Sensing, 10(3), 437.
https://doi.org/10.3390/rs10030437SSS

"""
import numpy as np

def calc_geom_kernel(solar_az,solar_zn,sensor_az,sensor_zn,kernel,b_r=1.,h_b =2.):
    """Calculate geometric scattering kernel.
       Constants b_r (b/r) and h_b (h/b) from Colgan et al. RS 2012
       Alternatives include MODIS specification:
           b/r : sparse: 1, dense: 2.5
           h/b : sparse, dense : 2

       All input geometry units must be in radians.

    Args:
        solar_az (numpy.ndarray): Solar azimuth angle.
        solar_zn (numpy.ndarray): Solar zenith angle.
        sensor_az (numpy.ndarray): Sensor view azimuth angle.
        sensor_zn (numpy.ndarray): Sensor view zenith angle.
        kernel (str): Li geometric scattering kernel type [li_dense,li_sparse, roujean].
        b_r (float, optional): Object height. Defaults to 10.
        h_b (float, optional): Object shape. Defaults to 2.

    Returns:
        numpy.ndarray: Geometric scattering kernel.

    """

    relative_az = sensor_az - solar_az

    # Eq. 37,52. Wanner et al. JGRA 1995
    solar_zn_p = np.arctan(b_r * np.tan(solar_zn))
    sensor_zn_p = np.arctan(b_r * np.tan(sensor_zn))
    # Eq 50. Wanner et al. JGRA 1995
    D = np.sqrt((np.tan(solar_zn_p)**2) + (np.tan(sensor_zn_p)**2) - 2*np.tan(solar_zn_p)*np.tan(sensor_zn_p)*np.cos(relative_az))
    # Eq 49. Wanner et al. JGRA 1995
    t_num = h_b * np.sqrt(D**2 + (np.tan(solar_zn_p)*np.tan(sensor_zn_p)*np.sin(relative_az))**2)
    t_denom = (1/np.cos(solar_zn_p))  + (1/np.cos(sensor_zn_p))
    t = np.arccos(np.clip(t_num/t_denom,-1,1))
    # Eq 33,48. Wanner et al. JGRA 1995
    O = (1/np.pi) * (t - np.sin(t)*np.cos(t)) * t_denom
    # Eq 51. Wanner et al. JGRA 1995
    cos_phase_p =  np.cos(solar_zn_p)*np.cos(sensor_zn_p) + np.sin(solar_zn_p)*np.sin(sensor_zn_p)*np.cos(relative_az)

    if kernel == 'li_sparse':
        # Eq 32. Wanner et al. JGRA 1995
        k_geom = O - (1/np.cos(solar_zn_p)) - (1/np.cos(sensor_zn_p)) + .5*(1+ cos_phase_p) * (1/np.cos(sensor_zn_p))
    elif kernel == 'li_dense':
        # Eq 47. Wanner et al. JGRA 1995
        k_geom = (((1+cos_phase_p) * (1/np.cos(sensor_zn_p)))/ (t_denom - O)) - 2
    elif kernel == 'li_sparse_r':
        # Eq 39. Lucht et al. TGRS 2000
        k_geom = O - (1/np.cos(solar_zn_p)) - (1/np.cos(sensor_zn_p)) + .5*(1+ cos_phase_p) * (1/np.cos(sensor_zn_p)) * (1/np.cos(solar_zn_p))
    elif kernel == 'li_dense_r':
        # Eq 5. Zhang et al. RS 2018 <-- Find a more original reference
        k_geom = (((1+cos_phase_p) * (1/np.cos(sensor_zn_p)) * (1/np.cos(solar_zn_p)))/ (t_denom - O)) - 2
    elif kernel == 'roujean':
        # Eq 2 Roujean et al. JGR 1992
        k_geom1 = (1/(2*np.pi)) * ((np.pi - relative_az)*np.cos(relative_az)+np.sin(relative_az)) *np.tan(solar_zn)*np.tan(sensor_zn)
        k_geom2 =  (1/np.pi) * (np.tan(solar_zn) + np.tan(sensor_zn) + np.sqrt(np.tan(solar_zn)**2 + np.tan(sensor_zn)**2 - 2*np.tan(solar_zn)*np.tan(sensor_zn)*np.cos(relative_az)))
        k_geom = k_geom1 -  k_geom2
    else:
        print("Unrecognized kernel type: %s" %  kernel)
        k_geom = None
    return k_geom


def calc_volume_kernel(solar_az,solar_zn,sensor_az,sensor_zn,kernel):
    """Calculate volume scattering kernel.
       All input geometry units must be in radians.

    Args:
        solar_az (numpy.ndarray): Solar azimuth angle.
        solar_zn (numpy.ndarray): Solar zenith angle.
        sensor_az (numpy.ndarray): Sensor view azimuth angle.
        sensor_zn (numpy.ndarray): Sensor view zenith angle.
        kernel (str): Volume scattering kernel type [ross_thick,ross_thin].

    Returns:
        numpy.ndarray: Volume scattering kernel.

    """

    relative_az = sensor_az - solar_az

    # Eq 2. Schlapfer et al. IEEE-TGARS 2015
    phase = np.arccos(np.cos(solar_zn)*np.cos(sensor_zn) + np.sin(solar_zn)*np.sin(sensor_zn)*  np.cos(relative_az))

    if kernel == 'ross_thin':
        # Eq 13. Wanner et al. JGRA 1995
        k_vol = ((np.pi/2 - phase)*np.cos(phase) + np.sin(phase))/ (np.cos(sensor_zn)*np.cos(solar_zn)) - (np.pi/2)
    elif kernel == 'ross_thick':
        # Eq 7. Wanner et al. JGRA 1995
        k_vol = ((np.pi/2 - phase)*np.cos(phase) + np.sin(phase))/ (np.cos(sensor_zn)*np.cos(solar_zn)) - (np.pi/4)
    elif kernel in ('hotspot','roujean'):
        # Eq 8 Roujean et al. JGR 1992
        k_vol1 = (4/(3*np.pi)) * (1/(np.cos(solar_zn) + np.cos(sensor_zn)))
        k_vol2 = (((np.pi/2) - phase) * np.cos(phase) + np.sin(phase))
        k_vol = k_vol1*(k_vol2- (1/3))
        if kernel == 'hotspot':
            # Eq. 12 Maignan et al. RSE 2004
            k_vol =  k_vol1* k_vol2 * (1 + (1 + (phase/np.radians(1.5)))**-1) - (1/3)
    else:
        print("Unrecognized kernel type: %s" % kernel)
        k_vol = None
    return k_vol
