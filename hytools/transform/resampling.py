# -*- coding: utf-8 -*-
""" Spectral resampling functions.

"""
import numpy as np

def gaussian(x,mu,fwhm):
    """
    Args:
        x (numpy.ndarray): Values along which to generate gaussian..
        mu (float): Mean of the gaussian function..
        fwhm (float): Full width half maximum..

    Returns:
        numpy.ndarray: Gaussian along input range.

    """

    c = fwhm/(2* np.sqrt(2*np.log(2)))
    return np.exp(-1*((x-mu)**2/(2*c**2)))

def resample_coeff(in_wave,in_fwhm,out_wave,out_fwhm, spacing = 1):
    """Given a set of source and destination wavelengths and FWHMs this
    function caculates the relative contribution or each input wavelength
    to the output wavelength. It assumes that both input and output
    response functions follow a gaussian distribution.

    All inputs shoud be provide in nanometers.


    Args:
        in_wave (list): Input wavelength centers.
        in_fwhm (list): Input full width half maxes.
        out_wave (list): Output wavelength centers.
        out_fwhm (list): Output full width half maxes.
        spacing (int, optional): Resolution at which to model the
                    spectral response functions. Defaults to 1.

    Returns:
        numpy.ndarray: Transform coeffiecients.

    """

    out_matrix = []
    min_spectrum = min(out_wave.min(),in_wave.min())//100 *100 - 100
    max_spectrum = 100 + max(out_wave.max(),in_wave.max())//100 *100
    oneNM = np.arange(min_spectrum,max_spectrum,spacing)

    for wave,fwhm, in zip(out_wave,out_fwhm):
        a =  gaussian(oneNM,wave,fwhm)
        out_matrix.append(np.divide(a,np.sum(a)))
    out_matrix = np.array(out_matrix)

    # For each source wavelength generate the gaussion response
    in_matrix = []
    for wave,fwhm in zip(in_wave,in_fwhm):
        in_matrix.append(gaussian(oneNM ,wave,fwhm))
    in_matrix = np.array(in_matrix)

    # Calculate the relative contribution of each source response function
    ratio =  in_matrix/in_matrix.sum(axis=0)
    ratio[np.isnan(ratio)] = 0
    ratio2 = np.einsum('ab,cb->acb',ratio,out_matrix)

    # Calculate the relative contribution of each input wavelength
    # to each destination wavelength
    coeffs = np.trapz(ratio2)

    return coeffs
