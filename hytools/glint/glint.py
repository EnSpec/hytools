import numpy as np
import ray
from scipy import stats
from scipy.linalg import inv, svd
from scipy.interpolate import interp1d
from scipy.optimize import leastsq 
from scipy.optimize import minimize
from ..masks import mask_create
from ..misc import progbar, pairwise, set_glint
from .lut import set_LUT 
from .ref_n import REFRACTIVE_INDICES


def set_glint_parameters(actors,config_dict):
    # Assign glint dict
    glint_dict = config_dict['glint']

    # Set Glint dict
    _ = ray.get(
        [a.do.remote(set_glint,glint_dict) for a in actors]
    )

    # If method is the sky-glint, optimize windspeed
    if glint_dict['type'] == 'SkyGlint':
        # Load Lookup table
        _ = ray.get(
            [a.do.remote(set_LUT) 
            for a in actors]
        )

    # Add glint correction
    _ = ray.get([a.do.remote(lambda x: x.corrections.append('glint')) for a in actors])


def apply_glint_correct(hy_obj,data,dimension,index):
    ''' Corrects glint based on the specified algorithm in the config.
        Options include: 
            Hochberg et al., 2003: Hochberg
            Gao et al., 2021: Gao
            Hedley + Sky glint: SkyGlint
            ...
    '''

    # Create water mask to guide calculation
    hy_obj.mask['water'] = mask_create(hy_obj,hy_obj.glint['calc_mask'])

    # Perform one of the corrections
    if hy_obj.glint['type'] == 'Hochberg':
        data = hochberg_2003_correction(hy_obj,data,dimension,index)

    elif hy_obj.glint['type'] == 'SkyGlint':
        # Optimize windspeed
        if isinstance(hy_obj.glint.get('correction'), type(None)):
            _ = optimize_windspeed_correction(hy_obj)
        data = sky_sun_glint_correction(hy_obj,data,dimension,index)

    elif hy_obj.glint['type'] == 'Gao':
        data = gao_2021_correction(hy_obj,data,dimension,index)
    # Can add more corrections here

    return data


def hochberg_2003_correction(hy_obj,data,dimension,index):
    """
    Glint correction algorithm following:

    Hochberg, EJ, Andréfouët, S and Tyler, MR. 2003. 
    Sea surface correction of high spatial resolution Ikonos images to 
    improve bottom mapping in near‐shore environments.. 
    IEEE Transactions on Geoscience and Remote Sensing, 41: 1724–1729.
    """

    # Get the Reference values (1467, 1600, 2190)
    corr_wave = hy_obj.glint['correction_wave']

    corr = np.copy(hy_obj.get_wave(corr_wave))
    corr[~hy_obj.mask['water']] = 0 

    corr_min = np.percentile(corr[corr > 0], 1)

    if dimension == 'line':
        # Get SWIR difference
        swir_diff = (
            corr - corr_min
        )[index, :]

        # Mask SWIR difference
        mask = hy_obj.mask['water'][index,:]
        swir_diff[~mask] = 0

        bandnums = len(hy_obj.wavelengths) 
        swir_diff = np.repeat(
            swir_diff[:, np.newaxis], 
            bandnums, 
            axis=1
        )

        # Calculate correction
        data = data - swir_diff 

    elif dimension == 'column':
        # Get SWIR difference
        swir_diff = (
            corr - corr_min
        )[:, index]

        # Get Correction
        mask = hy_obj.mask['water'][:, index]
        swir_diff[~mask] = 0

        bandnums = len(hy_obj.wavelengths) 
        swir_diff = np.repeat(
            swir_diff[:, np.newaxis], 
            bandnums, 
            axis=1
        )

        # Calculate correction
        data = data - swir_diff 

    elif (dimension == 'band'):
        # Get SWIR difference
        swir_diff = (
            corr - corr_min
        )

        # Get Correction
        mask = hy_obj.mask['water']
        swir_diff[~mask] = 0

        data = data - swir_diff 

    elif dimension == 'chunk':
        # Get Index
        x1, x2, y1, y2 = index

        # Get SWIR difference
        swir_diff = (
            corr - corr_min
        )[y1:y2,x1:x2]

        # Get Correction
        mask = hy_obj.mask['water'][y1:y2, x1:x2]
        swir_diff[~mask] = 0

        bandnums = len(hy_obj.wavelengths) 
        swir_diff = np.repeat(
            swir_diff[:, :, np.newaxis], 
            bandnums, 
            axis=2
        )

        # Correct
        data = data - swir_diff

    elif dimension == 'pixels':
        y, x = index

        # Get SWIR difference
        swir_diff = (
            corr - corr_min
        )[y,x]

        # Get correction
        mask = hy_obj.mask['water'][y, x]
        swir_diff[~mask] = 0

        bandnums = len(hy_obj.wavelengths) 
        swir_diff = np.repeat(
            swir_diff[:, np.newaxis], 
            bandnums, 
            axis=1
        )

        # correct
        data = data - swir_diff 

    return data


def optimize_windspeed_correction(hy_obj):
    """
    Glint correction algorithm following:

    Extended Hochberg correction to include sky glint correction
    This is the optimization step to find the effective wind speed
    """
    corr_wave = hy_obj.get_wave(hy_obj.glint['correction_wave'])
    
    # Get solar zenith array
    solar_zn = np.degrees(hy_obj.get_anc('solar_zn'))

    # Get LUT
    lut = hy_obj.glint['lut']

    # Initial gueses for wind optimization
    x0 = [10]
    # bounds = [(0, 20)]
    bounds = hy_obj.glint['bounds']

    # Get water mask
    masks = hy_obj.mask['water']

    # Initialize correction
    correction = np.zeros([*corr_wave.shape, len(hy_obj.wavelengths)])
    windspeeds = np.zeros([*corr_wave.shape, len(hy_obj.wavelengths)])
    # Initialize iterator
    it = np.nditer(corr_wave, flags=['multi_index'])
    for pixel in it:
        y, x = it.multi_index

        if not masks[y, x]:
            continue

        # Get solar zenith
        sol_zen = solar_zn[y, x]

        # Optimize wind speed
        res = minimize(
             err,
             x0,
             args=(float(pixel), lut, sol_zen),
             method='tnc',
             tol=1e-9,
             bounds=bounds
        )
        glint = glint_spectrum(lut, float(res.x), sol_zen)
        correction[y, x, :] = glint
        windspeeds[y, x, :] = res.x 

    hy_obj.glint['correction'] = correction 
    hy_obj.glint['windspeed'] = windspeeds 


def sky_sun_glint_correction(hy_obj,data,dimension,index):
    """
    Glint correction algorithm following:

    This applies the correction to the sample of data
    """

    # Get correction wave
    corrections = hy_obj.glint['correction']

    if dimension == 'line':

        # Get slice of correction 
        correction = corrections[index, :]

        # Correct
        data = data - correction
        data[(data < 0) & (data != hy_obj.no_data)] = 0

    elif dimension == 'column':

        # Get slice of correction 
        correction = corrections[:, index]

        # Correction
        data = data - correction
        data[(data < 0) & (data != hy_obj.no_data)] = 0

    elif dimension == 'band':

        # Get slice of correction 
        correction = corrections[:, :, index]

        # Correct
        data = data - correction
        data[(data < 0) & (data != hy_obj.no_data)] = 0

    elif dimension == 'chunk':

        # Get slice of correction 
        correction = corrections[y1:y2, x1:x2, :]

        # Correct
        data = data - correction
        data[(data < 0) & (data != hy_obj.no_data)] = 0

    elif dimension == 'pixels':
        y, x = index

        # Get slice of correction 
        correction = corrections[y, x, :]

        # Correct
        data = data - correction
        data[(data < 0) & (data != hy_obj.no_data)] = 0

    return data


def gao_2021_correction(hy_obj,data,dimension,index):
    """
    Glint correction algorithm following:

    Gao BC, Li RR. 
    Correction of Sunglint Effects in High Spatial Resolution 
    Hyperspectral Imagery Using SWIR or NIR Bands and Taking Account of 
    Spectral Variation of Refractive Index of Water. 
    Adv Environ Eng Res 2021;2(3):16; doi:10.21926/aeer.2103017.
    """

    # Get the Reference values (1467, 1600, 2190) and band number
    corr_wave = hy_obj.glint['correction_wave']
    corr_num = np.argmin(np.abs(hy_obj.wavelengths - corr_wave))

    # Get simulated 0 theta fresnel reflectence
    # REF_IND is a hardcoded array of n by wavelength
    B_Simu = fresnel_spectra(10**-5, hy_obj.wavelengths, REFRACTIVE_INDICES)
    B_Simu = np.reshape(B_Simu, (1, len(B_Simu)))

    # Get reflectence at reference wavelength
    B_Ref  = hy_obj.get_wave(corr_wave)

    # Get ratio between simulated fresnel and reference value
    RTO = B_Ref / B_Simu[0,:][corr_num]
    RTO[~hy_obj.mask['no_data']] = 0
    RTO[~hy_obj.mask['water']] = 0

    if dimension == 'line':
        # Get RTO for the line
        RTO_line = RTO[index,:]
        RTO_line = np.reshape(RTO_line, (len(RTO_line), 1))

        # Get glint spectrums for the whole line
        correction = RTO_line * B_Simu

        # Correct
        data = data - correction

    elif dimension == 'column':
        RTO_col = RTO[:, index]
        RTO_col = np.reshape(RTO_col, (len(RTO_col), 1))

        # Get glint spectrums for the whole line
        correction = RTO_col * B_Simu

        # Correct
        data = data - correction

    elif (dimension == 'band'):
        # Get glint spectrum for the band
        correction = B_Simu[0,:][index] * RTO

        # Correct
        data = data - correction

    elif dimension == 'chunk':
        RTO_chunk = RTO[y1:y2, x1:x2]
        RTO_chunk = np.reshape(
            RTO_chunk, 
            (
                RTO_chunk.shape[0],
                RTO_chunk.shape[1],
                1
            )
        )
        # Get correction
        correction = RTO_chunk * B_Simu

        # Correct
        data = data - correction

    elif dimension == 'pixels':
        y, x = index

        RTO_pixels = RTO[y, x]
        RTO_pixels = np.reshape(RTO_pixels, (len(RTO_pixels), 1))

        # Get Correction
        correction = RTO_pixels * B_Simu

        # Correct
        data = data - correction

    return data


def glint_spectrum(lut, windspeed, solzen):
    """
    Given the windspeed and solar zenith angle, produce a glint 
    spectrum in reflectance units
    """
    cloudfrac = 0 
    vector = np.array([windspeed,solzen,cloudfrac])

    return lut.data(vector) 


def err(x, corr_val, lut, sol_zen):
    """
    windspeed fit error, assuming the channel "ref_band" is pure glint.
    """
    mdl = glint_spectrum(lut, x, sol_zen)
    er = pow(mdl[lut.ref_band] - corr_val, 2)

    return er


def zenith_refracted(theta, n):
    """
    Find zenith of the outgoing reflected light
    n is the refractive index of water at a specific wavelength
    """
    theta_p = np.degrees(
        np.arcsin(np.sin(np.radians(theta)) / n)
    )

    return theta_p


def fresnel_reflectence(theta, theta_p):
    """
    Uses the fresnel equation to find the 
    percentege of incident light reflected
    """
    theta_rad = np.radians(theta)
    theta_p_rad = np.radians(theta_p)

    return (
        (
            (np.sin(theta_rad - theta_p_rad)**2) 
            / (np.sin(theta_rad + theta_p_rad)**2)
        ) + (
            (np.tan(theta_rad - theta_p_rad)**2) 
            / (np.tan(theta_rad + theta_p_rad)**2)
        )
    )  / 2


def fresnel_spectra(theta, xs, ns):
    """
    Solves for the spectrum of reflected light
    according to fresnels equations
    """
    spectra = []
    for x in xs:
        n = np.interp(x, ns[:,0], ns[:,1])
        theta_p = zenith_refracted(theta, n)
        spectra.append(fresnel_reflectence(theta, theta_p))

    return np.array(spectra)
