import numpy as np
from numba import jit


@jit(nopython=True)
def gaussian(x,mu,fwhm):
    """Return a gaussian distribution.

    Parameters
    ----------
    x : Numpy array of values along which to generate gaussian.
    mu : Mean of the gaussian function.
    fwhm : Full width half maximum.

    Returns
    -------
    Numpy array of gaussian along input range.
    """

    c = fwhm/(2* np.sqrt(2*np.log(2)))
    return np.exp(-1*((x-mu)**2/(2*c**2)))

def resample_coeff_single(srcWaves, dstWaves,dstFWHMs):
    """ Return a set of coeffiencients for spectrum resampling

    Given a set of source wavelengths, destination wavelengths and FWHMs this
    function caculates the relative contribution or each input wavelength
    to the output wavelength. It assumes that output
    response functions follow a gaussian distribution.

    Parameters
    ----------
    srcWaves : List of source wavelength centers.
    dstWaves : List of destination wavelength centers.
    dstFWHMs : List of destination full width half maxes.

    Returns
    -------
    m x n matrix of coeffiecients, where m is the number of source wavelengths
    and n is the number of destination wavelengths.
    """

    # For each destination band calculate the relative contribution
    # of each wavelength to the band response at source resolution
    dstMatrix = []
    #oneNM = np.arange(280,2600)
    for dstWave,dstFWHM in zip(dstWaves,dstFWHMs):
        a =  gaussian(srcWaves -.5,dstWave,dstFWHM)
        b =  gaussian(srcWaves +.5,dstWave,dstFWHM)
        areas = (a +b)/2
        dstMatrix.append(np.divide(areas,np.sum(areas)))
    dstMatrix = np.array(dstMatrix)

    return dstMatrix.T


@jit(nopython=True)
def resample_coeff(srcWaves,srcFWHMs,dstWaves,dstFWHMs, spacing = 1):
    """ Return a set of coeffiencients for spectrum resampling

    Given a set of source and destination wavelengths and FWHMs this
    function caculates the relative contribution or each input wavelength
    to the output wavelength. It assumes that both input and output
    response functions follow a gaussian distribution.

    Parameters
    ----------
    srcWaves : List of source wavelength centers.
    srcFWHMs : List of source full width half maxes.
    dstWaves : List of destination wavelength centers.
    dstFWHMs : List of destination full width half maxes.
    spacing : resolution at which to model the spectral resposnse functions

    Returns
    -------
    m x n matrix of coeffiecients, where m is the number of source wavelengths
    and n is the number of destination wavelengths.
    """

    # For each destination band calculate the relative contribution
    # of each wavelength to the band response at 1nm resolution

    all_waves = np.concatenate((srcWaves,dstWaves))

    min_spectrum = all_waves.min()/100 *100
    max_spectrum = 100 + all_waves.max()//100 *100
    oneNM = np.arange(min_spectrum,max_spectrum,spacing)

    dstMatrix = np.zeros((len(oneNM),len(dstWaves)))

    for i,(dstWave,dstFWHM) in enumerate(zip(dstWaves,dstFWHMs)):
        a =  gaussian(oneNM -.5,dstWave,dstFWHM)
        b =  gaussian(oneNM +.5,dstWave,dstFWHM)
        areas = (a +b)/2
        dstMatrix[:,i] = np.divide(areas,np.sum(areas))

    # For each source wavelength generate the gaussion response
    # function at 1nm resolution
    srcMatrix =
    for srcWave,srcFWHM in zip(srcWaves,srcFWHMs):
        srcMatrix.append( gaussian(oneNM ,srcWave,srcFWHM))
    srcMatrix = np.array(srcMatrix)

    # Calculate the relative contribution of each source response function
    ratio =  srcMatrix/srcMatrix.sum(axis=0)
    ratio[np.isnan(ratio)] = 0
    ratio2 = np.einsum('ab,cb->acb',ratio,dstMatrix)

    # Calculate the relative contribution of each input wavelength
    # to each destination wavelength
    coeffs = np.trapz(ratio2)

    return coeffs


srcWaves = np.arange(400,2500,5)
dstWaves = np.arange(400,2500,5)
dstFWHMs = np.array([5 for x in dstWaves])
srcFWHMs = np.array([5 for x in srcWaves])

test =resample_coeff(srcWaves,srcFWHMs,
        dstWaves,dstFWHMs, spacing = 1)










def est_fwhm(hyObj, dstWaves, dstFWHMs):
    """
    Acquire source wavelength information from source dataset, and estimate the list of destination full width half maxes if they are not specify by the user.

    Parameters
    ----------
    hyObj : hyTools data object

    dstWaves : List of destination wavelength centers.

    dstFWHMs : List of destination full width half maxes.


    Returns
    -------
    srcWaves : List of source wavelength centers.
    srcFWHMs : List of source full width half maxes.
    dstFWHMs : List of destination full width half maxes.
    """

    srcWaves = hyObj.wavelengths
    srcFWHMs = hyObj.fwhm

    if dstFWHMs is None:
      gap = 0.5 * (dstWaves[1:] - dstWaves[:-2])
      #print(gap,gap[1:], gap[:-1],gap[-1])
      dstFWHMs_middle = gap[1:] + gap[:-1]
      #print(dstFWHMs_middle)
      dstFWHMs = np.append(np.append(gap[0]*2, dstFWHMs_middle), gap[-1]*2)

    return (srcWaves, srcFWHMs, dstFWHMs)


def est_transform_matrix(srcWaves, dstWaves, srcFWHMs, dstFWHMs, method_code):

    if method_code==0:
        coeffs =  resample_coeff_single(srcWaves,dstWaves,dstFWHMs)
    elif method_code==1:
        coeffs =  resample_coeff(srcWaves,srcFWHMs,dstWaves,dstFWHMs, spacing = 1)
    else:
        coeffs =  matrix_inverse(srcWaves,srcFWHMs,dstWaves,dstFWHMs)

    return coeffs



def resample_img(hyObj, output_name, dstWaves, method="single_FWHM",  dstFWHMs = None):
    """

    Parameters
    ----------
    hyObj : hyTools data object

    output_name: str
        Path name for resampled file.

    dstWaves
        List of destination wavelength centers.

    dstFWHMs
        List of destination full width half maxes. Default is None.

    method
        single_FWHM (default)
            interpolation without srcFWHMs (List of source full width half maxes)
        two_FWHM
            interpolation with both srcFWHMs and dstFWHMs
        two_FWHM_minnorm
            Minimum-norm least squares problem (underdetermined case while resampling to 1nm) solved by pseudoinverse,
            with both srcFWHMs and dstFWHMs

        All methods require the list of source wavelength centers (srcWaves), which should be stored in hyObj.

    Returns
    -------
    None
    """

    # Dictionary of all methods
    methodDict = {"single_FWHM": 0,
                            "two_FWHM":1,
                            "two_FWHM_minnorm": 2
                           }

    srcWaves, srcFWHMs, dstFWHMs =  est_fwhm(hyObj, dstWaves, dstFWHMs)

    band_mask = [x==1 for x in hyObj.bad_bands]

    coeffs =  est_transform_matrix(srcWaves, dstWaves, srcFWHMs, dstFWHMs, methodDict[method])

    # update destination bad band list
    new_badband = np.dot(band_mask,coeffs)
    new_badband = (new_badband > 0.9).astype(np.uint8)

    # update header dictionary of the destination image
    new_headerDict = hyObj.header_dict
    new_headerDict["bands"] = len(dstWaves)
    new_headerDict["wavelength"] = '{'+','.join(['%g' % x for x in dstWaves])+'}'
    new_headerDict["fwhm"]  = '{'+','.join(['%g' % x for x in dstFWHMs])+'}'
    new_headerDict["bbl"] = '{'+','.join([str(x) for x in new_badband])+'}'
    new_bandnames = ['Band_'+str(x) for x in dstWaves]
    new_headerDict["band names"] = '{'+ ','.join(new_bandnames) +'}'


    if  hyObj.file_type == "ENVI":
        writer = writeENVI(output_name,new_headerDict)
    elif hyObj.file_type == "HDF":
        writer = None
    else:
        print("ERROR: File format not recognized.")

    iterator = hyObj.iterate(by = 'chunk')

    while not iterator.complete:
        chunk = iterator.read_next()
        #resampled_chunk = np.dot(coeffs,chunk[:,:,band_mask])
        resampled_chunk = np.dot(chunk[:,:,:], coeffs)
        #resampled_chunk[ resampled_chunk <1000 ] = hyObj.no_data
        writer.write_chunk(resampled_chunk,iterator.current_line,iterator.current_column)

    writer.close()

    return 1





