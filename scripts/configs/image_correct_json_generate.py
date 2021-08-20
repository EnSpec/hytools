'''Template script for generating image_correct configuration JSON files.

    These settiing are meant only as an example, are not appropriate for
    all situations and may need to be adjusted
'''

import os
import json
import glob
import numpy as np

#Output path for configuration file
config_file = "/data1/temp/ht_test/ic_test.json"
config_dict = {}

#Only coefficients for good bands will be calculated
config_dict['bad_bands'] =[[300,400],[1337,1430],[1800,1960],[2450,2600]]
config_dict['bad_bands'] =[[300,400],[900,2600]]  # Subset for testing

# Input data settings for NEON
#################################################################
# config_dict['file_type'] = 'neon'
# images= glob.glob("/data1/temp/ht_test/*.h5")
# images.sort()
# config_dict["input_files"] = images

# Input data settings for ENVI
#################################################################

''' Only difference between ENVI and NEON settings is the specification
of the ancillary datasets (ex. viewing and solar geometry). All hytools
functions assume that the ancillary data and the image date are the same
size, spatially, and are ENVI formatted files.

The ancillary parameter is a dictionary with a key per image. Each value
per image is also a dictionary where the key is the dataset name and the
value is list consisting of the file path and the band number.

'''

config_dict['file_type'] = 'envi'
aviris_anc_names = ['path_length','sensor_az','sensor_zn',
                    'solar_az', 'solar_zn','phase','slope',
                    'aspect', 'cosine_i','utc_time']
images= glob.glob("/data1/temp/ht_test/*img")
images.sort()
config_dict["input_files"] = images

config_dict["anc_files"] = {}
anc_files = glob.glob("/data1/temp/ht_test/*ort")
anc_files.sort()
for i,image in enumerate(images):
    config_dict["anc_files"][image] = dict(zip(aviris_anc_names,
                                                [[anc_files[i],a] for a in range(len(aviris_anc_names))]))

# Water masks for glint corrections
mask_files = glob.glob("/data1/temp/ht_test/*mask")
mask_files.sort()
for i, image in enumerate(images):
	config_dict["anc_files"][image]['water'] = [mask_files[i], 0]


# Export settings
#################################################################
''' Options for subset waves:
    1. List of subset wavelenths
    2. Empty list, this will output all good bands, if resampler is
    set it will also resample.
    - Currently resampler cannot be used in conjuction with option 1
'''

config_dict['export'] = {}
config_dict['export']['coeffs']  = True
config_dict['export']['image']  = True
config_dict['export']['masks']  = True
config_dict['export']['subset_waves']  = [660,550,440,850]
config_dict['export']['output_dir'] = "/data1/temp/ht_test/"
config_dict['export']["suffix"] = 'brdf'

#Corrections
#################################################################
''' Specify correction(s) to be applied, corrections will be applied
in the order they are specified.

Options include:

    ['topo']
    ['brdf']
    ['topo','brdf']
    ['brdf','topo']
    [] <---Export uncorrected images

'''

config_dict["corrections"]  = ['topo','brdf']

#Topographic Correction options
#################################################################
'''
Types supported:
    - 'cosine'
    - 'c'
    - 'scs
    - 'scs+c'
    - 'mod_minneart'
    - 'precomputed'

Apply and calc masks are only needed for C and SCS+C corrections. They will
be ignored in all other cases and correction will be applied to all
non no-data pixels.

'c_fit_type' is only applicable for the C or SCS+C correction type. Options
include 'ols' or 'nnls'. Choosing 'nnls' can limit overcorrection.

For precomputed topographic coefficients 'coeff_files' is a
dictionary where each key is the full the image path and value
is the full path to coefficients file, one per image.
'''
config_dict["topo"] =  {}
config_dict["topo"]['type'] =  'scs+c'
config_dict["topo"]['calc_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
                                             'min': 0.1,'max': 1.0}],
                                    ['ancillary',{'name':'slope',
                                                  'min': np.radians(5),'max':'+inf' }],
                                    ['ancillary',{'name':'cosine_i',
                                                  'min': 0.12,'max':'+inf' }],
                                    ['cloud',{'method':'zhai_2018',
                                              'cloud':True,'shadow':True,
                                              'T1': 0.01,'t2': 1/10,'t3': 1/4,
                                              't4': 1/2,'T7': 9,'T8': 9}]]

config_dict["topo"]['apply_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
                                             'min': 0.1,'max': 1.0}],
                                    ['ancillary',{'name':'slope',
                                                  'min': np.radians(5),'max':'+inf' }],
                                    ['ancillary',{'name':'cosine_i',
                                                  'min': 0.12,'max':'+inf' }]]
config_dict["topo"]['c_fit_type'] = 'nnls'

# config_dict["topo"]['type'] =  'precomputed'
# config_dict["brdf"]['coeff_files'] =  {}

#BRDF Correction options
#################################################################3
'''
Types supported:
    - 'universal': Simple kernel multiplicative correction.
    - 'local': Correction by class. (Future.....)
    - 'flex' : Correction by NDVI class
    - 'precomputed' : Use precomputed coefficients

If 'bin_type' == 'user'
'bins' should be a list of lists, each list the NDVI bounds [low,high]

Object shapes ('h/b','b/r') only needed for Li kernels.

For precomputed topographic coefficients 'coeff_files' is a
dictionary where each key is the full the image path and value
is the full path to coefficients file, one per image.
'''

config_dict["brdf"]  = {}

# Options are 'line','scene', or a float for a custom solar zn
# Custom solar zenith angle should be in radians
config_dict["brdf"]['solar_zn_type'] ='scene'

# Universal BRDF config
#----------------------
# config_dict["brdf"]['type'] =  'universal'
# config_dict["brdf"]['grouped'] =  True
# config_dict["brdf"]['sample_perc'] = 0.1
# config_dict["brdf"]['geometric'] = 'li_dense_r'
# config_dict["brdf"]['volume'] = 'ross_thick'
# config_dict["brdf"]["b/r"] = 2.5
# config_dict["brdf"]["h/b"] = 2
# config_dict["brdf"]['calc_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
#                                             'min': 0.1,'max': 1.0}]]
# config_dict["brdf"]['apply_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
#                                             'min': 0.1,'max': 1.0}]]
# config_dict["brdf"]['diagnostic_plots'] = True
# config_dict["brdf"]['diagnostic_waves'] = [440,550,660,850]

#----------------------
# ## Flex BRDF configs
# ##------------------
config_dict["brdf"]['type'] =  'flex'
config_dict["brdf"]['grouped'] =  True
config_dict["brdf"]['geometric'] = 'li_dense_r'
config_dict["brdf"]['volume'] = 'ross_thick'
config_dict["brdf"]["b/r"] = 2.5
config_dict["brdf"]["h/b"] = 2
config_dict["brdf"]['sample_perc'] = 0.1
config_dict["brdf"]['interp_kind'] = 'linear'
config_dict["brdf"]['calc_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
                                              'min': 0.1,'max': 1.0}],
                                    ['kernel_finite',{}],
                                    ['ancillary',{'name':'sensor_zn',
                                                  'min':np.radians(2),'max':'inf' }],
                                    ['neon_edge',{'radius': 30}],
                                    ['cloud',{'method':'zhai_2018',
                                              'cloud':True,'shadow':True,
                                              'T1': 0.01,'t2': 1/10,'t3': 1/4,
                                              't4': 1/2,'T7': 9,'T8': 9}]]

config_dict["brdf"]['apply_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
                                              'min': 0.05,'max': 1.0}]]

# ## Flex dynamic NDVI params
config_dict["brdf"]['bin_type'] = 'dynamic'
config_dict["brdf"]['num_bins'] = 18
config_dict["brdf"]['ndvi_bin_min'] = 0.05
config_dict["brdf"]['ndvi_bin_max'] = 1.0
config_dict["brdf"]['ndvi_perc_min'] = 10
config_dict["brdf"]['ndvi_perc_max'] = 95

# ## Flex fixed bins specified by user
# config_dict["brdf"]['bin_type'] = 'user'
# config_dict["brdf"]['bins']  = [[0.1,.25],[.25,.75],[.75,1]]
# ##-----------------

## Precomputed BRDF coefficients
##------------------------------
# config_dict["brdf"]['type'] =  'precomputed'
# config_dict["brdf"]['coeff_files'] =  {}
##------------------------------

#----------------------
# ## Glint correction configs
# ##------------------
config_dict["glint"]  = {}
# Included types are: hedley, hochberg, and gao
# Glint correction options
'''
Types supported:
    - hochberg
    - hedley
    - gao
'''
config_dict['glint']['type'] = 'hochberg'

# Reference band to guide corrections (nm)
'''
Common reference bands include:
    - 860nm (NIR)
    - 1650nm (SWIR)
    - 2190nm (SWIR)
'''
config_dict['glint']['correction_wave'] = 1650 

# If glint correction is Hedley method, you need column, row chnks for homogenous deep water.
'''
The Hedley-specific config would be in the form of:
[ImagePath]: [y1, y2, x1, x1]

e.g.:
config_dict["glint"]["deep_water_sample"] = {
     "/path_to_image1": [
        137, 574, 8034, 8470
     ],
     "/path_to_image2": [
        48, 393, 5780, 5925
     ],
}

#Wavelength resampling options
##############################
'''
Types supported:
   - 'gaussian': needs output waves and output FWHM
   - 'linear', 'nearest', 'nearest-up',
      'zero', 'slinear', 'quadratic','cubic': Piecewise
      interpolation using Scipy interp1d

config_dict["resampler"] only needed when resampling == True
'''
config_dict["resample"]  = False
# config_dict["resampler"]  = {}
# config_dict["resampler"]['type'] =  'cubic'
# config_dict["resampler"]['out_waves'] = []
# config_dict["resampler"]['out_fwhm'] = []

# Remove bad bands from output waves
# for wavelength in range(450,660,100):
#     bad=False
#     for start,end in config_dict['bad_bands']:
#         bad = ((wavelength >= start) & (wavelength <=end)) or bad
#     if not bad:
#         config_dict["resampler"]['out_waves'].append(wavelength)


config_dict['num_cpus'] = len(images)

with open(config_file, 'w') as outfile:
    json.dump(config_dict,outfile,indent=3)
