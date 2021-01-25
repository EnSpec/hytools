import os
import json
import glob
import numpy as np

home = os.path.expanduser("~")

config_file = "/home/adam/Dropbox/projects/hytools/configs/ic_config_neon.json"
config_dict = {}

#Only coefficients for good bands will be calculated
config_dict['bad_bands'] =[[300,400],[1337,1430],[1800,1960],[2450,2600]]
config_dict['bad_bands'] =[[300,400],[700,2600]]  # Subset for testing


# Input data settings for NEON
#################################################################
# config_dict['file_type'] = 'neon'
# images= glob.glob("/data/*.h5")
# images.sort()
# config_dict["input_files"] = images

# Input data settings for ENVI
#################################################################
# ''' Only difference between ENVI and NEON settings is the specification
# of the ancillary datasets (ex. viewing and solar geometry). All hytools
# functions assume that the ancillary data and the image date are the same
# size, spatially, and are ENVI formatted files.

# The ancillary parameter is a dictionary with a key per image. Each value
# per image is also a dictionary where the key is the dataset name and the
# value is list consisting of the file path and the band number.

# '''

config_dict['file_type'] = 'envi'
aviris_anc_names = ['path_length','sensor_az','sensor_zn',
                    'solar_az', 'solar_zn','phase','slope',
                    'aspect', 'cosine_i','utc_time']
images= glob.glob("/data/*rfl_v1a_img")
images.sort()
config_dict["input_files"] = images

config_dict["anc_files"] = {}
anc_files = glob.glob("/data/*obs_ort")
anc_files.sort()
for i,image in enumerate(images):
    config_dict["anc_files"][image] = dict(zip(aviris_anc_names,
                                                [[anc_files[i],a] for a in range(len(aviris_anc_names))]))

config_dict['num_cpus'] = len(images)


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
config_dict['export']['subset_waves']  = [660,550,440]
config_dict['export']['output_dir'] = '/data/output/yose_2013/'
config_dict['export']["suffix"] = 'user'

# Specify which correction to apply and order of application
# For no correction provide empty list: []
config_dict["corrections"]  = ['topo','brdf']

#Topographic Correction options
#################################################################
'''
Types supported:
    - 'scsc':  Sonnen et al. 2005
    - 'mod_minneart':
    - 'cosine :
    - 'precomputed'
    - 'none': No correction.

See documentation for options specific to each
correction type
'''
config_dict["topo"] =  {}
config_dict["topo"]['type'] =  'scs+c'
config_dict["topo"]['mask'] = "topo1"

# Precomputed topographic coefficients
# config_dict["topo"]  = {}
# config_dict["topo"]['type'] =  'precomputed'
# config_dict["brdf"]['coeff_files'] =  {}


#BRDF Correction options
#################################################################3
'''
Types supported:
    - 'standard/global': Simple kernel multiplicative correction.
    - 'class': Correction by class.
    - 'precomputed'
    - 'none': No correction.

For precomputed coefficients, input a dictionary with the
file name as the key and pathname to the coeffcients file
as the value.

See documentation for options specific to each
correction type.

If 'by_type' == 'user'
'bins' should be a list of lists, each list the NDVI bounds [low,high]

Object shapes ('h/b'.'b/r') only needed for Li kernels
'''

config_dict["brdf"]  = {}

# Standard BRDF coefficients
############################################
# config_dict["brdf"]['type'] =  'standard'
# config_dict["brdf"]['grouped'] =  True
# config_dict["brdf"]['sample_perc'] = 0.1
# config_dict["brdf"]['geometric'] = 'li_sparse_r'
# config_dict["brdf"]['volume'] = 'ross_thick'
# config_dict["brdf"]['mask'] = "brdf1"

# Precomputed BRDF coefficients
#-------------------------------
# config_dict["brdf"]['type'] =  'precomputed'
# config_dict["brdf"]['coeff_files'] =  {}

#Class BRDF coefficients
#------------------------
config_dict["brdf"]['type'] =  'class'
config_dict["brdf"]['grouped'] =  True
config_dict["brdf"]['geometric'] = 'li_dense_r'
config_dict["brdf"]['volume'] = 'ross_thick'
config_dict["brdf"]["b/r"] = 2.5
config_dict["brdf"]["h/b"] = 2
config_dict["brdf"]['sample_perc'] = 0.1
config_dict["brdf"]['interp_kind'] = 'linear'

# NDVI threshold for where to apply correction
config_dict["brdf"]
config_dict["brdf"]['ndvi_min'] = 0.00
config_dict["brdf"]['ndvi_max'] = 1.0
#Normalize to scene average solar zenith angle
config_dict["brdf"]['solar_zn_norm'] =True
# Buffer around edge of NEON images to remove artifact
# Exclude for non NEON images
config_dict["brdf"]['neon_buffer'] = True
#Mask sensor zenith angle less than this value
config_dict["brdf"]['sensor_zn_min'] = np.radians(0)

# Dynamic NDVI params
#---------------------
config_dict["brdf"]['bin_type'] = 'dynamic'
config_dict["brdf"]['num_bins'] = 18
config_dict["brdf"]['ndvi_bin_min'] = 0.05
config_dict["brdf"]['ndvi_bin_max'] = 1.0
config_dict["brdf"]['ndvi_perc_min'] = 10
config_dict["brdf"]['ndvi_perc_max'] = 95

# Fixed bins specified by user
#------------------------
# config_dict["brdf"]['bin_type'] = 'user'
# config_dict["brdf"]['bins']  = [[0,.15],[.15,1]]

#Wavelength resampling options
#################################################################
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


with open(config_file, 'w') as outfile:
    json.dump(config_dict,outfile,indent=3)














