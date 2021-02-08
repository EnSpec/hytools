import os
import json
import glob

home = os.path.expanduser("~")

config_file = "/home/chlus/dev_hytools/trait_config_krusty.json"

config_dict = {}
config_dict['file_type'] = 'envi'
config_dict["output_dir"] = '/home/chlus/dev_hytools/data/output/traits/'
config_dict['bad_bands'] =[[300,400],[1337,1430],[1800,1960],[2450,2600]]

# Input data settings for NEON
#################################################################
# config_dict['file_type'] = 'neon'
# images= glob.glob("/data/*.h5")
# images.sort()
# config_dict["input_files"] = images

# Input data settings for ENVI
#################################################################
''' Only differnce between ENVI and NEON settings is the specification
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
images= glob.glob("/home/chlus/dev_hytools/data/yose/*img")
images.sort()
config_dict["input_files"] = images

config_dict["anc_files"] = {}
anc_files = glob.glob("/home/chlus/dev_hytools/data/yose/*ort")
anc_files.sort()
for i,image in enumerate(images):
    config_dict["anc_files"][image] = dict(zip(aviris_anc_names,
                                                [[anc_files[i],a] for a in range(len(aviris_anc_names))]))

config_dict['num_cpus'] = len(images)

# Assign correction coefficients
##########################################################
config_dict['corrections'] = ['brdf']

topo_files = glob.glob("/home/chlus/dev_hytools/data/output/traits/*topo*full.json")
topo_files.sort()
config_dict["topo"] =  dict(zip(images,topo_files))

brdf_files = glob.glob("/home/chlus/dev_hytools/data/output/traits/*brdf*full.json")
brdf_files.sort()
config_dict["brdf"] =  dict(zip(images,brdf_files))

config_dict["resampling"]  = {}
config_dict["resampling"]['type'] =  'cubic'
# config_dict["resampling"]['out_waves'] = []
# config_dict["resampling"]['out_fwhm'] = []


# Define trait coefficients
##########################################################
config_dict["trait_models"]  = glob.glob('/home/chlus/dev_hytools/data/traits/*.json')

with open(config_file, 'w') as outfile:
    json.dump(config_dict,outfile)

















