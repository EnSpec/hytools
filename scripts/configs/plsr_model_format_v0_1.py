'''JSON format for PLSR trait models v0.1

This version assumes that all transforms are applied using all model wavelenghts.

TODO: Develop standard codes for spectrometer, both airborne/spaceborne and field.
TODO: Allow for more options for specifying wavelength subsets

'''

model_dict = {}

# Metadata
#####################################
'''
- Only wavelengths used in the model should be
    included in the list of wavelengths.
'''
model_dict["name"] = str
model_dict["units"] = str
model_dict["description"] = str
model_dict["wavelength_units"] =  str
model_dict["wavelengths"] = list
model_dict["fwhm"] = list
model_dict["spectrometer"] = str
model_dict["type"] =  'plsr'

# Diagnostics
#####################################
'''Currently the only required diagnostics are 'min'
and 'max', these are the min and max values of the
dataset used to build the model and are used to generate
the data range mask, which identifies pixels with predictions
outside of the model dataset range.
'''

model_dict["model_diagnostics"] = {}
model_dict["model_diagnostics"]["rmse"] = float
model_dict["model_diagnostics"]["r_squared"] = float
model_dict["model_diagnostics"]["min"] = float
model_dict["model_diagnostics"]["max"] = float

# Model
#####################################
'''
transform:  List of transforms to be applied in order of application.
            Options:
                - 'vector': vector norm using np.linalg.norm
                - 'mean' : Normalize to mean
                - 'absorb' : log(1/R)
                - .......
            Examples:
                ['vector','absorb']
            Empty list for no transforms ([])

coefficients: List of lists, sublists are the coefficients for
                mode iterations.
'''

model_dict['model'] = {}
model_dict['model']["components"] =  int
model_dict['model']["transform"] = list
model_dict['model']["intercepts"] = list
model_dict['model']["coefficients"] list

model_path = '/'
with open(model_path, 'w') as outfile:
    json.dump(model_dict,outfile)



















