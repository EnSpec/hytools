# -*- coding: utf-8 -*-
"""
HyTools:  Hyperspectral image processing library

Copyright (C) 2021 University of Wisconsin

Authors: Adam Chlus, Zhiwei Ye, Philip Townsend.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

JSON format for PLSR trait models v0.1

This version assumes that all transforms are applied using all model wavelenghts.

TODO: Develop standard codes for spectrometer, both airborne/spaceborne and field.
TODO: Allow for more options for specifying wavelength subsets

"""


import json

model_dict = {}

# Metadata
#####################################
'''
    trait : Trait name (str)
    units : Trait units (str)
    description: Model description (str)
    wavelength_units: Wavelength units (str)
    wavelengths : Model wavelengths (list)
        Only wavelengths used in the model should be
        included in the list of wavelengths.
    fwhm : Model fwhm (list)
    type : Model type (str)
'''
model_dict["name"] = ''
model_dict["units"] = ''
model_dict["description"] = ''
model_dict["wavelength_units"] = ''
model_dict["wavelengths"] = []
model_dict["fwhm"] = []
model_dict["spectrometer"] = ''
model_dict["type"] =  ''

# Diagnostics
#####################################
'''Currently the only required diagnostics are 'min'
and 'max', these are the min and max values of the
dataset used to build the model and are used to generate
the data range mask, which identifies pixels with predictions
outside of the model dataset range.

'''
model_dict["model_diagnostics"] = {}
model_dict["model_diagnostics"]["rmse"] = 0.0
model_dict["model_diagnostics"]["r_squared"] = 0.0
model_dict["model_diagnostics"]["min"] = 0.0
model_dict["model_diagnostics"]["max"] = 0.0

# Model
#####################################
'''
transform:  List of transforms to be applied in order of application.
            Options:
                - 'vector': vector norm using np.linalg.norm
                - 'mean' : Normalize to mean
                - 'absorb' : log(1/R)
            Examples:
                ['vector','absorb']
            Empty list for no transforms ([])
coefficients: List of lists, sublists are the coefficients for
                model iterations.

intercepts : Permuted model intercepts (list)
components : Number of model component (int)
'''

model_dict['model'] = {}
model_dict['model']["components"] =  0
model_dict['model']["transform"] = ['mean']
model_dict['model']["intercepts"] = []
model_dict['model']["coefficients"] =[[],[]]

model_path = '*.json'
with open(model_path, 'w') as outfile:
    json.dump(model_dict,outfile)





















