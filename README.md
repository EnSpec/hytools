# HyTools


[![DOI](https://zenodo.org/badge/315419247.svg)](https://zenodo.org/badge/latestdoi/315419247)


HyTools is a python library for processing airborne and spaceborne
imaging spectroscopy data. At it's
core it consists of functions for reading and writing
[ENVI](https://www.l3harrisgeospatial.com/docs/ENVIImageFiles.html)
formatted images and reading [NEON
AOP](https://www.neonscience.org/data-collection/airborne-remote-sensing)
HDF files along with a series of image processing functions including
spectral resampling, topographic, BRDF and glint correction, spectral
transforms, masking and more. We have also created a series of command
line tools which combine these functions and provide a streamlined
workflow for processing images.

For examples see the HyTools basics ipython notebook [here](./examples/hytools_basics_notebook.ipynb). For FlexBRDF correction see the tutorial [here](./examples/FlexBRDF_tutorial.md). An alternative of FlexBRDF correction without using Ray can be checked out [here](./examples/separated_flexbrdf.md). NetCDF (EMIT and AVIRIS-NG) and Geographic Lookup Table (GLT) support are added in the latest version [here](./examples/netcdf_glt.md).

# Installation

To install with pip run:
```bash
pip install hy-tools
```
or 
```bash
python -m pip install git+https://github.com/EnSpec/hytools.git
```
or clone
```bash
git clone https://github.com/EnSpec/hytools.git
```
and install with setuptools
```bash
python setup.py install
```

## Basic usage
```python

import hytools as ht

#Create a HyTools container object
hy_obj = ht.HyTools()

#Read and load ENVI file metadata
hy_obj.read_file('./envi_file')

#Calculate NDVI, retrieves closest wavelength to input wavelength
ir = hy_obj.get_wave(900)
red = hy_obj.get_wave(660)
ndvi = (ir-red)/(ir+red)

#or

# Calculate normalized difference index, NDVI by default
ndvi = hy_obj.ndi()

#Other options for retrieving data
band = hy_obj.get_band(10)
column = hy_obj.get_column(1)
line = hy_obj.get_line(234)
chunk = hy_obj.get_chunk(0,100,0,100)
pixels = hy_obj.get_pixels([102,434],[324,345])

# Create a writer object to write to new file
writer = ht.io.WriteENVI('output_envi',hy_obj.get_header())

#Create an iterator object to cycle though image
iterator = hy_obj.iterate(by = 'line')

# Cycle line by line, read from original data
while not iterator.complete:
   #Read next line
   line = iterator.read_next()

   #Do some calculations.......
   radiance = line * gain + offset

   #Write line to file
   writer.write_line(radiance,iterator.current_line)

writer.close()
```
