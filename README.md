# HyTools

HyTools is a python library for processing airborne and spaceborne
imaging spectroscopy data, with a focus on terrestrial scenes. At it's
core it consists of functions for reading and writing
[ENVI](https://www.l3harrisgeospatial.com/docs/ENVIImageFiles.html)
formatted images and reading [NEON
AOP](https://www.neonscience.org/data-collection/airborne-remote-sensing)
HDF files along with a series of image processing functions including
spectral resampling, topographic and BRDF correction, spectral
transforms, masking and more. We have also created a series of command
line tools which combine these functions and provide a streamlined
workflow for processing images.

For complete documentation see:
[hytools.readthedocs.io](https://hytools.readthedocs.io/en/latest/contents.html)


## Dependencies
- h5py
- matplotlib
- numpy
- pandas
- scikit-learn
- scipy

# Installation
To install run:

```python
python setup.py install
```
## Basic usage
```python

import hytools as ht

#Create a HyTools container object
envi = ht.HyTools()

#Read and load file metadata
envi.read_data('./envi_file.bin',file_type= 'envi')

#Calculate NDVI, retrieves closest wavelength to input wavelength
ir = hy_obj.get_wave(900)
red = hy_obj.get_wave(660)
ndvi = (ir-red)/(ir+red)

#Other options for retrieving data
band = hy_obj.get_band(10)
column = hy_obj.get_column(1)
line = hy_obj.get_line(234)
chunk = hy_obj.get_chunk(x1,x2,y1,y2)

# Create a writer object to write to new file
writer = ht.io.WriteENVI('output_envi.bin',hy_obj.header_dict)

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
