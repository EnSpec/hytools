# HyTools

HyTools is a python library for processing airborne and spaceborne
imaging spectroscopy data, with a focus on terrestrial scenes. At it's
core it consists of functions for reading and writing
[ENVI](https://www.l3harrisgeospatial.com/docs/ENVIImageFiles.html)
formatted images and reading [NEON
AOP](https://www.neonscience.org/data-collection/airborne-remote-sensing)
HDF files along with series of image processing functions including
spectral resampling, topographic and BRDF correction, spectral
transforms, masking and more. We have also created a series of command
line tools which combine these functions and provide a more
streamlined workflow for processing images.

For complete documentation see: [hytools.readthedocs.io](https://hytools.readthedocs.io)


## Dependencies
- numpy
- h5py
- scipy

# Installation
To install run:

```python
python setup.py install
```

## Basic usage
```python
import hytools as ht

#Read an ENVI file
hy_obj = ht.open_envi('envi_file.bin')
#Map image data to numpy memmap object
hy_obj.load_data()

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
