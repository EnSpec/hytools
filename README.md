# HyTools

HyTools is a python library for working with imaging spectroscopy data, with a focus on terrestrial scenes. 
At it's core it consists of a series of functions for reading and writing ENVI-formatted images in addition to 
functionalty for reading NEON and PRISMA-formatted AOP HDF files. Built on top of these functions are a series of higher
level processing tools for data analysis which include spectral resampling, topographic and BRDF correction, spectral transforms,
maskings and more.

## Dependencies
- numpy
- h5py
- gdal

## Basic usage
```python
import hytools as ht

#Read an ENVI file
hyObj = ht.openENVI('envi_file.bin')
#Map image data to numpy memmap object
hyObj.load_data()

#Calculate NDVI, retrieves closest wavelength to input lambda in nm
ir = hyObj.get_wave(900)
red = hyObj.get_wave(660)
ndvi = (ir-red)/(ir+red)

#Other options for retrieving data
band = hyObj.get_band(10)
column = hyObj.get_column(1)
line = hyObj.get_line(234)
chunk = hyObj.get_chunk(x1,x2,y1,y2)

# Create a writer object to write to new file
writer = ht.file_io.writeENVI('output_envi.bin',hyObj.header_dict)

#Create an iterator object to cycle though image
iterator = hyObj.iterate(by = 'line')

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
