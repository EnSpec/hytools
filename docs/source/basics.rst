.. _basics:

===========
 Basic use
===========


Loading images
==============

HyTools includes options for loading both ENVI formatted binary files, NASA NetCDF files,
and NEON AOP HDF files.

.. code-block:: python

   import hytools as ht

   #Create a HyTools container object 
   envi = ht.HyTools()

   #Read and load file metadata
   envi.read_data('./envi_file.bin',file_type= 'envi')

For reading NEON data the process is the same:

.. code-block:: python
		
   #Load an NEON HDF image
   neon = ht.HyTools()
   neon.read_data("./neon_file.h5",'neon')
     
  
Reading data
============

There are several ways to read data using a :class:`~hytools.base.HyTools` object. One option
is to use one of the 'get' methods:

.. code-block:: python

	wave = neon.get_wave(900)		
	band = neon.get_band(10)
	column = neon.get_column(1)
	line = neon.get_line(234)
	chunk = neon.get_chunk(x1,x2,y1,y2)
	pixels = neon.get_pixels([0,1,2],[3,4,5])
	
We can also retrieve masked data, where a binary mask is used to
return a subset of the data. Currently masking only works using the
:meth:`~hytools.base.HyTools.get_band` or
:meth:`~hytools.base.HyTools.get_wave` methods. First we need to
generate a mask, which can be done using the
:meth:`~hytools.base.HyTools.gen_mask` method.

.. code-block:: python

	# NDVI masking function	
	def masker(hy_obj):
	    ir = hy_obj.get_wave(900)
	    red = hy_obj.get_wave(660)
	    ndvi = (ir-red)/(ir+red)
	    return ndvi > .5	

	# Generate mask
	neon.gen_mask(masker)

	# Retrieve pixels where mask is True
	pixels = neon.get_band(100, mask_values = True)
	

Alternatively an :class:`~hytools.base.Iterator` can be used to cycle along a
specified axis of the dataset either by line, column, band or
chunk. This is useful for cycling through and image, applying
a function/algorithm and then writing to a file.

.. code-block:: python
		
   iterator = hy_obj.iterate(by = 'line')

Next cycle through the image line by line until complete:

.. code-block:: python
		
   while not iterator.complete:  
       line = iterator.read_next() 


Writing data
============

Currently writing is only supported for ENVI files and NetCDF files, however data from
NEON hdf files can be easy written to ENVI format using builtin
functions.

First an ENVI header dictionary needs to be generated to specify the
file size, datatype, interleave and other relevant metadata. This is
done using the :func:`~hytools.io.envi.envi_header_from_hdf` function.

.. code-block:: python

    header_dict = envi_header_from_hdf(neon)

In this case we are going to export an RGBI image so we need to update
the number of bands:

.. code-block:: python

   head_dict['bands'] = 4
    
Next we create an :class:`~hytools.io.envi.WriteENVI` object which
generates the header and image file using the specifications in the
header dictionary:

.. code-block:: python

    output_name = './neon.bin'
    writer = WriteENVI(output_name,header_dict)

Finally we can write the bands to file. First we retrieve the closest
wavelength to each input wavelength using the
:meth:`~hytools.base.HyTools.get_wave` method, next we write the band
to the new file with the :meth:`~hytools.io.envi.WriteENVI.write_band`
method.

.. code-block:: python

   for band_num,wavelength enumerate([660,550,440,880]):
       wave = neon.get_wave(wavelength)
       writer.write_band(wave,band_num)
   writer.close()
		










