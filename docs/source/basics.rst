.. _basics:

===========
 Basic use
===========


Loading images
==============

HyTools includes options for loading both ENVI formatted binary files
and NEON AOP HDF files.

.. code-block:: python

   import hytools as ht

   #Load an ENVI file
   envi = ht.open_envi('envi_file.bin')
   envi.load_data()
   
   #Load an NEON HDF image
   neon = ht.open_neon('./neon_file.h5')
   neon.load_data()
   

Reading data
============

There are several ways to read data using a hytools object. One option
is to use one of the 'get' methods:

.. code-block:: python

	wave = neon.get_wave(900)		
	band = neon.get_band(10)
	column = neon.get_column(1)
	line = neon.get_line(234)
	chunk = neon.get_chunk(x1,x2,y1,y2)


Alternatively an iterator can be used to cycle along a
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

Currently writing is only supported for ENVI files, however data from NEON hdf
files can be easy written to ENVI format using builtin functions.

First and ENVI header dictionary needs to be generated to specify the
file size, datatype, interleave and other relevant metadata.

.. code-block:: python

    header_dict = ENVI_header_from_hdf(neon)

In this case we are going to export an RGBI image so we need to update
the number of bands:

.. code-block:: python

   head_dict['bands'] = 4
    
Next we create an ENVI writer object which generates the header and image file
using specification in the header dictionary:

.. code-block:: python

    output_name = './neon.bin'
    writer = writeENVI(output_name,header_dict)

Finally we can write the bands to file. First we retrieve the closest
wavelength to each input wavelength using the ``get_wave()`` method, next
we write the band to the new file with the ``write_band()`` method.

.. code-block:: python

   for band_num,wavelength enumerate([660,550,440,880]):
       wave = neon.get_wave(wavelength)
       writer.write_band(wave,band_num)
   writer.close()
		










