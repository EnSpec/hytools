.. _basics:

===========
 Basic use
===========


Loading images
==============

HyTools includes options for loading both ENVI formatted binary file
and NEON AOP HDF files.

.. code-block:: python

   import hytools as ht

   #Load an ENVI file
   envi = ht.open_envi('envi_file.bin')

   #Load an NEON HDF image
   neon = ht.open_neon('./neon_file.h5')


Reading data
============

There are several ways to read data using a hytools object. 

.. code-block:: python

	band = hy_obj.get_band(10)
	column = hy_obj.get_column(1)
	line = hy_obj.get_line(234)
	chunk = hy_obj.get_chunk(x1,x2,y1,y2)	
		

HyTools also includes and iterator class which iterates along a
specified axis of the dataset either by line, column, band or
chunk. This is useful for cycle through and image, applying
a function/algorithm to an image and then writing to a file.

.. code-block:: python
		
   iterator = hy_obj.iterate(by = 'line')

Next cycle through the image line by line until complete:

.. code-block:: python
		
   while not iterator.complete:  
       line = iterator.read_next() 
       radiance = line * gain + offset


Writing data
============

Writing is only supported for ENVI files





