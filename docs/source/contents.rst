About HyTools
=====================

HyTools is a python library for processing airborne and spaceborne
imaging spectroscopy data, with a focus on terrestrial scenes. At it's
core it consists of functions for reading and writing `ENVI
<https://www.l3harrisgeospatial.com/docs/ENVIImageFiles.html>`_
formatted images and reading `NEON AOP
<https://www.neonscience.org/data-collection/airborne-remote-sensing>`_
HDF files along with a series of image processing functions including
spectral resampling, topographic and BRDF correction, spectral
transforms, masking and more. We have also created a series of command
line tools which combine these functions and provide a streamlined
workflow for processing images.

Examples
--------

BRDF correction
~~~~~~~~~~~~~~~

.. raw:: html
	 
  <embed>

  <link rel="stylesheet" href="/_static/css/slider.css">
  <script src="/_static/js/slider.js" type="text/javascript" ></script>

  <div id="slider" class="beer-slider" data-beer-label="">
     <img src="/_static/images/research/3d_rgb.jpg" alt="">
     <div class="beer-reveal" data-beer-label="">
        <img src="/_static/images/research/rgb_rgb.jpg" alt="">
     </div>
   </div>

   <script type="text/javascript">
       new BeerSlider(document.getElementById('slider'));
   </script>

   </embed>

  

.. image:: brdf_before_after.png

Topographic correction
~~~~~~~~~~~~~~~~~~~~~~
.. image:: topo_correct.gif	   

	   
