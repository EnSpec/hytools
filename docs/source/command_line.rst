.. _command:

====================
 Command line tools
====================


image_correct.py
================

The 'image_correct' script is a general purpose tool that utilizes
topographic, BRDF and wavelength resampling fuctions to modify/correct
an image.

Configuration file
------------------

Image correction options are specified using a JSON file. Example
configuration files can be found in the github repository here. We
have also supplied a script to automatically generate a JSON
configuration file.



trait_estimate.py
=================

The 'trait_estimate.py' script is a tool to generate maps of canopy
foliar traits given a set of model parameters. Currently only PLSR
models are supported, however support for other statiscal modeling
frameworks like Gaussian process regression is currently underway.



