.. molecularprofiles documentation master file, created by
   sphinx-quickstart on Thu May 30 16:30:40 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to molecularprofiles' documentation!
============================================

This is molecularprofiles, a Python package that will help in the analysis of molecular profile information obtained from global data assimilation systems, like GFS, GDAS or ECMWF.

This library works only with grib file format and is specifically designed for the analysis of molecular content above the CTA sites, at El Roque de los Muchachos in the island of La Palma, and at Paranal in Chile. Other locations can be used as well, but some functions may not work as intended.

The library helps the user to perform several tasks:

	* Download the global data assimilation system data
	* Extract these data and transform it from grib to text format
	* Analyze these data and make plots in order to visualize the content
	* Generate input cards for the CORSIKA software for simulation of air showers
	* It also has some utilities to compute times and dates


* Code: https://github.com/pmunar/molecularprofiles
* Issues: https://github.com/pmunar/molecularprofiles/issues
* Documentation: http://molecularprofiles.readthedocs.org/
* Contact: pere.munar@uab.cat

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   setup
   download GDAS
   download ECMWF
   extract grib file information
   the MolecularProfile class
   tutorial




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
