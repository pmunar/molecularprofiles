README grib_utils.py

This is a README file for the grib_utils.py script.

In this script there are a set of functions to work with ECMWF and GDAS grib2 files.

The main functions are:

- date2mjd
- mjd2date
- readgribfile2text
- readgribfile2magic

These functions are described in more detail below:

***********************
* - date2mjd          *
***********************

This function transforms a date in YYY-MM-DD HH:MM:SS format into MDJ

***********************
* - mjd2date          *
***********************

This function transforms a date in MDJ format into YYY-MM-DD HH:MM:SS format


***********************
* - readgribfile2text *
***********************

This function reads a grib2 file and gets the data to write it into a txt file with the following format:
Date, year, month, day, hour, Plevel, T_average, h, n, U, V, RH

where Plevel is a pressure level, h is the height, n the density, U and V the components of the wind and RH the relative humidity.

In order to read the grib2 file and transcript it into txt, the function calls other functions within the script:

	- get_gribfile_variables
	- get_observatory_coordinates
	- get_closest_gridpoint
	- GetAltitudeFromGeopotential

***********************
* - readgribfile2magic*
***********************

This function reads a grib2 file and gets the data to write it into a txt file with the following format, as was used in the MAGIC Collaboration within the MARS software:
year, month, day, hour, index, h, Temperature, U, V, RH

In order to read the grib2 file and transcript it into txt, the function calls other functions within the script:

	- get_gribfile_variables
	- get_observatory_coordinates
	- get_closest_gridpoint
	- GetAltitudeFromGeopotential

The script can be called outside the python interpreter.

Examples:

Example 1:

python grib_utils.py -r my_ecmwf_data.grib north 0.75

Example 2:

python grib_utils.py -mjd 57590

Example 3:

python grib_utils.py -rmagic my_gdas_data.grib north 1.0
