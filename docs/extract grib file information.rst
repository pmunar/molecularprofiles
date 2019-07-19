.. _extract grib file information:

Extract grib file information
=============================

In this section we will cover how to extract the information contained in the grib files (either from GDAS or from ECMWF) that we downloaded.
For this purpose we will use the grib_utils.py program.

grib_utils.py
-------------

This program contains a series of functions that will allow us to extract and interact with the information from grib files. 
In order to see what can be done with this program one just needs to type:

.. code-block:: bash

    > grib_utils.py --help

and the help will be displayed.

The first time we want to extract the information from a grib file we simply need to execute:

.. code-block:: bash

    > grib_utils.py -r -g the_grib_file_name.grb -gridstep gridstep -o observatory

where *gridstep* is a value that must be either 0.75 or 1.0, for ECMWF or GDAS, respectively, and *observatory* is *North* or *South*.
If you want to download data for other coordinates other than the ones of the CTA sites, you can specify your coordinates with the -c flag:

.. code-block:: bash

    > grib_utils.py -r -g the_grib_file_name.grb -gridstep gridstep -c lat lon

where *lat* and *lon* are the geographycal latitude and longitude, in degrees, respectively.


This will take some time, specially if the gribfile file is large. It will also **use a huge amount of RAM** (it is due to the pygrib library and how it access the data) so the **user should have a machine with at least 8 Gb of RAM**. The output of this task is a txt file containing the grib file information but in a simpler way, as a text table that is more easily accessible. The final table will contain the following fields or columns:


 -	**Date**: date of the entry, in YYYYMMDD format
 - 	**year**: the year of the format in YYYY format
 -	**month**: the month of the entry
 -	**day**: the day of the entry
 -	**hour**: the hour of the entry, in UTC
 -	**MJD**: the Modified Julian Day of the entry
 -	**P**: the barometric pressure of the entry, in hPa
 -	**Temp**: the temperature of the entry, in Kelvin
 -	**h**: the height of the entry, in meters
 -	**n**: the molecular density of the entry, in particles per cubic centimeter
 -	**U**: the U component of the wind, in m/s
 -	**V**: the V component of the wind, in m/s
 -	**RH**: the relative humidity, in values between 0 and 1
 -	**n_exp**: the molecular density divided by the standard density (Ns) and multiplied by the exponential of the height divided by the standard height (Hs)
 -	**wind_direction**: the direction of the wind, in degrees, where 0 degrees correspond to North direction
 -	**wind_speed**: the module of wind speed, in m/s

Then, this new file can be loaded with the MolecularProfiles class and analyzed.

The grib_utils.py program can also perform other tasks:
If one wants to process several grib files, the process can be speeded up by running it in parallel.
For this purpose, the user can create a file containing the list of grib files to process

.. code-block:: bash

    > ls files*grb >> list_of_girb_files.txt

and then call the program with another flag:

.. code-block:: bash

    > grib_utils.py -r -g list_of_grib_files.txt -p -gridstep 1.0 -o north

where the -p flag activates the parallel processing. If the computer where it is executed has N CPUs, the process will use N-1 CPUs. It will process one file per CPU untill all the files in the list_of_files.txt are processed.

Finally, if we have processed several grib files, either linearly or in parallel, we can merge the results to a single file, for more commodity while analyzing it:

.. code-block:: bash

    > grib_utils.py -m list_of_txt_files.txt

where list_of_txt_files.txt is a file containing all the txt files created in previous steps.

Once the final file is created it can be read in whichever program or app the user prefers. 
