.. _download ECMWF

download ECMWF
==============

In order to download ECMWF data first we need to `register at the ECMWF site <https://apps.ecmwf.int/registration/>`__ and follow `these instructions <https://confluence.ecmwf.int/display/WEBAPI/Accessing+ECMWF+data+servers+in+batch#AccessingECMWFdataserversinbatch-key>`__ in order to set the appropriate key and environment for the ECMWF api to work.

Once we have these steps covered, we need to execute the download_ecmwf script from within our working directory:

.. code-block:: bash

	>download_ecmwf.py date_start date_end latitude longitude

where date_stard and date_end must be in YYYY-MM-DD format and latitude and longitude must be in degrees.

This will prepare and download the data in files that contain one month of data as a maximum amount of data. For time periods longer than one month, separate files will be downloaded.

Once the download has finished you will have your data ready to be processed.