.. _download ECMWF

download ECMWF
==============

In order to download ECMWF data first we need to `register at the ECMWF site <https://apps.ecmwf.int/registration/>`__ and follow `these instructions <https://confluence.ecmwf.int/display/WEBAPI/Accessing+ECMWF+data+servers+in+batch#AccessingECMWFdataserversinbatch-key>`__ in order to set the appropriate key and environment for the ECMWF api to work.

Once we have these steps covered, we need to execute the download_ecmwf script from within our working directory:

.. code-block:: bash

	>download_ecmwf.py date_start date_end latitude longitude

where date_stard and date_end must be in YYYY-MM-DD format and latitude and longitude must be in degrees.

This will prepare and download the data in files that will contain by default 7 days of data. The user can specify other grouping, by introducing a new argument, called *days*:

.. code-block:: bash

	>download_ecmwf.py date_start date_end latitude longitude -days 10

In this example, now the files would contain 10 days of data each one. The final file may not contain this amount of data, because it will contain up to the end date selected.

** It is strongly recommended to keep the amount of dates in each data file low (7 days is enough) since it will improve the RAM usage in the reading step, with grib_utils.py**

Once the download has finished you will have your data ready to be processed.