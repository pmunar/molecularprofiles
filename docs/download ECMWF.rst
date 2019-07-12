.. _download ECMWF

download ECMWF
==============

In order to download ECMWF data we need to execute the download_ecmwf script:

.. code-block:: bash

	>download_ecmwf.py date_start date_end latitude longitude

where date_stard and date_end must be in YYYY-MM-DD format and latitude and longitude must be in degrees.

This will prepare and download the data in files that contain one month of data as a maximum amount of data. For time periods longer than one month, separate files will be downloaded.

Once the download has finished you will have your data ready to be processed.