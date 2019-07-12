.. _download GDAS

download GDAS
=============

In order to download GDAS data we need to modify the ds083.2_test file that resides in the /molecularprofile/gdas_scripts folder. 

In this file we have to modify the following lines:

* date=201607010000/to/201607150000

where the data format is YYYYMMDDhhmm

*nlat, slat, wlon, elon lines: this lines are the coordinates of our location of interest. As of now, the software only allows for processing single grid point data, meaning that the nlat and slat values must be the same, and wlon and elon must be the same. 

Once this is set we can execute the main download_gdas script:

.. code-block:: bash

    >download_gdas.py -submit ds083.2_test

This will print a lot of text and will send us an email when the data is ready to be downloaded. Alternatively, we can execute the command 

.. code-block:: bash

    >download_gdas.py -get_status requestindex

in order to see if the data is ready or not. requestindex is a unique identifier for your query, that appears at the end of the block of text printed in the first command execution.

Once the data is ready to be retrieved, we can download it with the following command:

.. code-block:: bash

    download_gdas.py -download [RequestIndex]

Once the download has finished you will have your data ready to be processed.