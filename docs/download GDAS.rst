.. _download GDAS

download GDAS
=============

If we want to download GDAS data the first thing we need to do is to `register at the Research Data Archive site <https://rda.ucar.edu/index.html?hash=data_user&action=register>`__ . Then you need to create a text file containing a single line with

.. code-block:: bash

    your_email, your_password

name it rdamspw.txt and place it within the folder $MOLECULARPROFILES_DIR/molecularprofile/gdas_scripts

Then we need to copy the ds083.2_test file, that resides in the $MOLECULARPROFILES_DIR/molecularprofiles/gdas_scripts folder, into our working directory. 

In the copy of this file we have to modify the following lines:

* date=201607010000/to/201607150000  (where the data format is YYYYMMDDhhmm)

* nlat, slat, wlon, elon lines: these lines are the geographycal coordinates of our location of interest, in degrees. As of now, the software only allows for processing single grid point data, meaning that the *nlat* and *slat* values must be the same, and *wlon* and *elon* must be the same as well. 

As an example, for La Palma location, we would have:

.. code-block:: bash

    # Optional, use for spatial subset requests 90 to -90
    nlat=28.8
    # Optional, use for spatial subset requests 90 to -90
    slat=28.8
    # Optional, use for spatial subset requests -180 to 180
    wlon=-17.8
    # Optioanal, use for spatial subset requests -180 to 180
    elon=-17.8


Once this is set we can execute the main download_gdas script:

.. code-block:: bash

    >download_gdas.py -submit ds083.2_test

This will print a lot of text and will send us an email when the data is ready to be downloaded. Alternatively, we can execute the command 

.. code-block:: bash

    >download_gdas.py -get_status requestindex

in order to see if the data is ready or not. *requestindex* is a unique identifier for your query, that appears at the end of the block of text printed in the first command execution and also in the email that is received once the data is ready. The user needs to look for *Index* at the last lines of the text block that was printed in the previous command. Here there is an example:

.. code-block:: bash
    :emphasize-lines: 2
    
    1:0.2.3,3!7-0.2-1:0.3.0,3!7-0.2-1:0.3.5,3!7-0.2-1:0.6.1;level=76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,98,100;nlat=28.8;slat=28.8;wlon=-17.8;elon=-17.8;product=1\n\nds083.2: 
    Request Index 380581 added for <User> user.email@email.address\nSubset/Format-Conversion 
    Data Request 380581:\n\nYour Subset/Format-Conversion Data request has been submitted successfully.\nA summary of your request is given below.\n\nYour request will be processed soon. You will be informed via email\nwhen the data is ready to be picked up.\n\nYou may check request status of data requests you have submitted via the web link\nhttps://rda.ucar.edu/#ckrqst\n\nIf the information is CORRECT no further action is need.\nIf the information is NOT CORRECT, or if you have additional comments\nyou may email to rpconroy@ucar.edu (Riley Conroy} with corrections or comments.\n\nRequest Summary:\nIndex    : 380581\nID       : USER380581\nCategory : Subset/Format-Conversion Data\nStatus   : Queue\nDataset  : ds083.2\nTitle    : NCEP FNL Operational Model Global Tropospheric Analyses, continuing from July 1999\nUser     : user name\nEmail    : user.email@email.address\nDate     : 2019-07-19\nTime     : 07:13:11\nCompress : GZ\nRequest Detail:\nDate Limits          :  2019-07-15 00:00 to 2019-07-19 00:00\nParameter            :  HGT/PRES/TMP/R H/P WAT/A PCP/U GRD/V GRD/T CDC/LANDN/TOZNE\nLevel Type           :  ISBL:1000/975/950/925/900/850/800/750/700/650/600/550/500/450/400/350/300/250/200/150/100/50/20\nLatitude Limits      :  28.8 N to 28.8 S\nLongitude Limits     :  -17.8 W to -17.8 E\nProduct              :  Analysis\n\n\n\n\nFri Jul 19 13:13:11 UTC 2019\n\n\n

Once the data is ready to be retrieved, we can download it with the following command:

.. code-block:: bash

    download_gdas.py -download [RequestIndex]

Once the download has finished you will have your data ready to be processed.