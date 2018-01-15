This is the README file for the ecmwf parf of these set of functions to check and plot the molecular profiles.

In this folder there are 3 python files:

- download_ecmwf.py
- molecularprofiles_class.py
- readecmwf.py

In the following lines those files are described in detail.

***********************
* - download_ecmwf.py *
***********************

This python script allows the user to download ecmwf grib2 format data for a selected location on the Earth by giving a latutide and longitude value. 

**Before starting**, and in order to be able to download the data, the user must have registered on the ecmwf site 

https://apps.ecmwf.int/registration/ 

Once registered, you  must login in the web: 

https://apps.ecmwf.int/auth/login/

Retrieve you key at 

https://api.ecmwf.int/v1/key/

Note that the key expires in 1 year. You will receive an email to the registered email address 1 month before the expiration date with the renewal instructions.

Copy the information in this page and paste it in the file 
	$HOME/.ecmwfapirc (Unix/Linux) 
or 
	%USERPROFILE%\.ecmwfapirc (Windows: usually in C:\Users\<USERNAME>\.ecmwfapirc ; see how to create a file with a leading dot) 

This is the content of $HOME/.ecmwfapirc (Unix/Linux) or %USERPROFILE%\.ecmwfapirc (Windows):

{
    "url"   : "https://api.ecmwf.int/v1",
    "key"   : "XXXXXXXXXXXXXXXXXXXXXX",
    "email" : "john.smith@example.com"
}

Once done, we continue explaining the functions.

This python script that allocates 3 python functions. To run them, you need to have installed numpy, calendar and ecmwfapi Python packages. The three functions are:

	-  retrieve_interim:
	
	A function to demonstrate how to iterate efficiently over several years and months etc., for a particular interim_request.
        Change the variables below to adapt the iteration to your needs.
        It calls the function "request_ecmwf"

       	  - inputs:
            *date_start : (string) starting date of the period to download in "YYYY-MM-DD" format
            *date_end   : (string) ending date of the period to download in "YYYY-MM-DD" format
            *latitude   : (float) geographical latitude of the place of interest in degrees
            *longitude  : (float) geographical longitude of the place of interest in degrees. If negative it will be transformed to positive.
            *outtag (optional) : (string) tag to append to the output file name

          - outputs:
            *outfile : grib2 format files, as many as months from date_start until date_end (e.g. if date_start = "2017-01-01" and date_end = "2017-03-05" it will produce three files, one containing all January data, one containing all February data and another one containing the first 5 days of March data)

	- request_ecwmf:
	
	This function downloads ERA-Interim forecasts, on pressure levels. The data volume for all pressure level data is about 5GB per day, and all pressure level data for more than a single day will exceed the WebAPI limit of 600.000 fields. Thus please restrict the download to what you really need.
   	It is likely you need to split your request, this is best done by time periods, ie. first download for month 1, then for month 2, etc.
	It calls the function "find_nearest".

	  - inputs:
            *date_i : string, initial date in string format as "YYYY-MM-DD"
            *date_f : string, final date in string format as "YYYY-MM-DD"
            *lat : float, latidude of interest in degrees
            *lon : float, longitude of interest in degrees
            *outfile : string, output file name

	  - output:
   	    *outfile : grib2 format file containing the downloaded data

	- find_nearest:
	This function is to find the nearest grid position to a given latitude or longitude. It compares the desired value to an array of values and find the closest one (in value).
	  - input:
	    *a: array-like. Values to compare
	    *num: float. Desired value to compare to array "a".
	  - output:
	    *value in "a" that is closest to the value "num".


The script can be called outside the python interpreter.
Example 1:

python download_ecmwf.py 2015-01-01 2016-02-07 28.5 -17.8

Example 2:

python download_ecmwf.py 2015-01-01 2016-02-07 28.5 -17.8 grib_files_downloaded

The second example will create the files with the tag name "grib_files_downloaded" on their file name.


***********************
* - download_ecmwf.py *
***********************

