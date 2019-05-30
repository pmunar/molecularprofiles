.. _setup:

Setup
=====

In order to have molecularprofiles ready to run, the steps you need to do are:

* Install or check for a Python installation in your system
* Install or check for installation of needed Python libraries
* Install molecularprofiles


Install or check for a Python installation in your system
---------------------------------------------------------

If you are using Linux, it is probable that you already have a Python installation ready in your system. However, I recommend you to make a parallel installation using `Anaconda <https://www.anaconda.com/>`__ distribution. This is an open-source Python distribution which allows the user to easily install packages in a clean and transparent way.

I will assume you have this distribution installed from now on. For other ways to install Python and associated packages, there is a lot of information around in the net.

Install or check for installation of needed Python libraries
------------------------------------------------------------

The additional needed Python libraries that molecularprofiles needs to run properly are:

* os
* sys
* gc
* configparser
* datetime
* contextlib
* argparse
* multiprocessing
* tqdm
* pandas
* scipy
* matplotlib

You can install them either by using `pip <https://pypi.org/project/pip/>`__ or `conda <https://docs.conda.io/en/latest/>`__, depending on your Python installation. In order to install the missing library (sys, in the next example) you can install it by typing:

.. code-block:: bash

    pip install sys 

or

.. code-block:: bash

    conda install sys 

Once you have all the needed Python libraries you can proceed to install molecularprofiles.

Install molecularprofiles
-------------------------

You can get the molecularprofiles package from the GitHub repository. In order to download it, just type:

.. code-block:: bash

    > git clone https://github.com/pmunar/molecularprofiles

Once downloaded, put the molecularprofiles folder in the directory where you want to 
install it (if you downloaded it elsewhere):

.. code-block:: bash

    > mv molecularprofiles /path/where/you/want/to/install/it

With the package in place, cd to the molecularprofiles directory:

.. code-block:: bash

    > cd molecularprofiles

and install it with the pip command:

.. code-block:: bash

    > pip install .

Note the "." at the end of the order. It is important!

Once installed (it takes a few seconds) it is almost ready to run.

Before running, the init-molecularprofiles.sh script must be executed. It sets some
usefull and important environment variables. But before running this script
there is one environment variable that needs to be set:

1- export the molecularprofiles_DIR variable. You can do it from the terminal:

.. code-block:: bash

    > export MOLECULARPROFILES_DIR=/example/path/molecularprofiles

In order to make the process more confortable for you, we recommend you to put this export within your .bashrc file.

Once it is done, the script can be executed:

.. code-block:: bash

    > ./molecularprofiles-init.sh

This script must be executed every time that the molecularprofiles package wants to be
used. An easy solution is to make an alias and put it into the .bashrc file.
An example of the line that would go into the .bashrc file:

.. code-block:: bash

    alias init-molecularprofiles=". /path/where/you/installed/it/molecularprofiles/init-molecularprofiles.sh"

After that, before using the software, type:

.. code-block:: bash

    > init-molecularprofiles

from wherever directory and the package will be ready.