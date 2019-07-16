.. _other functionalities:

Other functionalities of the package
====================================

The molecularprofiles package comes with a set of functions that may be useful to interact with the meteorological data for the CTA.

In order to load them one needs to import these functions within a Python interpreter:

.. code-block:: python

    from molecularprofiles.utils.dataframe_ops import *

This will load all the functions contained within the dataframe_ops.py script.

For instance, once we have successfully loaded some data using the steps described in :doc:`the MolecularProfile class`. 

Once we have loaded the data we can start interacting with the dataframe that is created and, for instance, use a function to plot a wind-rose graphic of the direction and speed of the wind