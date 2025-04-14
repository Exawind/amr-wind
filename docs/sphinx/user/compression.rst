.. _compression:

Compressing field output files
==============================

When writing out field output files, sometimes called plot files or plt files, it is possible to output the files with HDF5 and compress these with ZFP. To enable this capability, AMR-Wind must be compiled with ``AMR_WIND_ENABLE_HDF5`` and ``AMR_WIND_ENABLE_HDF5_ZFP``. The input file :ref:`documentation <inputs_io>` documents the input options to output these files. The `AMReX HDF5 documentation <https://amrex-codes.github.io/amrex/docs_html/IO.html#hdf5-plotfile-compression>`_ further documents valid options.

While generating the compressed files is straightforward, it takes a little work to read these in various post-processing utilities.

ParaView
--------

To read ZFP compressed HDF5 files into ParaView, ParaView must have the environment variable ``HDF5_PLUGIN_PATH`` populated and pointing to the HDF5 ZFP plugin location. Within a spack environment with the ``h5z-zfp`` package installed:

.. code-block:: console

    $ export HDF5_PLUGIN_PATH=$(spack location -i h5z-zfp)/plugin
    $ paraview

Then the user can load the file as usual.

yt
--

To read HDF5 files into `yt <https://yt-project.org>`_, the HDF5 file attributes need to be modified. We have provided a script to do so in `tools` directory of this repository. It requires the ``h5py`` and ``hdf5plugin`` Python packages.

.. code-block:: console

    $ python modify_hdf5_attributes.py -f plt00000 plt00010

The resulting file can then be opened in yt.
