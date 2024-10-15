Walkthrough
===========

This section demonstrates a typical AMR-Wind workflow, walking
through the steps required to simulate wind turbines in a turbulent atmospheric
boundary flow. The compilation instructions outline how to use Exawind-Manager,
a Spack-based package management tool customized for the ExaWind software suite,
of which AMR-Wind is an integral part. The turbulent flow conditions are established
through precursor simulations, and then turbines are placed in the flow. 

.. note::

   This tutorial is intended to provide an example of how AMR-Wind is often used, but there
   are many variations and alternative workflows that AMR-Wind provides. Please 
   consult the :doc:`../user/user`, especially the :ref:`capabilities list<capabilities>`
   and :ref:`input file reference <input-file-ref>`, for additional details on other
   AMR-Wind features and options.

.. toctree::
   :glob:
   :maxdepth: 2

   compiling.rst
   precursor.rst
   turbine.rst
   terrain.rst