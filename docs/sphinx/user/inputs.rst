.. _inputs:

AMR-Wind inputs file
=====================

To run :program:`amr_wind`, the user must provide a text file containing inputs
describing the problem and any additional command-line arguments that override
the parameters in the input file for that particular invocation of the
executable.

.. code-block:: console

   # Parse input parameters from `inputs.abl` but change max_step to 20
   $ ./amr_wind inputs.abl time.max_step=20


The input file is a simple text file containing ``key = value`` entries for the
input parameters. The text file can include comments, any text beginning with
``#`` till the end of line (EOF) is interpreted as comments and ignored by the
parser. Input file processing is handled by `AMReX ParmParse library
<https://amrex-codes.github.io/amrex/docs_html/Basics.html#parmparse>`_. This
section documents the various input file parameters and their default values (if
available). In :program:`amr_wind`, the input file is broken into *sections*
indicated by a namespace prefix. For example, all inputs related to the problem
domain are prefixed with ``geometry.`` and so on. A sample input file is shown below

.. literalinclude:: ./amr_wind_inputs.txt
   :linenos:


Input file reference
---------------------

The AMR-Wind input file is organized in the following sections

======================= ============================================================
Section                 Description
======================= ============================================================
``geometry``            Computational domain information
``amr``                 Mesh refinement controls
``time``                Simulation time controls
``io``                  Input/Output controls
``incflo``              CFD algorithm and physics controls
``transport``           Transport equation controls
``turbulence``          Turbulence model controls 
``ABL``                 Atmospheric boundary layer (ABL) controls
``SyntheticTurbulence`` Inject turbulence using body forces
``Momentum sources``    Activate Momentum source terms and their parameters
``Boundary conditions`` Boundary condition types and gradients
``MLMG options``        Multi-Level Multi-Grid Linear solver options
``Tagging``             Static and dynamic refinement options
``Sampling``            Data probes to sample field data during simulations
``Averaging``           Time averaging and correlations
======================= ============================================================

This section documents the parameters available within each section.

.. note::

   - Boolean flags (``true/false``) can also be indicated using integers in the text file
     and uses the convention ``0 = False`` and ``1 = True``.

   - Quotes around strings are optional.

   - If an input is repeated only the last one is used.


.. toctree::
   :maxdepth: 2

   inputs_geometry.rst
   inputs_amr.rst
   inputs_time.rst
   inputs_io.rst
   inputs_incflo.rst
   inputs_transport.rst
   inputs_turbulence.rst
   inputs_Momentum_Sources.rst
   inputs_ABL.rst
   inputs_SyntheticTurbulence.rst
   inputs_Static_Refinement.rst
   inputs_Boundary_conditions.rst
   inputs_MLMG.rst
   inputs_AMR_Tagging.rst
   inputs_Sampling.rst
   inputs_Averaging.rst
   inputs_KineticEnergy.rst
   inputs_Enstrophy.rst
   inputs_Actuator.rst
