.. _run:

Running AMR-Wind
=================

To run :program:`amr_wind`, the user must provide a text file containing inputs
describing the problem and any additional command-line arguments that override
the parameters in the input file for that particular invocation of the
executable.

.. code-block:: console

   # Parse input parameters from `inputs.abl` but change max_step to 20
   $ ./amr_wind inputs.abl time.max_step=20


See :ref:`input file description <inputs>` for a list of AMR-Wind input flags and proper syntax.

Restarting AMR-Wind
--------------------

To restart a case, an initial simulation must be run with AMR-Wind checkpoints. These checkpoints 
are a slice of the solution, and solution settings, that will provide initial input for a restart. Check the
:ref:`time inputs <inputs_time>` and :ref:`i/o inputs <inputs_io>` for input flags to set up checkpoints.

#. Change your ``time.start_time`` to the restart time (optional)
#. Set ``io.restart_file`` to your AMR-Wind checkpoint directory
#. Re-submit case as above

Restarting AMR-Wind with OpenFAST Turbines
-------------------------------------------

.. note::
   Currently, AMR-Wind automatically creates an OpenFAST .chkp file for every AMR-Wind checkpoint. These .chkp files 
   are found within each OpenFAST turbine directory.

#. Change your ``time.start_time`` to the restart time (optional)
#. Set ``io.restart_file`` to your AMR-Wind checkpoint directory
#. Set ``Actuator.T1.openfast_sim_mode = restart`` for each turbine (i.e. T1 in this example)
#. Set ``Actuator.T1.openfast_restart_file`` as the OpenFAST checkpoint file. Note that this must be a relative path from your amr-wind case root, and the ".chkp" must be removed from the filename
#. Set ``Actuator.T1.openfast_start_time`` to the restart time and double check that ``Actuator.T1.openfast_stop_time`` is ok
#. Re-submit case as above