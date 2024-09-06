.. _turbine:

Turbine simulation walkthrough
==============================

Now that we have run our precursor simulation and saved inflow boundary condition data, we can run a turbine simulation. Here is the input file:

.. literalinclude:: ./turbines_inp.txt
   :linenos:

This file looks like the precursor input file, except the following changes:

* We are now reading in boundary condition data, not writing it out. Similarly, the x- and y- boundaries are no longer periodic, requiring 
  us to specify some extra characteristics about the inflow.
* ``incflo.physics`` now includes ``Actuator`` so that turbine-related calculations can take place
* Similarly, we remove ``ABLForcing`` from ``ICNS.source_terms``, and we replace it with three new forces: 
  ``BodyForce ABLMeanBoussinesq ActuatorForcing``. We then provide extra information about these forces.
  More on how to calculate these values down below.
* We also add details under ``ABL`` to inform the wall function.
* We add information about the turbines and surrounding mesh refinements into the input file. Because of the finer mesh, the small dt is
  required for this simulation.

.. collapse:: Further details about these source terms and wf specifications

    Now that this simulation is not periodic and will feature turbine wakes, the previous approach to forcing the velocity field
    (ABLForcing) is not viable. In addition to the inflow information from the precursor simulation, we also want to mimic
    the forcing terms. This is the purpose of BodyForce.

    ABLMeanBoussinesq modifies the BoussinesqBuoyancy term to make it reference a profile of average temperature values instead of
    a single reference temperature. Instead of modifying the physics, this term is intended to modify the pressure distribution
    to make it more compatible with the outflow condition.

    ActuatorForcing indicates to the momentum equation that terms from the actuator turbine model should be included. Both the
    Actuator physics instance and the source_term specification need to be present for the simulation to run as intended.

    Similar to how ABLForcing is not a viable option in an inflow-outflow simulation with turbines,
    wall function implementations that rely on planar averages are not viable either. Instead, time averages of planar-averaged
    precursor data from is passed to the inflow-outflow simulation to inform that wall model.

When running an inflow-outflow simulation, information must be provided for the new source terms and the wall function. These
values are intended to be statistical data calculated from the precursor simulation. During the precursor, ABL Stats
output planar-averaged data at discrete times. To average this data in time and obtain values for the input arguments,
applying a script is the most direct approach. For this step, run the ``calc_inflowoutflow_stats.py`` script, provided in
the tools directory of the AMR-Wind repo, as follows:

.. code-block:: console

    python <path to amr-wind repo>/tools/calc_inflowoutflow_stats.py -sf <path to precursor directory>/post_processing/abl_statistics14400.nc -ts 7200 -te 9000


The ``-ts`` and ``-te`` arguments communicate the start time and end time for the time averaging, respectively. Running the
script provides two important things:

* Complete input lines to be copy-pasted into the turbine input file, which feature time-averaged information
* An ``avg_theta.dat`` file with an average temperature profile for ``ABLMeanBoussinesq``. It is critical that this file 
  is available where the job will run; otherwise, AMR-Wind will fail by not being able to find the file.

Although the script provides a value for the ``BodyForce.magnitude``, this has been commented out after being pasted into the
input file. This is because we prefer to use a timetable for the body force, which was output by the precursor simulation.
If the ``forcing_timetable_output_file`` argument is not used in a precursor simulation setup, this timetable file will not be available.
In that case, using a constant, uniform body force vector is the standard approach, and this script provides the
time-averaged value from the precursor to populate ``BodyForce.magnitude``.

Because we are simulating turbines, OpenFAST files need to be included for each of the turbines. In this walkthrough, we use
2.8 MW turbines from NREL's `open source turbine repository <https://github.com/NREL/openfast-turbine-models/tree/master/IEA-scaled/NREL-2.8-127/OpenFAST>`_.
When simulating OpenFAST turbines through AMR-Wind instead of directly through OpenFAST, it is important to make the
following changes to the OpenFAST files:

* AeroDyn: Make sure ``WakeMod`` is ``0``
* ElastoDyn: Set the initial RPM ``RotSpeed`` and initial yaw angle ``NacYaw`` to reasonable values
* ``*.fst``: Set ``CompInflow`` to be ``2`` and ``OutFileFmt`` to be ``1``
* ServoDyn: Make sure ``DLL_FileName`` points to a ``libdiscon.so`` file from ROSCO. If you compiled using exawind-manager,
  see the :ref:`section of the documentation <rosco-dyn-lib>` that discusses how to determine the correct path to this ROSCO library file.

When all of these steps are complete, the job directory for running the turbine simulation includes

.. code-block:: console

    avg_theta.dat  T0_AMRWind_AWAKEN/  T1_AMRWind_AWAKEN/  T2_AMRWind_AWAKEN/  turbine.inp

and the turbine inflow-outflow simulation is ready to be submitted.