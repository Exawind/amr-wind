.. _turbine:

Turbine Simulation Walkthrough
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
values are intended to be statistical values that are calculated from the precursor simulation. During the precursor, ABL Stats
output planar-averaged data at discrete times. To average this data in time and obtain values for the input arguments,
applying a script is the most direct approach. For this step, run the ``calc_inflowoutflow_stats.py`` script, provided in
the tools directory of the AMR-Wind repo, as follows:

.. code-block:: console

    python <path to amr-wind repo>/tools/calc_inflowoutflow_stats.py -sf <inflow-outflow job directory in scratch>/../precursor/post_processing/abl_statistics14400.nc -ts 7200 -te 9000

This script provides two important things:
* Information that you should copy-paste into the input file
* An `avg_theta.dat` file with important information for `ABLMeanBoussinesq`. It is critical that you copy or symbolically link this file into the same place as where you input file will run (usually /scratch/), otherwise AMR-Wind will fail with a mysterious error message

Because we're simulating turbines, we also need to include OpenFAST files for each of the turbines. In this demo, we use 2.8 MW turbines from NREL's open source [turbine repo](https://github.com/NREL/openfast-turbine-models/tree/master/IEA-scaled/NREL-2.8-127/OpenFAST). When simulating OpenFAST turbines through AMR-Wind instead of directly through OpenFAST, it is important to make the follow changes to the OpenFAST files:
* AeroDyn: Make sure `WakeMod` is 0
* ElastoDyn: Set the initial RPM `RotSpeed` and initial yaw angle `NacYaw` to reasonable values
* `*.fst`: Set `CompInflow` to be 2 and `OutFileFmt` to be 1
* ServoDyn: Make sure `DLL_FileName` points to a `libdiscon.so` file from ROSCO

When all is said and done, the /scratch/ directory where I run the turbine simulation looks like

.. code-block:: console

    amr_wind  avg_theta.dat	T0_AMRWind_AWAKEN/  T1_AMRWind_AWAKEN/	T2_AMRWind_AWAKEN/  turbine.i