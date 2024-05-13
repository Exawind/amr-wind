.. _turbine:

Turbine Simulation Walkthrough
==============================

Now that we have run our precursor simulation and saved inflow boundary condition data, we can run a turbine simulation. Here is the input file:

.. literalinclude:: ./turbines_inp.txt
   :linenos:

This file looks like the precursor input file, except the following changes:
* We're now reading in boundary condition data, not writing it out. Similarly, the x- and y- boundaries are no longer periodic, and we specify some extra characteristics about the inflow.
* `incflo.physics` now includes `Actuator`
* Similarly, we drop `ABLForcing` from `ICNS.source_terms`, and we replace it with three new forces: `BodyForce ABLMeanBoussinesq ActuatorForcing`. We then provide extra information about these forces. More on how to calculate these values down below.
* We add information about the turbines into the config file

When running an inflow-outflow simulation, you need to calculate information about `BodyForce` and `ABLMeanBoussinesq` from the precursor simulation. Do do this for this example, run the `calc_inflowoutflow_stats.py` as follow:

.. code-block:: console

    python calc_inflowoutflow_stats.py -sf /scratch/orybchuk/wakedynamics/amr-wind-tutorial/02_atmosphere/precursor/post_processing/abl_statistics14400.nc -ts 7200 -te 9000

This script gives you two important things:
* Information that you should copy-paste into the config file
* An `avg_theta.dat` file with important information for `ABLMeanBoussinesq`. It is critical that you copy or symbolically link this file into the same place as where you config file will run (usually /scratch/), otherwise AMR-Wind will fail with a mysterious error message

Because we're simulating turbines, we also need to include OpenFAST files for each of the turbines. In this demo, we use 2.8 MW turbines from NREL's open source [turbine repo](https://github.com/NREL/openfast-turbine-models/tree/master/IEA-scaled/NREL-2.8-127/OpenFAST). When simulating OpenFAST turbines through AMR-Wind instead of directly through OpenFAST, it is important to make the follow changes to the OpenFAST files:
* AeroDyn: Make sure `WakeMod` is 0
* ElastoDyn: Set the initial RPM `RotSpeed` and inital yaw angle `NacYaw` to reasonable values
* `*.fst`: Set `CompInflow` to be 2 and `OutFileFmt` to be 1
* ServoDyn: Make sure `DLL_FileName` points to a `libdiscon.so` file from ROSCO

When all is said and done, the /scratch/ directory where I run the turbine simulation looks like

.. code-block:: console

    amr_wind  avg_theta.dat	T0_AMRWind_AWAKEN/  T1_AMRWind_AWAKEN/	T2_AMRWind_AWAKEN/  turbine.i