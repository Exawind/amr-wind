.. _precursor:

Precursor (ABL) Walkthrough
===========================

Now that we have AMR-Wind compiled, we can start to run simulations. In this section, we'll run two simulations of a weakly convective atmospheric boundary layer. The first simulation will be of a transient "spinup" phase, and the second simulation will serve as a "precursor" for a turbine simulation.

Before we run anything, we first need to set up the [input](https://exawind.github.io/amr-wind/user/inputs.html) file (aka the "configuration"/"config" file). This is a text file, and its filename traditionally ends with `.i` or `.inp`. There are two general approaches to set one up: manually, or with the use of [amrwind-frontend](https://github.com/lawrenceccheung/amrwind-frontend). If I am setting up a big simulation, with many turbines and many refinement zones, I use amrwind-frontend. It shows the locations of all those objects, which is a great sanity check for expensive simulations. If I am running a simpler simulation, I will copy-paste text from a simulation that I have successfully run in the past. 

Spinup
------
Here is the content of our spinup config file:

.. literalinclude:: ./spinup_inp.txt
   :linenos:

To get an idea of what each of these lines means, see [here](https://exawind.github.io/amr-wind/user/inputs.html). Sometimes new input variables are added but aren't located in the docs, and in that case, you will need to look into the guts of the AMR-Wind code to deduce what those variables do.

We are going to run the spinup simulation for a duration of two hours. In this tutorial, we are simulation a fairly small domain with a coarse uniform grid resolution of 20 m. The spinup simulation does not employ any grid refinement.

I am running this simulation on NREL's supercomputer Eagle, and I have this file placed in the directory `/projects/wakedynamics/orybchuk/amr-wind-tutorial/02_atmosphere/spinup`. The contents of this directory are 

.. code-block:: console

    analyze_spinup.ipynb  logs/  runfile.sbatch  spinup.i

I like executing simulations in Eagle's scratch directory, so I also created `/scratch/orybchuk/wakedynamics/amr-wind-tutorial/02_atmosphere/spinup`. I `cd` into that `/scratch/` directory and then symbolically link to the config file in projects, with the command

.. code-block:: console

    ln -sf /projects/wakedynamics/orybchuk/amr-wind-tutorial/02_atmosphere/spinup/spinup.i .

Now, the contents of the scratch directory are 

.. code-block:: console

    spinup.i

To kick off the job on the Eagle's compute nodes, I use `runfile.sbatch`. Its contents are

.. literalinclude:: ./runfile.sbatch
   :linenos:

In this script, I symbolically link my AMR-Wind executable into the /scratch/ directory. I also run on two nodes, each of which have 36 cores, hence `-n 72`. Some of the other variables in this script, like `EXAWIND_DIR`, are specific to Eagle. You might not need to specify these in order to get AMR-Wind to run.

To kick off the job, go into the /projects/ directory, and execute

.. code-block:: console

    sbatch runfile.sbatch

You can view the job log by opening up `logs/job_output_filename.###.out`. Note: you need to manually create `logs/` before kicking off your simulation, otherwise your job will silently fail because of SLURM issues.

If your job is successful, the log file will look something like this at the end:

.. code-block:: console
    
    Writing plot file       plt14402 at time 7201
    Writing checkpoint file chk14402 at time 7201
    Time spent in InitData():    0.816465397
    Time spent in Evolve():      2404.252191
    Unused ParmParse Variables:
       [TOP]::incflo.delp(nvals = 3)  :: [0., 0., 0.]

    AMReX (22.05-29-g1305eb3d364d) finalized

And the /scratch/ directory will now look something like

.. code-block:: console

    amr_wind@   chk05400/	chk09000/  chk12600/  chk14402/   core.27800  core.9246   plt05400/  plt09000/  plt12600/  plt14402/	      spinup.i@
    chk00000/  chk01800/		    chk03600/		     chk07200/	chk10800/  chk14400/  core.27159  core.8529   plt00000/  plt01800/		  plt03600/		   plt07200/  plt10800/  plt14400/  post_processing/

Once the spinup simulation is done, it is important to sanity check that the fields make sense. One quick test to do this would be to look at the evolution of horizontally averaged vertical profiles. The Jupyter notebook `analyze_spinup.ipynb` contains code to help with this.

You can also open the `plt#####` files using Paraview or other software that AMReX is compatible with. These files show a volume of the instantaneous fields at that timestep.

Precursor simulation
--------------------
After sufficiently spinning up turbulence, I kick off a "precursor simulation". 

In the context of wind turbine LES, a precursor is a simulation that is run without a turbine for the explicit purpose of generating inflow boundary conditions. Due to the way LES works, spinup and precursor simulations are almost always run with cyclic boundary conditions. This means wind that exits the outflow simulation is then recirculated back into the inflow. If we want to simulate a statistically homogeneous atmosphere, that's fine. However, this is problematic if you have a wind turbine---turbines generate wakes, and we don't want wakes recirculating back into the inflow. So we run wind turbine simulations with a prescribed "inflow boundary condition" (where the wind data comes from the precursor) and an "outflow boundary condition" (usually a pressure BC).

Here is the content of the precursor simulation.

.. literalinclude:: ./precursor_inp.txt
   :linenos:

This file is almost identical to the spinup config file, except there are a few differences:
* I run this simulation with pre-defined regions where the mesh is refined. The refinement locations are specified in the `tagging` section, and I tell AMR-Wind to use `amr.max_level=2` levels of refinement. This means my finest grid cell is now 20 / 2 / 2 = 5 m wide.
* Because my finest mesh is now smaller, I also reduced my timestep by a factor of 4 to 0.125 seconds. 
* I am starting this simulation from the last timestep of the spinup simulation, using the `io.restart_file` line
* I am also now saving out boundary condition data, using the `ABL.bndry*` lines
* I plan to make some videos of hub-height planes and cross sections, so I have increased the frequency of `sampling.output_frequency`

Just like before, I create a projects directory `/projects/wakedynamics/orybchuk/amr-wind-tutorial/02_atmosphere/precursor` and a corresponding scratch directory `/scratch/orybchuk/wakedynamics/amr-wind-tutorial/02_atmosphere/precursor`  . After linking the config file, I kick off a job with `sbatch runfile.sbatch`, and I confirm that this job successfully ran to completion by checking the log file.

Just like before, you can visualize the volume files using Paraview. Here, we're going to visualize the `sampling.xy-domain` and `sampling.xz-domain` data as well. That data is saved to `/scratch/orybchuk/wakedynamics/amr-wind-tutorial/02_atmosphere/precursor/post_processing/sampling14400.nc`. You can analyze data directly out of that file, but I like to reformat the data so that it is spatially sorted. To reformat data, run the code in `reformat_precursor_planes.ipynb`. To visualize the reformatted data, run `viz_precursor_planes.ipynb`.