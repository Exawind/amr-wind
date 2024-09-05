#!/bin/bash
#SBATCH --nodes=8
#SBATCH --time=4:00:00        # Wall clock time (HH:MM:SS) - once the job exceeds this time, the job will be terminated (default is 5 minutes)
#SBATCH --account=fy210193        # WC ID
#SBATCH --job-name=GE2p8    # Name of job
#SBATCH --partition=batch,short # partition/queue name: short or batch
#SBATCH --qos=normal           # Quality of Service: long, large, priority or normal 
# Number of nodes - the number of nodes you have requested (for a list of SLURM environment variables see "man sbatch")
export nodes=$SLURM_JOB_NUM_NODES

module purge; module load cde/v3/gcc/10.3.0 cde/v3/openmpi/4.1.2-gcc-10.3.0 cde/v3/hdf5/1.10.6-gcc-10.3.0-openmpi-4.1.2 cde/v3/netcdf-c/4.8.1-gcc-10.3.0-openmpi-4.1.2 cde/v3/cmake/3.23.1

# Number MPI processes to run on each node (a.k.a. PPN)
# CTS1 has 36 cores per node and Ghost
#Chama has 16 cores per node
export cores=16
export ncpus=$((nodes * cores))
export OMP_PROC_BIND=spread 
export OMP_PLACES=threads

mpiexec --bind-to core --npernode $cores --n $ncpus /projects/wind_uq/lcheung/AMRWindBuilds/hfm.20240102/amr-wind/build/amr_wind Calibrate_dx2p5_1_EPS5.00_WS_9.0.inp
    
