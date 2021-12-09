#!/bin/bash

#SBATCH --job-name=amr_abl
#SBATCH --account=hfm
#SBATCH --nodes=29
##SBATCH --nodes=15
#SBATCH --time=12:15:00
#SBATCH --output=out.%x_%j
##SBATCH --qos=high

source ~/.bash_profile
amr_env

ranks_per_node=36
mpi_ranks=$(expr $SLURM_JOB_NUM_NODES \* $ranks_per_node)
export OMP_NUM_THREADS=1  # Max hardware threads = 4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#amr_exec=/home/lmartine/amr-wind-tony/build/amr_wind
#amr_exec=/home/lmartine/amr-wind-marc/build/amr_wind
amr_exec=/home/lmartine/amr-wind/build/amr_wind

echo "Job name       = $SLURM_JOB_NAME"
echo "Num. nodes     = $SLURM_JOB_NUM_NODES"
echo "Num. MPI Ranks = $mpi_ranks"
echo "Num. threads   = $OMP_NUM_THREADS"
echo "Working dir    = $PWD"

cp ${amr_exec} $(pwd)/amr_wind

# Case with epsilon/c=4
mkdir ec4
cd ec4
cp ${amr_exec} $(pwd)/amr_wind
cp ../fixed_wing.i ../static_box.txt ../NACA64_A17.dat .

srun -n 64  -c 1 --cpu_bind=cores $(pwd)/amr_wind fixed_wing.i Actuator.FixedWingLine.epsilon_chord=4.0 4.0 4.0  amr.max_level=1 > out.log 2>&1
cd ..

# Case with epsilon/c=2
mkdir ec2
cd ec2
cp ${amr_exec} $(pwd)/amr_wind
cp ../fixed_wing.i ../static_box.txt ../NACA64_A17.dat .

srun -n 64  -c 1 --cpu_bind=cores $(pwd)/amr_wind fixed_wing.i Actuator.FixedWingLine.epsilon_chord=2.0 2.0 2.0 amr.max_level=1 > out.log 2>&1
cd ..

# Case with epsilon/c=1
mkdir ec1
cd ec1
cp ${amr_exec} $(pwd)/amr_wind
cp ../fixed_wing.i ../static_box.txt ../NACA64_A17.dat .

srun -n 64  -c 1 --cpu_bind=cores $(pwd)/amr_wind fixed_wing.i Actuator.FixedWingLine.epsilon_chord=1.0 1.0 1.0 amr.max_level=2 > out.log 2>&1
cd ..

# Case with epsilon/c=0.5
mkdir ec05
cd ec05
cp ${amr_exec} $(pwd)/amr_wind
cp ../fixed_wing.i ../static_box.txt ../NACA64_A17.dat .

srun -n 64  -c 1 --cpu_bind=cores $(pwd)/amr_wind fixed_wing.i Actuator.FixedWingLine.epsilon_chord=0.5 0.5 0.5 amr.max_level=2 > out.log 2>&1
cd ..

# Case with epsilon/c=0.25
mkdir ec025
cd ec025
cp ${amr_exec} $(pwd)/amr_wind
cp ../fixed_wing.i ../static_box.txt ../NACA64_A17.dat .

srun -n 64  -c 1 --cpu_bind=cores $(pwd)/amr_wind fixed_wing.i Actuator.FixedWingLine.epsilon_chord=0.25 0.25 0.25 amr.max_level=3 > out.log 2>&1
cd ..










