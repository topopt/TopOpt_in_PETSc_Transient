#!/bin/sh
#SBATCH --job-name=tot_red
#SBATCH --output=Run_%j.out
#SBATCH --error=Run_%j.err

#SBATCH --partition=workq
#SBATCH --ntasks-per-core 1
#SBATCH --ntasks-per-node 32

# in general use a node range
#SBATCH --nodes=16

# is needed when allocating all cores
#SBATCH --exclusive
## the two sbatch options below are commented out!
#SBATCH --time=40:00:00

module purge
module load hpcx/2.5.0/ompi
module load GCC/9.2.0-2.32

# Create ouput directory
DIR=RESULT_Fig5c_$SLURM_JOB_ID
RESTARTDIR=

RESTARTFILE=${RESTARTDIR}/Restart00.dat
RESTARTITR=${RESTARTDIR}/Restart00_itr_f0.dat
RESTARTSOL=${RESTARTDIR}/RestartSol00.dat

# CREATE JOB FOLDER
mkdir -p ${DIR}

# COPY JOB SCRIPT FOR REFERENCE 
cp Job_fig5c.sh ${DIR}/Job_fig5c.sh

#Static re-evaaluation in transient setting
srun ./topopt -nlvls 4 -volfrac 0.3 -V0 1.0 -xcmax 40 -ycmax 20 -zcmax 20  -nx 193 -ny 97  -nz 97  -maxItr 600 -ncp 1 -tsc 50  -dalpha 0.01 -dbeta 0.01 -filter 2 -Rmin 2.4 -t1 50  -viewer_binary_skip_info -workdir ${DIR} -penalK 5 -sInt 2  -projectionFilter true -eta 0.15  -ti 1 -movlim 0.2 -ksp_initial_guess_nonzero true -Emax 1e3 -Rmax 4.16666666667e-4 -objective 3  -optimizationRegion 2 -reduction true -omega0 1.3 -boundaries 1 -loadloc 1 -loadsignal 1 \

exit
