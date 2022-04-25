#!/bin/sh
#SBATCH --job-name=RAMP_red_V0_500ts
#SBATCH --output=Run_bracket_big_%j.out
#SBATCH --error=Run_bracket_big_%j.err

#SBATCH --partition=workq
#SBATCH --ntasks-per-core 1
#SBATCH --ntasks-per-node 1

# in general use a node range
#SBATCH --nodes=1
#SBATCH --exclude=sn002,sn250,sn252
##SBATCH --nodelist=sn[081-090,218-223]

# is needed when allocating all cores
#SBATCH --exclusive
## the two sbatch options below are commented out!
#SBATCH --time=40:00:00
##SBATCH --mem-per-cpu=6000

#Load modules
module purge
module load hpcx/2.5.0/ompi
module load GCC/9.2.0-2.32

# Create ouput directory
DIR=RESULT_$SLURM_JOB_ID

# Restart functionaloty
RESTARTDIR=
ONLYLOADDESIGN=true
LOADDESIGN=0 #0: load x, 2: load xP
RESTARTFILE=${RESTARTDIR}/Restart01.dat
RESTARTITR=${RESTARTDIR}/Restart01_itr_f0.dat
RESTARTSOL=${RESTARTDIR}/RestartSol01.dat

# CREATE JOB FOLDER
mkdir -p ${DIR}

# COPY JOB SCRIPT FOR REFERENCE 
cp job.sh ${DIR}/job.sh

# Make 
make topopt -j

# Run program
srun ./topopt -nlvls 1   -xcmax 30 -ycmax 10 -zcmax 10  -nx 31 -ny 11  -nz 11  -maxItr 10   -viewer_binary_skip_info -workdir ${DIR}  -ksp_type fcg -mg_coarse_ksp_type preonly -mg_coarse_pc_type lu -ksp_initial_guess_nonzero true -restartFileVec $RESTARTFILE -restartFileItr $RESTARTITR -restartFileVecSol $RESTARTSOL -onlyLoadDesign $ONLYLOADDESIGN -loadDesign $LOADDESIGN -reduction true \
