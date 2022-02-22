#!/bin/sh
##SBATCH -N 1
#SBATCH --ntasks=2
#SBATCH -t 6-00:00:00                                
#SBATCH -J lightcurve
#SBATCH --partition=neutron-star
#SBATCH --mem 252000
#SBATCH --nodelist=helios-cn[039]


# ####SBATCH --mem 252000
# #SBATCH --tasks-per-node=2                                                                                                                              


module purge
module load anaconda3/2021-05
module load openmpi/3.1.6

#export PATH="$HOME/.poetry/bin:$PATH"
source $HOME/.poetry/env

echo $PATH
#srun -n $SLURM_NTASKS --mpi=pmix_v2 
sh test_parameter_estimation.sh 2

#poetry run python lightcurve.py
