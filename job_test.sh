#!/bin/bash
#SBATCH --job-name=testjob
#SBATCH --output=/home/bdorsma/my_outputs/test_lightcurve_plot_%j.out
#SBATCH --error=/home/bdorsma/my_outputs/test_lightcurve_plot_%j.err
#SBATCH -N 1
#SBATCH --cpus-per-task=2  
#SBATCH -t 6:00:00
##SBATCH --partition=neutron-star
##SBATCH --nodelist=helios-cn[039]

# Create work folder
#mkdir -p /hddstore/$USER
#export work_folder=$(mktemp -d -p /hddstore/$USER)
export work_folder=/hddstore/$USER/workfolder
mkdir -p $work_folder

echo $work_folder $SLURMD_NODENAME

# Copy files over to work folder
project=kilonova_bayesian_analysis
cp -r /home/$USER/$project $work_folder

# Go to output folder
cd $work_folder/$project

# Run job
module purge
module load anaconda3/2021-05
#module load openmpi/3.1.6

#curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | POETRY_HOME=/hddstore/bdorsma/poetry python -
#source /hddstore/bdorsma/poetry/env
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python -
source $HOME/.poetry/env


#srun -n $SLURM_NTASKS --mpi=pmix_v2 
sh test_parameter_estimation.sh


# Copy results
export results_folder=/home/$USER/my_outputs/$SLURM_JOB_ID
mkdir -p $results_folder
cp -r $work_folder/$project/output_files/plots $results_folder

# remove working folder
#rm -rf $work_folder
