#PBS -l select=1:ncpus=8:model=san 
#PBS -l walltime=03:00:00
#PBS -N log_par_est_kilonova
#PBS -j oe
#PBS -m bae
#PBS -M bas.dorsman@student.uva.nl
#PBS -k oed
#PBS -o /output_files/jobs/$PBS_JOBID.txt

cd $PBS_O_WORKDIR
poetry run python parameter_estimation.py
