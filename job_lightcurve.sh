#PBS -q devel
#PBS -l select=1:model=san
#PBS -l walltime=00:01:00
#PBS -N lightcurve
#PBS -j oe

cd $PBS_O_WORKDIR
poetry run python lightcurve.py
