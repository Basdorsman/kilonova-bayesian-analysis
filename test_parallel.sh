module load anaconda3/2021-05
#module load openmpi/3.1.6

source $HOME/.poetry/env

poetry run python test_parallel.py
