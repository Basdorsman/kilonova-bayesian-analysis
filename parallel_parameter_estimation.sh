# This shell script uses GNU parallel. https://www.gnu.org/software/parallel/parallel_tutorial.html


# ---- arguments ---- #
# If an argument will remain the same throughout all jobs it is a constant argument. If it changes, it is a parallel argument. In the latter case the variable must be formatted as: variable='variable=firstvalue variable=secondvalue'.

# constant arguments
method='timeout' #test, timeout
delay=0
print_progress='0' #'0'=False
max_time=360000
include_optical='False'
include_uv='D1'

# parallel arguments
# include_uv='include_uv=D1 include_uv=D1,D2' #"include_uv=D1 include_uv=D1,D2"
read_data='read_data=shock read_data=kilonova read_data=kilonova_uvboost read_data=shock read_data=kilonova read_data=shock read_data=kilonova_uvboost' #"read_data=shock read_data=kilonova read_data=kilonova_uvboost"
model='model=shock model=shock model=shock model=kilonova model=kilonova model=kilonova_uvboost model=kilonova_uvboost' 
dist='dist=40 dist=100 dist=160'

export read_data model method dist delay print_progress max_time include_optical include_uv


# ---- jobs ---- #

# produce data
#parallel poetry run python ./produce-data/produce-data.py ::: $read_data ::: $dist # produce necessary data (perhaps not necessary)

# dry run (flag: --dry-run)
parallel --dry-run poetry run python parameter_estimation.py ::: $read_data :::+ $model ::: $dist

# run (useful flag: --progress)
parallel --progress poetry run python parameter_estimation.py ::: $read_data :::+ $model ::: $dist
