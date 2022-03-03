# This shell script uses GNU parallel. https://www.gnu.org/software/parallel/parallel_tutorial.html


# ---- arguments ---- #
# If an argument will remain the same throughout all jobs it is a constant argument. If it changes, it is a parallel argument. In the latter case the variable must be formatted as: variable='variable=firstvalue variable=secondvalue'.

# constant arguments
method='sample' #test, sample
delay=0
print_progress=0 #'0'=False
sample='auto'
#include_optical='r'
#include_uv='NUV_D'
parallel=10
resume_previous='True'
save_after_seconds='False'
dlogz_threshold=0.5
dist=40


# parallel arguments
include_optical='include_optical=r include_optical=False'
include_uv='include_uv=NUV_D include_uv=NUV_D,D2' #"include_uv=D1 include_uv=D1,D2"
# read_data='read_data=shock read_data=kilonova read_data=kilonova_uvboost read_data=kilonova read_data=kilonova_uvboost' #"read_data=shock read_data=kilonova read_data=kilonova_uvboost"
read_data='read_data=kilonova read_data=kilonova_uvboost'
model='model=kilonova model=kilonova_uvboost'
# model='model=shock model=shock model=shock model=kilonova model=kilonova_uvboost' 
# model='model=kilonova model=kilonova_uvboost'
# dist='dist=40 dist=100 dist=160'
# dist='dist=100 dist=160'



export read_data model method dist delay print_progress include_optical include_uv sample resume_previous save_after_seconds parallel dlogz_threshold


# ---- jobs ---- #
# ::: means combine all possible combinations, :::+ means matched arguments.

# produce data
#parallel poetry run python ./produce-data/produce-data.py ::: $read_data ::: $dist # produce necessary data (perhaps not necessary)

# dry run (flag: --dry-run)
parallel --dry-run poetry run python parameter_estimation.py ::: $model :::+ $read_data ::: $include_optical :::+ $include_uv

# run (useful flag: --progress)
parallel --progress poetry run python parameter_estimation.py ::: $model :::+ $read_data ::: $include_optical :::+ $include_uv
