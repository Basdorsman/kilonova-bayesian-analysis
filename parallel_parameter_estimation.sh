# If an argument will remain the same throughout all jobs it is a constant argument. If it changes, it is a parallel argument. Denote this as a space separated string.

# constant arguments
model='shock'
method='test' #test, timeout
dist=40
delay=0
print_progress='0' #'0'=False
max_time=360000
include_optical='False'
include_uv='D1'

# parallel arguments
# include_uv='include_uv=D1 include_uv=D1,D2' #"include_uv=D1 include_uv=D1,D2"
read_data='read_data=shock read_data=kilonova' #"read_data=shock read_data=kilonova"

export read_data model method dist delay print_progress max_time include_optical include_uv


# jobs
parallel poetry run python ./produce-data/produce-data.py ::: $read_data # produce necessary data (may not be necessary)
parallel --dry-run poetry run python parameter_estimation.py ::: $read_data # --dry-run
parallel --progress poetry run python parameter_estimation.py ::: $read_data



