# constant arguments
model='shock'
method='test' 
dist=40
delay=0
print_progress='0' #False
max_time=360000
include_optical='r'
export model method dist delay print_progress max_time include_optical

# parallel arguments
include_uv="D1 D1,D2"
read_data="shock kilonova"

# jobs
parallel poetry run python ./produce-data/produce-data.py ::: $read_data # produce necessary data (may not be necessary)
parallel --dry-run poetry run python parameter_estimation.py ::: $include_uv ::: $read_data # --dry-run
parallel --progress poetry run python parameter_estimation.py ::: $include_uv ::: $read_data



