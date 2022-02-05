read_data=shock
model=kilonova_uvboost
method=test
dist=40
delay=0
print_progress=True
max_time=360000
include_optical=False
include_uv=D1
sample=auto
export read_data model method dist delay print_progress max_time include_optical include_uv sample

## example of nested loop
#for read_data in 'shock' 'kilonova'; do
#	export read_data
#	for include_optical in 'r' 'False'; do
#		export include_optical
#		poetry run python ./produce-data/produce-data.py
#       	poetry run python ./parameter_estimation.py
#    done
#done

poetry run python parameter_estimation.py
#poetry run python ./produce-data/produce-data.py
