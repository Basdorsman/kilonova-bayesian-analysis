model='shock'
method='timeout' 
dist=40
delay=0
print_progress=True
max_time=360000
include_uv='D1'
export model method dist delay print_progress max_time include_uv

# example of nested loop
for read_data in 'shock' 'kilonova'; do
	export read_data
	for include_optical in 'r' 'False'; do
		export include_optical
		poetry run python ./produce-data/produce-data.py
        	poetry run python ./parameter_estimation.py
    done
done
