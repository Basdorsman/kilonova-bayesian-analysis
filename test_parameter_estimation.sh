read_data=shock
model=kilonova_uvboost
method=timeout
dist=40
delay=0
print_progress=True
max_time=360000
include_optical=False
include_uv=NUV_D
sample=auto
run_mode=external
export read_data model method dist delay print_progress max_time include_optical include_uv sample run_mode

poetry run python parameter_estimation.py
#poetry run python ./produce-data/produce-data.py
