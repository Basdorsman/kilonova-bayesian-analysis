read_data=shock
model=kilonova_uvboost
method=sample
dist=40
delay=0
print_progress=True
include_optical=False
include_uv=NUV_D
sample=auto
resume_previous=False
save_after_seconds=1000
parallel=True
dlogz_threshold=10000
export read_data model method dist delay print_progress include_optical include_uv sample resume_previous save_after_seconds parallel dlogz_threshold

poetry run python parameter_estimation.py
#poetry run python ./produce-data/produce-data.py
