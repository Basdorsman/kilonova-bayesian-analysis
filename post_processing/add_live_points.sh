read_data=shock
model=kilonova_uvboost
dist=40
delay=0
print_progress=True
include_optical=False
include_uv=NUV_D
sample=auto
resume_previous=True
save_after_seconds=1200
dlogz_threshold=1000
parallel=8
method=sample
export read_data model method dist delay print_progress include_optical include_uv sample resume_previous save_after_seconds parallel dlogz_threshold

cd ../
poetry run python parameter_estimation.py
