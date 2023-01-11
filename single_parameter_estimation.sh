read_data=kilonova
model=kilonova
dist=40
delay=24
print_progress=True
include_optical=False
include_uv=NUV_D
sample=auto
resume_previous=True
save_after_seconds=1200
dlogz_threshold=0.5
parallel=8
method=sample
redden=True
optical_delay=12 # no optical delay corresponds to default 12
export read_data model method dist delay print_progress include_optical include_uv sample resume_previous save_after_seconds parallel dlogz_threshold redden optical_delay

# produce data if necessary
poetry run python produce-data/produce-data.py

# run parameter estimation
poetry run python parameter_estimation.py

