read_data=shock
model=shock
dist=40
delay=0
print_progress=True
include_optical=False
include_uv=NUV_D
sample=auto
resume_previous=True
save_after_seconds=1200
dlogz_threshold=0.5
parallel=8
method=plot
export read_data model method dist delay print_progress include_optical include_uv sample resume_previous save_after_seconds parallel dlogz_threshold

#env=dorado-paramest-7nIrBfvE-py3.8
#home_env=$HOME/.cache/pypoetry/virtualenvs/$env
#path_to_envs=$HOME/.cache/pypoetry/virtualenvs
#cp -r $home_env $path_to_envs
#echo "cp -r $home_env $path_to_envs"
#poetry config virtualenvs.path $path_to_envs
#poetry env list

#poetry install
#poetry env info
#poetry run python ./produce-data/produce-data.py 
poetry run python parameter_estimation.py

