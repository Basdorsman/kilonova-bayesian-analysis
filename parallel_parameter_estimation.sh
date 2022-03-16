# This shell script uses GNU parallel. https://www.gnu.org/software/parallel/parallel_tutorial.html


# ---- arguments ---- #
# If an argument will remain the same throughout all jobs it is a constant argument. If it changes, it is a parallel argument. In the latter case the variable must be formatted as: variable='variable=firstvalue variable=secondvalue'.

# constant arguments
method='sample' # test, sample
# delay=0
print_progress=0 #'0'=False
sample='auto'
include_optical='False' # include_optical='r'
include_uv='NUV_D'
parallel=8
resume_previous='True'
save_after_seconds=1200
dlogz_threshold=0.5
#dist=40
#read_data='shock'
# model='kilonova'

# parallel arguments
# include_optical='include_optical=r include_optical=False include_optical=r'
# include_uv='include_uv=NUV_D include_uv=NUV_D,D2 include_uv=False' #"include_uv=D1 include_uv=D1,D2"
 #"read_data=shock read_data=kilonova read_data=kilonova_uvboost"
# read_data='read_data=kilonova read_data=shock'
# model='model=kilonova model=kilonova_uvboost'
# dist='dist=40 dist=100 dist=160'
# dist='dist=100 dist=160'
delay='delay=2 delay=4 delay=8'
dist='dist=40 dist=100 dist=160'

# for a full heatmap use matched arguments read data and model:
read_data='read_data=shock read_data=kilonova read_data=kilonova_uvboost read_data=shock read_data=kilonova read_data=shock read_data=kilonova_uvboost'
# model='model=kilonova model=shock'
model='model=shock model=shock model=shock model=kilonova model=kilonova model=kilonova_uvboost model=kilonova_uvboost' 


export read_data model method dist delay print_progress include_optical include_uv sample resume_previous save_after_seconds parallel dlogz_threshold


# ---- jobs ---- #
# ::: means combine all possible combinations, :::+ means matched arguments.

# produce data
parallel poetry run python ./produce-data/produce-data.py ::: $delay ::: $dist ::: $read_data :::+ $model # produce necessary data (perhaps not necessary)

# dry run (flag: --dry-run)
parallel --dry-run poetry run python parameter_estimation.py ::: $delay ::: $dist ::: $read_data :::+ $model

# run (useful flag: --progress)
parallel --progress poetry run python parameter_estimation.py ::: $delay ::: $dist ::: $read_data :::+ $model
