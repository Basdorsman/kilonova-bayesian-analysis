# This shell script uses GNU parallel. https://www.gnu.org/software/parallel/parallel_tutorial.html

module load gnuparallel/20220322


# ---- arguments ---- #

# If an argument will remain the same throughout all jobs it is a constant argument. If it changes, it is a parallel argument. In the latter case the variable must be formatted as: variable='variable=firstvalue variable=secondvalue'.

# ---- constant arguments ---- #
method=sample # plot, sample
print_progress=0 #'0'=False
sample=auto
parallel=8
resume_previous=True
save_after_seconds=1800
dlogz_threshold=0.5
redden=True
optical_delay=12 # No optidal delay corresponds to default 12
#dist=70
#read_data=kilonova_uvboost
#model=kilonova_uvboost
#include_optical=r
#include_uv=False
#delay=0


# ---- parallel arguments ---- #

# ---- distances for figure 2 and 6 ---- #
dist='dist=40 dist=70 dist=100 dist=130 dist=160'

# ---- read data and model for figure 2 and 6 ---- #
read_data='read_data=shock read_data=kilonova read_data=kilonova_uvboost read_data=shock read_data=kilonova read_data=shock read_data=kilonova_uvboost'
model='model=shock model=shock model=shock model=kilonova model=kilonova model=kilonova_uvboost model=kilonova_uvboost' 

# ---- bands for figure 2 ---- #
include_optical='include_optical=r include_optical=False include_optical=r'
include_uv='include_uv=NUV_D include_uv=NUV_D include_uv=False' #"include_uv=D1 include_uv=D1,D2"
delay=0

# ---- delays for figure 6 ---- #
#delay='delay=0 delay=2 delay=4 delay=12'
#include_optical=False  # include_optical='r' 'False'
#include_uv=NUV_D  # 'NUV_D' 'False'

export read_data model method dist delay print_progress include_optical include_uv sample resume_previous save_after_seconds parallel dlogz_threshold redden optical_delay

# ::: means combine all possible combinations, :::+ means matched arguments.
# ---- jobs figure 2 ---- #

# produce data
parallel poetry run python ./produce-data/produce-data.py ::: $dist ::: $read_data :::+ $model ::: $include_optical :::+ $include_uv # produce necessary data (perhaps not necessary)
# dry run (flag: --dry-run)
parallel --dry-run poetry run python parameter_estimation.py ::: $dist ::: $read_data :::+ $model ::: $include_optical :::+ $include_uv
# run (useful flag: --progress)
parallel --progress poetry run python parameter_estimation.py ::: $dist ::: $read_data :::+ $model ::: $include_optical :::+ $include_uv

# ---- jobs figure 6 ---- #

# produce data
# parallel poetry run python ./produce-data/produce-data.py ::: $dist ::: $read_data :::+ $model ::: $delay # produce necessary data (perhaps not necessary)
# dry run (flag: --dry-run)
# parallel --dry-run poetry run python parameter_estimation.py ::: $dist ::: $read_data :::+ $model ::: $delay
# run (useful flag: --progress)
# parallel --progress poetry run python parameter_estimation.py ::: $dist ::: $read_data :::+ $model ::: $delay
