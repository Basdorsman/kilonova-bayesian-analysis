model='shock'
read_data='shock'
method='test' 
dist=40
delay=0
print_progress=True
max_time=360000
export model read_data method dist delay print_progress max_time
poetry run python ./produce-data/produce-data.py

export include_optical='r'
export include_uv='D1'
poetry run python ./parameter_estimation.py 


# example of nested loop
# for include_optical in 'False' 'r'; do
#    export include_optical
#     for include_uv in 'D1,D2' 'D1'; do
#         export include_uv
#         poetry run python parameter_estimation.py
#     done
# done