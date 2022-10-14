delay=0
dist=40
read_data='shock'
redden='True'
export delay dist read_data redden

poetry run python produce-data/produce-data.py
