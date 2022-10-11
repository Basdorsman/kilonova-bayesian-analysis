delay=0
dist=40
read_data=shock
export delay dist read_data

poetry run python produce-data/produce-data.py
