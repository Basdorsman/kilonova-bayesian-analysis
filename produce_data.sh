delay=0
dist=40
read_data='kilonova'
redden='True'
optical_delay=12
export delay dist read_data redden optical_delay


poetry run python produce-data/produce-data.py
