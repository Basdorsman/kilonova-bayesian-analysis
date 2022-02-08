import pickle

#date = '22-02-08 1051'
date='22-02-08 1054'
file = f'kilonovamodel_shockdata_40Mpc/{date} priorlims 12'

with open(file,'rb') as priorlims:
    prior=pickle.load(priorlims)
