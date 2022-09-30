import pickle

#date = '22-02-08 1051'

file = f'results/kilonova_uvboostmodel_shockdata_0h_delay/40Mpc_no_opticalband_NUV_Dband_parameters'

with open(file,'rb') as data:
    information=pickle.load(data)
