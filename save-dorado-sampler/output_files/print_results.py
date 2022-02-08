import pickle
import numpy as np


def read_results(resultsstring):
    with open(resultsstring,'rb') as analysis_results: 
        results = pickle.load(analysis_results)
    logz = results.logz[-1]
    log10z = logz*np.log10(np.e)
    print()
    print(resultsstring)
    results.summary()
    print('log10z',log10z)

# testing
model = 'kilonova'
read_data = 'shock'
dist = '40'
timestamp = '22-02-07 2027'


resultsstring = f"{model}model_{read_data}data_{dist}Mpc/{timestamp} results 12"
read_results(resultsstring)
