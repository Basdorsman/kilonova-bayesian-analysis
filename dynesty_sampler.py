from dynesty.dynesty import NestedSampler
import dill as pickle
import time
import numpy as np
import os, fnmatch
from schwimmbad import MultiPool

def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    try:
        lowest_result = min(result)
        return lowest_result
    except:
        return []
        
def initiateSampler(loglikelihood, priortransform, ndim, parallel=True, sample='auto'):
    if parallel:
        print('initiating multipool sampler')
        with MultiPool() as pool:
            sampler = NestedSampler(loglikelihood, priortransform, ndim, pool=pool, sample=sample)
            print('poolsize = ',pool.size)
    elif parallel==False:
        print('initiating sampler')
        sampler = NestedSampler(loglikelihood, priortransform, ndim, sample=sample)
    return sampler

def getSampler(ndim, folderstring, filestring, loglikelihood=None, priortransform=None, parallel=True, sample='auto', intermediate_outputs=True):
    if intermediate_outputs:
        if isinstance(find(filestring+'_sampler_dlogz=*', folderstring), str):
            intermediate_output = find(filestring+'_sampler_dlogz=*', folderstring)
            print('opened file: '+intermediate_output)
            with open(intermediate_output,'rb') as samplerfile:
                sampler = pickle.load(samplerfile)
            with open(intermediate_output+'_rstate','rb') as rstatefile:
                sampler.rstate = pickle.load(rstatefile)
        else:
            sampler = initiateSampler(loglikelihood, priortransform, ndim, parallel=parallel, sample=sample)
    else:
        sampler = initiateSampler(loglikelihood, priortransform, ndim, parallel=parallel, sample=sample)
    return sampler

def wrappedSampler(sampler, loglikelihood, priortransform, ndim, folderstring, filestring, sample='auto', intermediate_outputs=True, save_after_seconds=60, print_progress=True, parallel=True, dlogz_threshold=0.5):
    if intermediate_outputs:
        previous_dlogz=False
        sample_start = time.time()
        with MultiPool() as pool:
            sampler.pool = pool
            sampler.queue_size = pool.size
            sampler.use_pool = {}
            sampler.use_pool_evolve = True
            sampler.use_pool_logl = True
            sampler.use_pool_ptform = True
            sampler.use_pool_update = True
            sampler.M = pool.map
            for it, res in enumerate(sampler.sample(dlogz=dlogz_threshold)):
                print(f'it: {it}, nc: {res[9]} ,delta_logz: {res[-1]}')
                if int(np.ceil(time.time()-sample_start))>save_after_seconds:
                    with open(folderstring+'/'+filestring+f'_sampler_dlogz={res[-1]}','wb') as samplerfile :
                        pickle.dump(sampler,samplerfile)
                    with open(folderstring+'/'+filestring+f'_sampler_dlogz={res[-1]}_rstate','wb') as rstatefile :
                        pickle.dump(sampler.rstate,rstatefile)
                    print(f'saved sampler at dlogz = {res[-1]}')
                    if previous_dlogz:
                        os.remove(folderstring+'/'+filestring+f'_sampler_dlogz={previous_dlogz}')
                        os.remove(folderstring+'/'+filestring+f'_sampler_dlogz={previous_dlogz}_rstate')
                        print(f'removed old sampler at dlogz = {previous_dlogz}')
                    previous_dlogz=res[-1]
                    sample_start = time.time()
            for it_final, res in enumerate(sampler.add_live_points()):
                pass
            print('added live points')
    elif intermediate_outputs==False:    
        print('internal sampler')
        sampler_start = time.time()
        with MultiPool() as pool:
            sampler.pool = pool
            sampler.queue_size = pool.size
            sampler.use_pool = {}
            sampler.use_pool_evolve = True
            sampler.use_pool_logl = True
            sampler.use_pool_ptform = True
            sampler.use_pool_update = True
            sampler.M = pool.map
            sampler.run_nested(print_progress=print_progress,dlogz=dlogz_threshold)
        sampler_time = int(np.ceil(time.time()-sampler_start))
        print(f'sampling time was {sampler_time} seconds')
    with open(folderstring+'/'+filestring+'_results','wb') as resultsfile :
        pickle.dump(sampler.results,resultsfile)
    print('saved results')
    
    
if __name__ == '__main__':
    filestring = '40Mpc_no_opticalband_NUV_Dband'
    folderstring = 'output_files/results/kilonova_uvboostmodel_shockdata_0h_delay'
    intermediate_output = find(filestring+'_sampler_dlogz=*', folderstring)
    print('opened file: '+intermediate_output)
    with open(intermediate_output,'rb') as samplerfile:
        sampler = pickle.load(samplerfile)