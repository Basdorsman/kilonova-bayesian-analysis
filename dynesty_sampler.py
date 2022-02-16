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
        print('Initiating new multipool sampler')
        with MultiPool() as pool:
            sampler = NestedSampler(loglikelihood, priortransform, ndim, pool=pool, sample=sample)
            print('poolsize = ',pool.size)
    else:
        print('Initiating new sampler')
        sampler = NestedSampler(loglikelihood, priortransform, ndim, sample=sample)
    return sampler

def getSampler(ndim, folderstring, filestring, loglikelihood=None, priortransform=None, parallel=True, sample='auto', resume_previous=True):
    if resume_previous:
        if isinstance(find(filestring+'_sampler_dlogz=*', folderstring), str):
            intermediate_output = find(filestring+'_sampler_dlogz=*', folderstring)
            print('Using previous sampler: '+intermediate_output)
            with open(intermediate_output,'rb') as samplerfile:
                sampler = pickle.load(samplerfile)
            with open(intermediate_output+'_rstate','rb') as rstatefile:
                sampler.rstate = pickle.load(rstatefile)
            previous_dlogz=intermediate_output.split('=')[1]
        else:
            sampler = initiateSampler(loglikelihood, priortransform, ndim, parallel=parallel, sample=sample)
            previous_dlogz=False
    else:
        sampler = initiateSampler(loglikelihood, priortransform, ndim, parallel=parallel, sample=sample)
        previous_dlogz=False
    return sampler, previous_dlogz

def externalSamplingLoop(sampler, folderstring, filestring, previous_dlogz=False, sample='auto', save_after_seconds=600, print_progress=True, dlogz_threshold=0.5):
    sample_start = time.time()
    print('sampler running...')
    for it, res in enumerate(sampler.sample(dlogz=dlogz_threshold)):
        if print_progress:
            print(f'it: {it}, nc: {res[9]} ,delta_logz: {res[-1]}')
        if save_after_seconds:
            if int(np.ceil(time.time()-sample_start))>save_after_seconds:
                with open(folderstring+'/'+filestring+f'_sampler_dlogz={res[-1]}','wb') as samplerfile :
                    pickle.dump(sampler,samplerfile)
                with open(folderstring+'/'+filestring+f'_sampler_dlogz={res[-1]}_rstate','wb') as rstatefile :
                    pickle.dump(sampler.rstate,rstatefile)
                if print_progress:
                    print(f'saved sampler at dlogz = {res[-1]}')
                if previous_dlogz:
                    os.remove(folderstring+'/'+filestring+f'_sampler_dlogz={previous_dlogz}')
                    os.remove(folderstring+'/'+filestring+f'_sampler_dlogz={previous_dlogz}_rstate')
                    if print_progress:
                        print(f'removed old sampler at dlogz = {previous_dlogz}')
                previous_dlogz=res[-1]
                sample_start = time.time()
    for it_final, res in enumerate(sampler.add_live_points()):
        pass
    print('added live points')
    with open(folderstring+'/'+filestring+'_results_test','wb') as samplerfile :
        pickle.dump(sampler,samplerfile)
    print('saved results')

def wrappedSampler(sampler, folderstring, filestring, previous_dlogz=False, sample='auto', save_after_seconds=60, print_progress=True, parallel=True, dlogz_threshold=0.5):
    if parallel:
        with MultiPool() as pool:
            sampler.pool = pool
            sampler.queue_size = pool.size
            sampler.use_pool = {}
            sampler.use_pool_evolve = True
            sampler.use_pool_logl = True
            sampler.use_pool_ptform = True
            sampler.use_pool_update = True
            sampler.M = pool.map
            externalSamplingLoop(sampler, folderstring, filestring, previous_dlogz=previous_dlogz, sample=sample, save_after_seconds=save_after_seconds, print_progress=print_progress,dlogz_threshold=dlogz_threshold)
    else:
        externalSamplingLoop(sampler, folderstring, filestring, previous_dlogz=previous_dlogz, sample=sample, save_after_seconds=save_after_seconds, print_progress=print_progress,dlogz_threshold=dlogz_threshold)

        
        
if __name__ == '__main__':
    filestring = '40Mpc_no_opticalband_NUV_Dband'
    folderstring = 'output_files/results/kilonova_uvboostmodel_shockdata_0h_delay'
    intermediate_output = find(filestring+'_sampler_dlogz=*', folderstring)
    print('opened file: '+intermediate_output)
    with open(intermediate_output,'rb') as samplerfile:
        sampler = pickle.load(samplerfile)
