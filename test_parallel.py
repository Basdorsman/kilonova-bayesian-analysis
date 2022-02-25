from schwimmbad import MultiPool
from dynesty.dynesty import NestedSampler
import numpy as np

# Define the dimensionality of our problem.
ndim = 3

# Define our 3-D correlated multivariate normal log-likelihood.
C = np.identity(ndim)
C[C==0] = 0.95
Cinv = np.linalg.inv(C)
lnorm = -0.5 * (np.log(2 * np.pi) * ndim +
                np.log(np.linalg.det(C)))

def loglike(x):
    return -0.5 * np.dot(x, np.dot(Cinv, x)) + lnorm

# Define our uniform prior via the prior transform.
def ptform(u):
    return 20. * u - 10.


print('initialization')
with MultiPool(processes=10) as pool:
    print(pool)    
    sampler = NestedSampler(loglike, ptform, ndim,pool=pool)

print('running')
with MultiPool(processes=10) as pool:
    print(pool)
    sampler.pool = pool
    sampler.queue_size = pool.size
    sampler.use_pool = {}
    sampler.use_pool_evolve = True
    sampler.use_pool_logl = True
    sampler.use_pool_ptform = True
    sampler.use_pool_update = True
    sampler.M = pool.map
    sampler.run_nested()
