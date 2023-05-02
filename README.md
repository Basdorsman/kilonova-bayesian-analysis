# kilonova-bayesian-analysis
Perform a bayesian analysis on kilonova light curves. 

The kilonova model is based on Hotokezaka & Nakar (2019) (https://arxiv.org/abs/1909.02581, https://github.com/hotokezaka/HeatingRate), see also: https://github.com/Basdorsman/kilonova-heating-rate.
The shock model is based on Piro & Kollmeier (2018) (https://ui.adsabs.harvard.edu/abs/2018ApJ...855..103P/abstract).

Related publication: https://iopscience.iop.org/article/10.3847/1538-4357/acaa9e
Related zenodo repository: https://zenodo.org/record/7540160

This package uses poetry (https://python-poetry.org/).

.sh files are PBS job files, which could come in handy if using a cluster. 

# Quick guide:
1. Run produce-data for either shock or kilonova model, to create a set of Dorado UV, and/or Optical data points.
2. As a foolproof, run parameter_estimation "test" mode, and check the output plots to see if the light curves are correctly calculated. The data points you see here are the data points that will be used for the Bayesian analysis.
3. Run parameter_estimation "timeout" or "pool" mode, to run the bayesian analysis. "timeout" mode runs on a single core but allows to set a timer after which the process will be cut off, and data saved. "pool" mode does not feature a timer, but can run on multiple cores. Doesn't work on all computers, and the efficiency gain from pooling is quite limited.

# Dependencies
* [Astropy]
* [Dorado-scheduling]
* [Dorado-sensitivity]
* [Dynesty]
* [Kilonova-heating-rate]
* [Schwimmbad]
* [GNU Parallel]


[Astropy]: https://www.astropy.org
[Dorado-scheduling]: https://dorado-scheduling.readthedocs.io/en/latest/
[Dorado-sensitivity]: https://pypi.org/project/dorado-sensitivity/
[Dynesty]: https://dynesty.readthedocs.io/en/latest/
[Kilonova-heating-rate]: https://pypi.org/project/kilonova-heating-rate/
[Schwimmbad]: https://schwimmbad.readthedocs.io/en/latest/
[GNU Parallel]: https://www.gnu.org/software/parallel/
