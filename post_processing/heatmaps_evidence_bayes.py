import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
import matplotlib.colors as colors
from heatmaps import heatmap, annotate_heatmap
import dill as pickle
import sys
sys.path.append('../')
from dynesty_sampler import find

analysisModels = ['shock','kilonova','kilonova_uvboost']
dataModels = ['shock','kilonova','kilonova_uvboost']
dist = 160

def getLogZ(model,data,dist):
    folderstring = f'../output_files/results/{model}model_{data}data_0h_delay'
    filestring=f'{dist}Mpc_no_opticalband_NUV_Dband'
    try: # Runs with dlogz_threshold=0.5
        with open(folderstring+'/'+filestring+'_results_dlogz=False','rb') as resultsfile:
            results = pickle.load(resultsfile)
        dlogz=0.5
        try:
            logz=results['logz'][-1]
        except:
            logz=results.results['logz'][-1]
    except FileNotFoundError: # Exception to include any newer results with custom dlogz_thresholds.
        if isinstance(find(filestring+'_results_dlogz=*', folderstring), str):
            intermediate_results = find(filestring+'_results_dlogz=*', folderstring)
            #print(intermediate_results)
            with open(intermediate_results,'rb') as file:
                sampler = pickle.load(file)
            dlogz=intermediate_results.split('=')[1]
            logz=sampler.results['logz'][-1]
        else: # No results available.
            dlogz=np.NaN
            logz=np.NaN
    log10z = logz/np.log(10)
    return log10z, dlogz, model, data


logz = np.asarray([[getLogZ(model,data,dist)[0] for model in analysisModels]
                   for data in dataModels])
#logz = np.asarray([[19.67748891, -500, -545.47659009],[-334.40136605, 8.34835139, np.NaN],[-574.52080744, np.NaN, 8.13502999]])
logb = np.asarray([[np.NaN,logz[0,0]-logz[0,1],logz[0,0]-logz[0,2]],
                   [logz[1,1]-logz[1,0],np.NaN,np.NaN],
                   [logz[2,2]-logz[2,0],np.NaN,np.NaN]])
datas = (logz,logb)
norms = (colors.SymLogNorm(linthresh=10, linscale=1, 
                          vmin=-10**4,
                          vmax=10**1, base=10),
        colors.LogNorm(vmin=1,
                        vmax=10**4))

yticklabelvisibles = (True, False)
analysisLabels=['shock','kilonova','kilonova uvboost']
dataLabels=['actual data:\nshock','actual data:\nkilonova','actual data:\nkilonova uvboost']

cmaps= ("RdYlGn","Greens")
fig, axes = plt.subplots(1,2,figsize=(10,5))


for ax, data, norm, cmap, visible in zip(axes, datas, norms, cmaps, yticklabelvisibles):
    im, cbar = heatmap(data, dataLabels, analysisLabels, ax=ax, cmap=cmap,
                       norm=norm,
                    cbar_kw={'drawedges':False, 'pad':0.01, 'shrink':0.75},
                    cbarlabel="Evidence: Log$_{10}$($\mathcal{Z}$)",
                    yticklabelvisible=visible)
    annotate_heatmap(im, data, threshold=0.5, valfmt="{x:.1f}")

fig.tight_layout()
plt.show()
fig.savefig(f'./plots/heatmap_{dist}Mpc.png',dpi=300,pad_inches=0.3,bbox_inches='tight')
