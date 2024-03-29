# https://matplotlib.org/3.4.3/gallery/images_contours_and_fields/image_annotated_heatmap.html#sphx-glr-gallery-images-contours-and-fields-image-annotated-heatmap-py
# https://matplotlib.org/stable/tutorials/colors/colormapnorms.html

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import math
import dill as pickle
import sys
sys.path.append('../')
from dynesty_sampler import find

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", yticklabelvisible=True, **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
    cbar.outline.set_visible(False)

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=0, ha="center",
             rotation_mode="anchor")
    plt.setp(ax.get_yticklabels(), rotation=80, va="bottom", ha="center",
             rotation_mode="anchor",visible=yticklabelvisible)

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar

def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=0.5, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # # Normalize the threshold to the images color range.
    # if threshold is not None:
    #     threshold = im.norm(threshold)
    # else:
    #     threshold = im.norm(data[~np.isnan(data)].max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if not np.isnan(data[i,j]):  # Skip if no data.
                kw.update(color=textcolor(data[i,j],im.norm,threshold=threshold))
                #kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
                text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
                texts.append(text)
    return texts

def textcolor(value, norm, threshold=(0.75, 0.75), textcolors=("black", "white")):
    if type(norm)==colors.SymLogNorm:
        base=norm._scale.base
        if norm.vmax>0:
            if(value < -base**(threshold[0]*math.log(-norm.vmin, base)) or value > base**(threshold[0]*math.log(norm.vmax, base))): 
                return textcolors[1]
            else:
                return textcolors[0]
        elif norm.vmax<0:    
            if(value < -base**(threshold[0]*math.log(-norm.vmin, base)) or value > -base**(threshold[0]*math.log(-norm.vmax, base))):
                return textcolors[1]
            else:
                return textcolors[0]
    elif type(norm)==colors.LogNorm:
        base=norm._scale.base
        if(value < base**(math.log(norm.vmax+norm.vmin, base)*threshold[1])):
           return textcolors[0]
        else:
           return textcolors[1]
       
def getLogZ(model,data,dist,optical_band='no_optical', uv_band='NUV_D', delay=0, optical_delay=12, redden=False):
    folderstring = f'../output_files/results/{model}model_{data}data_{delay}h_delay'
    filestring=f'{dist}Mpc_{optical_band}band_{uv_band}band'
    if redden: # Just a band-aid. There are old versions and new versions of data that may or may not include these bits of string in their names.
        filestring+=f'_redden_{redden}'
    if optical_delay:
        filestring+=f'_optical_delay_{optical_delay}'
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
        else: 
            print(f'No results available for {model}model {data}data dist={dist} delay={delay} redden={redden}')
            dlogz=np.NaN
            logz=np.NaN
    log10z = logz/np.log(10)
    return log10z, dlogz, model, data

def bayesPlot(fig, ax, models=['shock','kilonova_uvboost','kilonova'], datas=['shock','kilonova_uvboost','kilonova'], delay=0, redden=False, dists=[40, 100, 160], optical_band='no_optical', uv_band='NUV_D', legend_labels = ['kilonova L, shock data','kilonova D, shock data','shock, kilonova L data','shock, kilonova D data'], linestyles=['-',':','--','-.'], **kwargs):
    logz = np.asarray([[[getLogZ(model,data,dist,optical_band=optical_band,uv_band=uv_band, delay=delay, redden=redden)[0] for dist in dists] for model in models] for data in datas])
    if len(models)==3 and len(datas)==3: 
        logb = np.asarray([[np.NaN,logz[0,0]-logz[0,1],logz[0,0]-logz[0,2]],[logz[1,1]-logz[1,0],np.NaN,np.NaN],[logz[2,2]-logz[2,0],np.NaN,np.NaN]],dtype=object)
        logb_fields = [(0,1),(0,2),(1,0),(2,0)]
    elif len(models)==2 and len(datas)==2:
        logb = np.asarray([logz[0,0]-logz[0,1], logz[1,1]-logz[1,0]])
        logb_fields = [0, 1]
    elif len(models)==2 and len(datas)==1:
        if datas == ['shock']:
            logb = np.asarray([logz[0,0]-logz[0,1]])
        else:
            logb = np.asarray([logz[0,1]-logz[0,0]])
        logb_fields = [0]
    for field, label, linestyle in zip(logb_fields, legend_labels, linestyles):
        ax.plot(dists,logb[field], label=label, linestyle=linestyle, **kwargs) #this is bad, since now I can't make different types of plots using this function. instead: return datapoints for plotting and plot in a separate function.
    return fig, ax

def BayesPlotRainbow(fig, ax, kilonovamodel, data, delays, redden=True, dists=[40, 100, 160], optical_band='no_optical', uv_band='NUV_D', colors=['red','orange','yellow','blue']):
    datas=[data]
    models=['shock', kilonovamodel]
    logbs = []
    for delay in delays:
        logz = np.asarray([[[getLogZ(model,data,dist,optical_band=optical_band,uv_band=uv_band, delay=delay, redden=redden)[0] for dist in dists] for model in models] for data in datas])
        if datas == ['shock']:
            logb = logz[0,0]-logz[0,1]
        else:
            logb = logz[0,1]-logz[0,0]
        logbs.append(logb)
        uppers = logbs
        lowers = logbs[1:]
        lowers.append(-1)
    for upper, lower, color in zip(uppers, lowers, colors):
        ax.fill_between(dists, upper, lower, color=color)

def getBayesData(models=['shock','kilonova_uvboost','kilonova'], datas=['shock','kilonova_uvboost','kilonova'], delay=0, optical_delays=[3, 6, 12], redden=False, dists=70, optical_band='no_optical', uv_band='NUV_D', legend_labels = ['kilonova L, shock data','kilonova D, shock data','shock, kilonova L data','shock, kilonova D data'], linestyles=['-',':','--','-.'], **kwargs):
    logz = np.asarray([[[getLogZ(model,data,dist,optical_band=optical_band,uv_band=uv_band, optical_delay=optical_delay, delay=delay, redden=redden)[0] for optical_delay in optical_delays] for model in models] for data in datas])
    if len(models)==3 and len(datas)==3:
        logb = np.asarray([[np.NaN,logz[0,0]-logz[0,1],logz[0,0]-logz[0,2]],[logz[1,1]-logz[1,0],np.NaN,np.NaN],[logz[2,2]-logz[2,0],np.NaN,np.NaN]],dtype=object)
        logb_fields = [(0,1),(0,2),(1,0),(2,0)]
    elif len(models)==2 and len(datas)==2:
        logb = np.asarray([logz[0,0]-logz[0,1], logz[1,1]-logz[1,0]])
        logb_fields = [0, 1]
    elif len(models)==2 and len(datas)==1:
        if datas == ['shock']:
            logb = np.asarray([logz[0,0]-logz[0,1]])
        else:
            logb = np.asarray([logz[0,1]-logz[0,0]])
        logb_fields = [0]
    return logb, logb_fields, legend_labels, linestyles


if __name__ == '__main__':
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"
    
    
    
    # Bayesplot optical delay
    fontsize_bayesplot=15
    optical_band='r'
    uv_band='no_uv'
    delay=0
    optical_delays = [3, 6, 12]
    redden = True
    first_datas=['3','6','12']
    markers = ['o', '*', 'P', 'D']
    bayesColors=['tab:orange','tab:blue','tab:green', 'tab:red']
    dist=70
    kilonovamodels = ['kilonova','kilonova_uvboost']
    bayesTitles = ['Default Kilonova Model', 'Lower Early Opacity Kilonova Model']
    legends = [False, True]
    ylabels = [True, False]

    fig, axes = plt.subplots(1,2, sharey=True, figsize=(11.5,3.5))
    for ax, kilonovamodel, title, legend, ylabel in zip(axes, kilonovamodels, bayesTitles, legends, ylabels):
        logb, logb_fields, legend_labels, linestyles = getBayesData(models=['shock', kilonovamodel], datas=['shock',kilonovamodel], delay=delay, optical_delays=optical_delays, redden=redden, dists=dist, optical_band=optical_band, uv_band=uv_band, legend_labels=['shock data, 70 Mpc','kilonova data, 70 Mpc'])
        for field, label, linestyle in zip(logb_fields, legend_labels, linestyles):
            ax.plot(optical_delays,logb[field], label=label, linestyle=linestyle)
        ax.set_title(f'{title}',fontsize=fontsize_bayesplot)
        ax.set_yscale('symlog', linthresh=2)
        ax.hlines(2,0,12,color='k',label='$\log_{10}(\mathcal{B})=2$')
        ax.set_xlabel('optical delay',fontsize=fontsize_bayesplot)
        if legend:
            ax.legend(bbox_to_anchor=(1,1),loc='upper left',fontsize=fontsize_bayesplot-2)
        if ylabel:
            ax.set_ylabel('$\log_{10}(\mathcal{B})$',fontsize=fontsize_bayesplot-1)
            ax.tick_params(axis='y', labelsize=fontsize_bayesplot-1)
    
    
    
    # https://matplotlib.org/stable/gallery/color/colorbar_basics.html#sphx-glr-gallery-color-colorbar-basics-py
    # Maybe even better idea: back to color maps. x-axis: distance, y-axis: delay, color: evidence.
    
    # BayesPlot rainbow
    # fontsize_bayesplot=15
    # optical_band='no_optical'
    # uv_band='NUV_D'
    # delays = [0, 2, 4, 12]
    # redden = True
    # first_datas=['1.2','3.2','5.2','13.2']
    # markers = ['o', '*', 'P', 'D']
    # bayesColors=['tab:orange','tab:blue','tab:green', 'tab:red']
    # dists=[40, 70, 100, 130, 160]
    # datass = [['shock', 'kilonova'], ['shock', 'kilonova_uvboost']] 
    # kilonovamodels = ['kilonova','kilonova_uvboost']
    # #bayesTitles = ['Default Kilonova Model', 'Lower Early Opacity Kilonova Model']
    # legendss = [[False, True], [False, False]]
    # ylabels = [True, False]

    # fig, axess = plt.subplots(2,2, sharey=True, figsize=(11.5,7))
    # for axes, kilonovamodel, datas, legends, ylabels in zip(axess, kilonovamodels, datass, legendss, ylabels):
    #     for ax, data, legend in zip(axes, datas, legends):
    #         # for delay, color in zip(delays, bayesColors):
    #         #     fig, ax = bayesPlot(fig, ax, models=['shock', kilonovamodel], datas=[data], delay=delay, redden=redden, dists=dists, optical_band=optical_band, uv_band=uv_band, color=color)
    #         #     ax.set_title(f'data={data}, kilonovamodel={kilonovamodel}',fontsize=fontsize_bayesplot)
    #         #     ax.set_yscale('symlog', linthresh=2)
    #         #     ax.hlines(2,40,160,color='k',label='$\log_{10}(\mathcal{B})=2$')
    #         BayesPlotRainbow(fig, ax, kilonovamodel, data, delays, redden=redden, dists=dists, optical_band=optical_band, uv_band=uv_band)
    #         ax.set_title(f'data={data}, kilonovamodel={kilonovamodel}',fontsize=fontsize_bayesplot)
    #         ax.set_yscale('symlog', linthresh=2)
    #         ax.hlines(2,40,160,color='k',label='$\log_{10}(\mathcal{B})=2$')
    #         ax.set_ylim(bottom=0)
            
    
    # BayesPlot
    # fontsize_bayesplot=15
    # optical_band='no_optical'
    # uv_band='NUV_D'
    # delays = [0, 2, 4, 12]
    # redden = True
    # first_datas=['1.2','3.2','5.2','9.2']
    # markers = ['o', '*', 'P', 'D']
    # bayesColors=['tab:orange','tab:blue','tab:green', 'tab:red']
    # dists=[40, 70, 100, 130, 160]
    # kilonovamodels = ['kilonova','kilonova_uvboost']
    # bayesTitles = ['Default Kilonova Model', 'Lower Early Opacity Kilonova Model']
    # legends = [False, True]
    # ylabels = [True, False]

    # fig, axes = plt.subplots(1,2, sharey=True, figsize=(11.5,3.5))
    # for ax, kilonovamodel, title, legend, ylabel in zip(axes, kilonovamodels, bayesTitles, legends, ylabels):
    #     for delay, first_data, marker, color in zip(delays, first_datas, markers, bayesColors):
    #         fig, ax = bayesPlot(fig, ax, models=['shock', kilonovamodel], datas=['shock',kilonovamodel], delay=delay, redden=redden, dists=dists, optical_band=optical_band, uv_band=uv_band, legend_labels=[f'Shock Data, Data @ {first_data}h',f'Kilonova Data, Data @ {first_data}h'], marker=marker, color=color)
    #     ax.set_title(f'{title}',fontsize=fontsize_bayesplot)
    #     ax.set_yscale('symlog', linthresh=2)
    #     ax.hlines(2,40,160,color='k',label='$\log_{10}(\mathcal{B})=2$')
    #     ax.set_xlabel('Luminosity Distance (Mpc)',fontsize=fontsize_bayesplot)
    #     if legend:
    #         ax.legend(bbox_to_anchor=(1,1),loc='upper left',fontsize=fontsize_bayesplot-2)
    #     if ylabel:
    #         ax.set_ylabel('$\log_{10}(\mathcal{B})$',fontsize=fontsize_bayesplot-1)
    #         ax.tick_params(axis='y', labelsize=fontsize_bayesplot-1)
    
    # Heatmaps
    # cmaps= ("RdYlBu","Greens")
    # yticklabelvisibles = (True, False)
    
    # analysisModels = ['shock','kilonova_uvboost','kilonova']
    # analysisLabels=['Shock','Kilonova Lower\nEarly Opacity','Kilonova Default']
    
    # dataModels = ['shock','kilonova_uvboost','kilonova']
    # dataLabels=['Data: Shock','Data: Kilonova\nLower Early Opacity','Data: Kilonova']
    
    # cbarlabels=["Evidence: Log$_{10}$($\mathcal{Z}$)","Bayes' factor: Log$_{10}$($\mathcal{B}$)"]
    # dist = 100
    # logz = np.asarray([[getLogZ(model,data,dist,optical_band=optical_band,uv_band=uv_band)[0] for model in analysisModels]
    #                    for data in dataModels])
    # logb = np.asarray([[np.NaN,logz[0,0]-logz[0,1],logz[0,0]-logz[0,2]],
    #                    [logz[1,1]-logz[1,0],np.NaN,np.NaN],
    #                    [logz[2,2]-logz[2,0],np.NaN,np.NaN]])
    # datas = (logz,logb)
    # norms = (colors.SymLogNorm(linthresh=10, linscale=1, 
    #                           vmin=-10**4,
    #                           vmax=10**1, base=10),
    #         colors.LogNorm(vmin=1,
    #                         vmax=10**4))
    
    # fig, axes = plt.subplots(1,2,figsize=(10,5))
    # for ax, data, norm, cmap, visible, cbarlabel in zip(axes, datas, norms, cmaps, yticklabelvisibles, cbarlabels):
    #     im, cbar = heatmap(data, dataLabels, analysisLabels, ax=ax, cmap=cmap,
    #                        norm=norm,
    #                     cbar_kw={'drawedges':False, 'pad':0.01, 'shrink':0.75},
    #                     cbarlabel=cbarlabel,
    #                     yticklabelvisible=visible)
    #     annotate_heatmap(im, data, threshold=(0.75, 0.75), valfmt="{x:.1f}", fontsize=14)
    
    fig.tight_layout()
    plt.show()