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

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data[~np.isnan(data)].max())/2.

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

def textcolor(value, norm, threshold=0.5, textcolors=("black", "white")):
    if type(norm)==colors.SymLogNorm:
        base=norm._scale.base
        if norm.vmax>0:
            if(value < -base**(threshold*math.log(-norm.vmin, base)) or value > base**(threshold*math.log(norm.vmax, base))): 
                return textcolors[1]
            else:
                return textcolors[0]
        elif norm.vmax<0:    
            if(value < -base**(threshold*math.log(-norm.vmin, base)) or value > -base**(threshold*math.log(-norm.vmax, base))):
                return textcolors[1]
            else:
                return textcolors[0]
    elif type(norm)==colors.LogNorm:
        base=norm._scale.base
        if(value < (norm.vmax+norm.vmin)/2):
           return textcolors[0]
        else:
           return textcolors[1]
       
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


if __name__ == '__main__':
    norm=colors.SymLogNorm(linthresh=1, linscale=1, vmin=-10, vmax=10, base=10)
    value = -11
    print(textcolor(value, norm))