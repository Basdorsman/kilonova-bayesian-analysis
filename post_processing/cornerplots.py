import corner
import dill as pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

def cornerPlot(model, samples, legend_texts, colors, linestyles, linewidths, smoother=2, plot_density=False, plot_datapoints=False, no_fill_contours=False, fill_contours=False, quantiles=(0.16,0.5,0.84), medians_visible=False, quantilelines_visible=False, levels=(1-np.exp(-0.5),1-np.exp(-2)), bins=20):
    '''
    Produces cornerplot using specified list of samples.
    
    This function produces a cornerplot with user specified samples, 
    legend_texts, and colors. Various kwargs that corner uses are also
    passed along. First, this function applies the correct fiducial parameter
    values for each of the three models. Second, it correctly stacks multiple
    plots and labels in my prefered way.
    
    Parameters:
        model (str): kilonova, kilonova_uvboost or shock model.
        samples (list): a single set or multiple sets of samples.
        legend_texts (list or str): texts to go in the legen.d
        colors (list or str): colors for each set of samples.
        smoother (int): variable for Corner.
        plot_density (bool): Color density inside 2-D plots.
        plot_datapoints (bool): variable for Corner. 
        no_fill_contours (bool): variable for Corner. 
        fill_contours (bool): variable for Corner.
        quantiles (tuple): variable for Corner. The quantiles for verticle 
                           lines.
        levels (miscellaneous): variable for Corner. The surface outlines.
        
    Returns:
        figure (figure object): the figure.
    
    '''
    # kilonova model
    if model == 'kilonova' or model == 'kilonova_uvboost':
        mass = 0.05
        if model == 'kilonova':
            velocities = np.asarray([0.1, 0.2, 0.4])
            opacities = np.asarray([3.0, 0.5])
        elif model == 'kilonova_uvboost':
            velocities = np.asarray([0.1, 0.2, 0.23])
            opacities = np.asarray([3.0, 0.04])
        n = 4.5
        theta_truths = np.concatenate((mass, velocities, opacities, n),axis=None)
        labels = [r'$M_\mathrm{ejecta}$',r'$v_\mathrm{min}$',r'$v_\mathrm{transition}$',
                  r'$v_\mathrm{max}$',r'$\kappa_\mathrm{high}$',r'$\kappa_\mathrm{low}$',
                  r'$n$']
        fontsize=26
        labelpad=0.09
        label_kwargs = {'fontsize': fontsize}
        title_kwargs = {'fontsize': fontsize-2}
        legend_kwargs = {'fontsize': fontsize}
        axislabel_font = fontsize-5

    # shock model
    elif model == 'shock':
        k = 5 # 0.1 cm^2/g
        m = 1 #0.01 solar masses
        v = 2 #0.1c
        r = 5 #10^10 cm #Initial radius for shock
        theta_truths = k, m, v, r
        labels = [r'$\kappa$', r'$M_\mathrm{shock}$', r'$v_\mathrm{shock}$', r'$R_\mathrm{shock}$']
        fontsize=22
        labelpad = 0.02
        label_kwargs = {'fontsize': fontsize}
        title_kwargs = {'fontsize': fontsize-2}
        legend_kwargs = {'fontsize': fontsize-2}
        axislabel_font = fontsize

    # contour kwargs
    hist_kwargs_list=[{},{},{}]
    contour_kwargs_list=[{},{},{}]
    for color, linestyle, linewidth, hist_kwargs, contour_kwargs in zip(colors, linestyles, linewidths, hist_kwargs_list, contour_kwargs_list):
        hist_kwargs['color']=color
        hist_kwargs['linestyle']=linestyle
        hist_kwargs['linewidth']=linewidth
        contour_kwargs['colors']=color
        contour_kwargs['linestyles']=linestyle
        contour_kwargs['linewidths']=linewidth

    # first figure
    figure = corner.corner(samples[0], bins=bins, labels=labels,label_kwargs=label_kwargs, labelpad=labelpad, smooth=smoother,truths=theta_truths,truth_color='k', plot_datapoints=plot_datapoints,plot_density=plot_density,no_fill_contours=no_fill_contours,fill_contours=fill_contours,levels=levels,show_titles=True, title_kwargs=title_kwargs,quantiles=quantiles,hist_kwargs=hist_kwargs_list[0], contour_kwargs=contour_kwargs_list[0])
    titledplot_indices = [i*(len(labels)+1) for i in range(samples[0].shape[1])]
    labels_strip = [labels[i] + ' = ' for i in range(samples[0].shape[1])]
    titles = []
    titles.append([figure.axes[titledplot_index].title._text.replace(label,'') for titledplot_index,label in zip(titledplot_indices, labels_strip)])
    
    # stack remaining figures
    for sample,hist_kwargs, contour_kwargs in zip(samples[1:],hist_kwargs_list[1:], contour_kwargs_list[1:]):
        corner.corner(sample, fig=figure, smooth=smoother, plot_datapoints=plot_datapoints,plot_density=plot_density,no_fill_contours=no_fill_contours,fill_contours=fill_contours,levels=levels,show_titles=True,title_kwargs=title_kwargs,quantiles=quantiles, hist_kwargs=hist_kwargs, contour_kwargs=contour_kwargs)
        titles.append([figure.axes[titledplot_index].title._text for titledplot_index in titledplot_indices])

    # color and thickness vlines
    truth_line = 3
    medians = [1]
    for i in range(len(samples[1:])):
        medians.append(len(colors)*i+5)
    for titledplot_index in titledplot_indices:
        for lines in range(len(figure.axes[titledplot_index].lines)):
            if quantilelines_visible:
                figure.axes[titledplot_index].lines[lines].set_visible(True)
                figure.axes[titledplot_index].lines[lines].set_linewidth(0.75)
            else:
                figure.axes[titledplot_index].lines[lines].set_visible(False)
        for median in medians:
            if medians_visible:
                figure.axes[titledplot_index].lines[median].set_visible(True)
                figure.axes[titledplot_index].lines[median].set_linestyle('-')
            else:
                figure.axes[titledplot_index].lines[median].set_visible(False)
        # truth line
        figure.axes[titledplot_index].lines[truth_line].set_visible(True)
        figure.axes[titledplot_index].lines[truth_line].set_linewidth(1.5)

    # set ax.titles
    for parameter_index in range(samples[0].shape[1]):
        titledplot_index = titledplot_indices[parameter_index]
        quantile_strings = ''
        for title, linestyle in zip(titles,linestyles):
            quantile_strings += '\n'+title[parameter_index]+' ('+linestyle+')'
        #figure.axes[titledplot_index].set_title(label = labels[parameter_index] + ' = ' + str(theta_truths[parameter_index]) + quantile_strings, loc='left')
        figure.axes[titledplot_index].title._text = labels[parameter_index] + ' = ' + str(theta_truths[parameter_index]) + quantile_strings
        figure.axes[titledplot_index].title._horizontalalignment = 'left'
        figure.axes[titledplot_index].title.set_position((0,0))
        
    legend_handles = []
    
    for legend_text, color, linestyle, linewidth in zip(legend_texts, colors, linestyles, linewidths):
        legend_handles.append(mlines.Line2D([], [], color=color, linestyle=linestyle, linewidth=linewidth, markersize=fontsize, label=legend_text))
    figure.legend(handles=legend_handles,loc='upper center', **legend_kwargs)
    
    for ax in figure.get_axes():
        ax.tick_params(axis='both', labelsize=axislabel_font)
        
    return figure

def getSamples(models, datas, delays, distances, bands, legend_texts=None):
    '''Fetches results from various samplings.
    
    This function looks up requested sampling results and fetches results. All
    parameters are strings or list of strings. But only one in total should be
    a list of strings. I combined the optical and uv bands into a single parameter.
    
    Parameters:
        models (str or list): analysis models.
        datas (str or list): data models.
        delays (str or list): hours of delay in dataset.
        distances (str or list): distances to event in Mpc.
        bands (str or list): choices for optical and uv bands

    Returns:
        samples: (list of numpy arrays): samples.
        legend_texts (list): whichever parameter that was a list.
        colors (str): colors, same length as the parameter that was a list, up 
                      to three. 
    '''
    variablesForString = [models, datas, delays, distances, bands]
    variableCount = 0
    for variablesList in models, datas, delays, distances, bands:
        if isinstance(variablesList,list) and len(variablesList)>1:
            files = []
            if legend_texts==None:
                legend_texts = variablesList
            for variable in variablesList:
                variablesForString[variableCount] = variable
                model, data, delay, distance, band = variablesForString
                files.append(f'../output_files/results/{model}model_{data}data_{delay}_delay/{distance}Mpc_{band}_results_dlogz=False')
        variableCount += 1
    samples = []
    for file in files:
        with open(file,'rb') as analysis_results:
            samples.append(pickle.load(analysis_results).results.samples)
    return(samples,legend_texts)

if __name__ == '__main__':
    model = 'shock' # kilonova, kilonova_uvboost, shock
    datas = 'shock'
    delay = '0h'
    distance = ['160', '100', '40']
    band = 'no_opticalband_NUV_Dband'
    legend_texts= ['160 Mpc','100 Mpc','40 Mpc']
    colors = ['blue','orange','green']
    linestyles = ['solid','dashdot','dashed']
    linewidths = [3,3,3]

    samples,legend_texts = getSamples(model, datas, delay, distance, band, legend_texts=legend_texts)
    figure = cornerPlot(model, samples, legend_texts, colors, linestyles, linewidths)
    
    plt.show()
    figure.savefig(f'plots/cornerplot_{band}_{model}.png',dpi=300,pad_inches=0.3,bbox_inches='tight')