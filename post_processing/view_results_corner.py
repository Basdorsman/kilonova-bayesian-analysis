import corner
import dill as pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

def cornerplot(model, samples, legend_texts, colors, smoother=2, plot_density=True, plot_datapoints=False, no_fill_contours=False, fill_contours=False ,quantiles=(0.16,0.5,0.84), levels=(1-np.exp(-0.5),1-np.exp(-2)),bins=20):
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
        plot_density (bool): varialbe for Corner.
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
        fontsize=19
    # shock model
    elif model == 'shock':
        k = 5 # 0.1 cm^2/g
        m = 1 #0.01 solar masses
        v = 2 #0.1c
        r = 5 #10^10 cm #Initial radius for shock
        theta_truths = k, m, v, r
        labels = [r'$\kappa$', r'$M_\mathrm{shock}$', r'$v_\mathrm{shock}$', r'$R_\mathrm{shock}$']
        fontsize=15

    label_kwargs = {'fontsize': fontsize}
    title_kwargs = label_kwargs
    # first figure
    figure = corner.corner(samples[0], bins=bins, labels=labels,label_kwargs=label_kwargs, smooth=smoother, color=colors[0],truths=theta_truths,truth_color='k', plot_datapoints=plot_datapoints,plot_density=plot_density,no_fill_contours=no_fill_contours,fill_contours=fill_contours,levels=levels,show_titles=True, title_kwargs=title_kwargs,quantiles=quantiles)
    titledplot_indices = [i*(len(labels)+1) for i in range(samples[0].shape[1])]
    labels_strip = [labels[i] + ' = ' for i in range(samples[0].shape[1])]
    titles = []
    titles.append([figure.axes[titledplot_index].title._text.replace(label,'') for titledplot_index,label in zip(titledplot_indices, labels_strip)])
    
    # stack remaining figures
    for sample,color in zip(samples[1:],colors[1:]):
        corner.corner(sample, fig=figure, smooth=smoother, color=color, plot_datapoints=plot_datapoints,plot_density=plot_density,no_fill_contours=no_fill_contours,fill_contours=fill_contours,levels=levels,show_titles=True,title_kwargs=title_kwargs,quantiles=quantiles)
        titles.append([figure.axes[titledplot_index].title._text for titledplot_index in titledplot_indices])

    # color and thickness vlines
    truth_line = 3
    medians = [1]
    for i in range(len(samples[1:])):
        medians.append(3*i+5)
    for titledplot_index in titledplot_indices:
        for median in medians:
            figure.axes[titledplot_index].lines[median].set_linestyle('-')
        for lines in range(len(figure.axes[titledplot_index].lines)):
            figure.axes[titledplot_index].lines[lines].set_linewidth(0.75)
        figure.axes[titledplot_index].lines[truth_line].set_linewidth(1.5)

    # set ax.titles
    for parameter_index in range(samples[0].shape[1]):
        titledplot_index = titledplot_indices[parameter_index]
        quantile_strings = ''
        for title,color in zip(titles,colors):
            quantile_strings += '\n'+title[parameter_index]+' ('+color+')'
        figure.axes[titledplot_index].title._text = labels[parameter_index] + ' = ' + str(theta_truths[parameter_index]) + quantile_strings

    legend_handles = []
    
    for legend_text, color in zip(legend_texts,colors):
        legend_handles.append(mlines.Line2D([], [], color=color, markersize=fontsize, label=legend_text))
    figure.legend(handles=legend_handles,loc='upper center',fontsize=fontsize)
    
    for ax in figure.get_axes():
        ax.tick_params(axis='both', labelsize=fontsize)
        
    return figure

def get_samples(models, datas, delays, distances, opticalbands, uvbands):
    '''Fetches results from various samplings.
    
    This function looks up requested sampling results and fetches results. All
    parameters are strings or list of strings. But only one in total should be
    a list of strings.
    
    Parameters:
        models (str or list): analysis models.
        datas (str or list): data models.
        delays (str or list): hours of delay in dataset.
        distances (str or list): distances to event in Mpc.
        opticalbands (str or list): choices for optical bands.
        uvbands (str or list): choices for uv bands.
    
    Returns:
        samples: (list of numpy arrays): samples.
        legend_texts (list): whichever parameter that was a list.
        colors (str): colors, same length as the parameter that was a list, up 
                      to three. 
    '''
    variablesForString = [models, datas, delays, distances, opticalbands, uvbands]
    variableCount = 0
    for variablesList in models, datas, delays, distances, opticalbands, uvbands:
        if isinstance(variablesList,list):
            files = []
            legend_texts = variablesList
            colors = ['blue','gold','red'][:len(variablesList)]
            for variable in variablesList:
                variablesForString[variableCount] = variable
                model, data, delay, distance, opticalband, uvband = variablesForString
                files.append(f'../output_files/results/{model}model_{data}data_{delay}_delay/{distance}Mpc_{opticalband}band_{uvband}band_results_dlogz=False')
        variableCount += 1
    samples = []
    for file in files:
        with open(file,'rb') as analysis_results:
            samples.append(pickle.load(analysis_results).results.samples)
    return(samples,legend_texts,colors)


model = 'kilonova_uvboost' # kilonova, kilonova_uvboost, shock
datas = 'kilonova_uvboost'
delay = '0h'
distance = ['160','100','40']
opticalband = 'no_optical'
uvband= 'NUV_D'
samples,legend_texts,colors = get_samples(model,datas,delay,distance,opticalband,uvband)
legend_texts = ['160 Mpc','100 Mpc','40 Mpc']
figure = cornerplot(model, samples, legend_texts, colors)

plt.show()
figure.savefig(f'plots/cornerplot_{model}.png',dpi=300,pad_inches=0.3,bbox_inches='tight')