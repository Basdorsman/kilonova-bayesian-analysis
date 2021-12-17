import corner
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

def cornerplot(model, samples, legend_texts, colors, smoother=2, plot_density=True, plot_datapoints=False, no_fill_contours=False, fill_contours=False ,quantiles=(0.16,0.5,0.84), levels=(1-np.exp(-0.5),1-np.exp(-2)),bins=20):
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

def get_samples(files):
    samples = []
    for file in files:
        with open(file,'rb') as analysis_results:
            samples.append(pickle.load(analysis_results).samples)
    return(samples)


model = 'kilonova' # kilonova, kilonova_uvboost, shock
plottype = 'uvoptical' #uvoptical, uv


# model presets
if plottype == 'uv':
    data = 'no_opticalband_uv'    
elif plottype == 'uvoptical':
    distance = '160'
    data = f'uv_vs_uvoptical_{distance}mpc'# no_opticalband_uv, uv_vs_uvoptical_{distance}mpc
modeldata = f'{model}_{data}'


if modeldata == 'kilonova_no_opticalband_uv':
    files = ["kilonovamodel_kilonovadata_40Mpc/21-05-26 1442 results",
             "kilonovamodel_kilonovadata_100Mpc/21-05-26 1526 results",
             "kilonovamodel_kilonova_opticaldata_160Mpc_no_opticalband_uv/21-10-04 1224 results"][::-1]
    legend_texts = [ 'Data at 40 Mpc','Data at 100 Mpc','Data at 160 Mpc'][::-1]
    colors = ['blue','gold','red']
elif modeldata == 'kilonova_uvboost_no_opticalband_uv':
    files = [f"{model}model_{model}_opticaldata_40Mpc_{data}/21-10-21 0959 results",
             f"{model}model_{model}_opticaldata_100Mpc_{data}/21-10-21 0905 results",
             f"{model}model_{model}_opticaldata_160Mpc_{data}/21-10-04 1443 results"][::-1]
    legend_texts = [ 'Data at 40 Mpc','Data at 100 Mpc','Data at 160 Mpc'][::-1]
    colors = ['blue','gold','red']
elif modeldata == 'kilonova_uvboost_uv_vs_uvoptical_40mpc':
    files = [f"{model}model_{model}_opticaldata_{distance}Mpc_no_opticalband_uv/21-10-21 0959 results",
             f"{model}model_{model}_opticaldata_{distance}Mpc_rband_uv/21-10-21 1228 results"]
    legend_texts = [ f'Data at {distance} Mpc UV',f'Data at {distance} Mpc UV + Optical']
    colors = ['blue','red']
elif modeldata == 'kilonova_uvboost_uv_vs_uvoptical_100mpc':
    files = [f"{model}model_{model}_opticaldata_{distance}Mpc_no_opticalband_uv/21-10-22 1431 results",
             f"{model}model_{model}_opticaldata_{distance}Mpc_rband_uv/21-11-01 1104 results"]
    legend_texts = [ f'Data at {distance} Mpc UV',f'Data at {distance} Mpc UV + Optical']
    colors = ['blue','red']
elif modeldata == 'kilonova_uvboost_uv_vs_uvoptical_160mpc':
    files = [f"{model}model_{model}_opticaldata_{distance}Mpc_no_opticalband_uv/21-10-04 1443 results",
             f"{model}model_{model}_opticaldata_{distance}Mpc_rband_uv/21-11-01 1212 results"]
    legend_texts = [ f'Data at {distance} Mpc UV',f'Data at {distance} Mpc UV + Optical']
    colors = ['blue','red']
elif modeldata == 'kilonova_uv_vs_uvoptical_40mpc':
    files = [f"{model}model_{model}_opticaldata_{distance}Mpc_no_opticalband_uv/21-09-10 1409 results",
             f"{model}model_{model}_opticaldata_{distance}Mpc_rband_uv/21-08-30 1032 results"]
    legend_texts = [ f'Data at {distance} Mpc UV',f'Data at {distance} Mpc UV + Optical']
    colors = ['blue','red']
elif modeldata == 'kilonova_uv_vs_uvoptical_100mpc':
    files = [f"{model}model_{model}_opticaldata_100Mpc_no_opticalband_uv/21-10-22 1351 results",
             f"{model}model_{model}_opticaldata_100Mpc_rband_uv/21-10-22 0910 results"]
    legend_texts = [ f'Data at {distance} Mpc UV',f'Data at {distance} Mpc UV + Optical']
    colors = ['blue','red']
elif modeldata == 'kilonova_uv_vs_uvoptical_160mpc':
    files = [f"{model}model_{model}_opticaldata_{distance}Mpc_no_opticalband_uv/21-10-04 1224 results",
             f"{model}model_{model}_opticaldata_{distance}Mpc_rband_uv/21-10-22 1042 results"]
    legend_texts = [ f'Data at {distance} Mpc UV',f'Data at {distance} Mpc UV + Optical']
    colors = ['blue','red']
elif modeldata == 'shock_no_opticalband_uv':
    files = [f"{model}model_{model}_opticaldata_40Mpc_{data}/21-11-01 1036 results",
             f"{model}model_{model}_opticaldata_100Mpc_{data}/21-11-01 1026 results",
             f"{model}model_{model}_opticaldata_160Mpc_{data}/21-11-01 1041 results"][::-1]
    legend_texts = [ 'Data at 40 Mpc','Data at 100 Mpc','Data at 160 Mpc'][::-1]
    colors = ['blue','gold','red']
elif modeldata == 'shock_uv_vs_uvoptical_100mpc':
    files = [f"{model}model_{model}_opticaldata_{distance}Mpc_no_opticalband_uv/21-11-01 1026 results",
             f"{model}model_{model}_opticaldata_{distance}Mpc_rband_uv/21-11-01 1324 results"]
    legend_texts = [ f'Data at {distance} Mpc UV',f'Data at {distance} Mpc UV + Optical']
    colors = ['blue','red']    
elif modeldata == 'shock_uv_vs_uvoptical_160mpc':
    files = [f"{model}model_{model}_opticaldata_{distance}Mpc_no_opticalband_uv/21-11-01 1041 results",
             f"{model}model_{model}_opticaldata_{distance}Mpc_rband_uv/21-11-01 1403 results"]
    legend_texts = [ f'Data at {distance} Mpc UV',f'Data at {distance} Mpc UV + Optical']
    colors = ['blue','red']  


samples = get_samples(files)
figure = cornerplot(model,samples, legend_texts, colors)



plt.show()
figure.savefig(f'plots/{modeldata}.png',dpi=300,pad_inches=0.3,bbox_inches='tight')
#figure.savefig('plots/test.png',dpi=300,pad_inches=0.3,bbox_inches='tight')