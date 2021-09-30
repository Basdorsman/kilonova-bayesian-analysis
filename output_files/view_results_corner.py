import corner
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

#with open("shockmodel_shockdata_40Mpc/21-05-11 0821 results",'rb') as analysis_results: 
#with open("kilonovamodel_kilonova_uvboostdata_40Mpc/21-05-24 1259 results",'rb') as analysis_results: 
with open("kilonovamodel_kilonovadata_40Mpc/21-05-26 1442 results",'rb') as analysis_results: 
#with open("kilonovamodel_kilonova_opticaldata_40Mpc_grband_uv/21-08-17 1956 results",'rb') as analysis_results: 
#with open("kilonovamodel_kilonova_opticaldata_40Mpc_rband_uv/21-08-30 1032 results",'rb') as analysis_results: 
#with open("shockmodel_shock_opticaldata_40Mpc_rband_uv/21-09-06 1201 results",'rb') as analysis_results: 
    samples1 = pickle.load(analysis_results).samples

#with open("shockmodel_shockdata_100Mpc/21-05-12 0050 results",'rb') as analysis_results: 
#with open("kilonovamodel_kilonova_uvboostdata_100Mpc/21-05-25 0945 results",'rb') as analysis_results: 
with open("kilonovamodel_kilonovadata_100Mpc/21-05-26 1526 results",'rb') as analysis_results: 
#with open("kilonovamodel_kilonova_opticaldata_40Mpc_grband_no_uv/21-08-18 0854 results",'rb') as analysis_results: 
#with open("kilonovamodel_kilonova_opticaldata_40Mpc_rband_no_uv/21-08-27 1012 results",'rb') as analysis_results: 
#with open("kilonovamodel_kilonova_opticaldata_40Mpc_no_opticalband_uv/21-09-10 1409 results",'rb') as analysis_results: 
#with open("shockmodel_shock_opticaldata_40Mpc_no_opticalband_uv/21-09-06 1212 results",'rb') as analysis_results: 
    samples2 = pickle.load(analysis_results).samples

if samples1.shape[1] == 7:
    model = 'kilonova'
    ########### KILONOVA TRUTHS ##################
    mass = 0.05
    velocities = np.asarray([0.1, 0.2, 0.4])
    opacities = np.asarray([3.0, 0.5])
    n = 4.5
    theta_truths = np.concatenate((mass, velocities, opacities, n),axis=None)
    labels = [r'$M_\mathrm{ejecta}$',r'$v_\mathrm{min}$',r'$v_\kappa$',
              r'$v_\mathrm{max}$',r'$\kappa_\mathrm{high}$',r'$\kappa_\mathrm{low}$',
              r'$n$']
    fontsize=19
elif samples1.shape[1] == 4: #kilonova
    model = 'shock'
    ########### SHOCK TRUTHS ##############
    k = 5 # 0.1 cm^2/g
    m = 1 #0.01 solar masses
    v = 2 #0.1c
    r = 5 #10^10 cm #Initial radius for shock
    theta_truths = k, m, v, r
    labels = [r'$\kappa$', r'$M_\mathrm{shock}$', r'$v_\mathrm{shock}$', r'$R_\mathrm{shock}$']
    fontsize=15

label_kwargs = {'fontsize': fontsize}
title_kwargs = {'fontsize': fontsize}
smoother = 2
plot_density = True
plot_datapoints = False
no_fill_contours = True
quantiles=(0.16,0.5,0.84)
levels = (1-np.exp(-0.5),1-np.exp(-2)) #1 sigma and 2 sigma levels
figure_blue = corner.corner(samples1, bins=20, labels=labels,label_kwargs=label_kwargs, smooth=smoother, color='blue',truths=theta_truths,truth_color='k', plot_datapoints=plot_datapoints,plot_density=plot_density,no_fill_contours=no_fill_contours,levels=levels,show_titles=True, title_kwargs=title_kwargs,quantiles=quantiles)
title_indices = [i*(len(labels)+1) for i in range(len(labels))]
labels_strip = [labels[i] + ' = ' for i in range(len(labels))]
titles1 = [figure_blue.axes[i].title._text.replace(label,'') for i,label in zip(title_indices, labels_strip)]
figure_red = corner.corner(samples2, fig=figure_blue, smooth=smoother, color='red', plot_datapoints=plot_datapoints,plot_density=plot_density,no_fill_contours=no_fill_contours,levels=levels,show_titles=True,title_kwargs=title_kwargs,quantiles=quantiles)
titles2 = [figure_blue.axes[i].title._text for i in title_indices]

#color quantile lines
for j in title_indices:
    for i in [1,5]:
        figure_blue.axes[j].lines[i].set_linestyle('-')
    for i in range(7):
        figure_blue.axes[j].lines[i].set_linewidth(0.75)
    figure_blue.axes[j].lines[3].set_linewidth(1.5)
    
# set ax.titles
for i in range(len(labels)):
    j = title_indices[i]
    figure_blue.axes[j].title._text = labels[i] + ' = ' + str(theta_truths[i]) + '\n' + titles2[i]+' (red)\n'+titles1[i]+' (blue)'

#legend = ['Data at 40 Mpc, uv only','Data at 40 Mpc, uv and optical']
legend = ['Data at 100 Mpc', 'Data at 40 Mpc']

red_line = mlines.Line2D([], [], color='red',
                          markersize=15, label=legend[0])
blue_line = mlines.Line2D([], [], color='blue',
                          markersize=15, label=legend[1])
figure_blue.legend(handles=[red_line, blue_line],loc='upper center',fontsize=fontsize)

for ax in figure_blue.get_axes():
    ax.tick_params(axis='both', labelsize=fontsize)

#figure_blue.savefig(f'plots/{model}_40mpc_UV_vs_UVO.png',dpi=300,pad_inches=0.3,bbox_inches='tight')
figure_blue.savefig(f'plots/{model}_40_and_100_mpc.png',dpi=300,pad_inches=0.3,bbox_inches='tight')

plt.show()

#figure_blue.savefig('plots/test.png',dpi=300)