'''

Script to generate slider plots of LLNL EBIT 
microcalorimeter spectra

cjperks
June 7th, 2024

'''

# Modules
from atomic_world.read_LLNL_EBIT.main import read_histos as rh

# Enables automatic reloading of modules
#%reload_ext autoreload
#%autoreload 2

# Make histogram plots
specs, scram = rh.plot_histos(
    x_line = [2179],
    #dE = 5, # [eV]
    scram_method = 'cs',
    num_ax = 4,
    spec_init = 15
    )

from matplotlib.ticker import ScalarFormatter
import scipy.constants as cnt
hc = cnt.h*cnt.c/cnt.e*1e10 # [eV*AA]

formatter = ScalarFormatter()
formatter.set_powerlimits((-2, 2))  
keys = [
    #'7.15 kV', '7.28 kV', '7.41 kV', '7.54 kV', '7.67(highCB) kV', '7.67(lowCB) kV', '7.93 kV', '8.19 kV', '8.27 kV', '8.45 kV',
    '8.71 kV', '8.97 kV', '9.23 kV', '9.49 kV', '9.75 kV',
    '15(lowCB) kV', '15(highCB) kV'
    ]

fig, ax = plt.subplots(figsize = (12,6))
vmax = 0
xlim = [6920, 8380]#[4500, 4600]
lw = 2

for kk in keys:
    #ax.step(
    #    specs[kk]['binlow'],
    #    specs[kk]['hist'],
    #    label = kk
    #    )
    ax.plot(
        #0.5*(specs[kk]['binlow']+specs[kk]['binhigh']),
        hc/(0.5*(specs[kk]['binlow']+specs[kk]['binhigh'])),
        specs[kk]['hist'],
        label = kk,
        linewidth = lw
        )
    ind = np.where(
        (specs[kk]['binlow'] >= xlim[0])
        & (specs[kk]['binlow'] <= xlim[1])
        )   

    vmax = np.max(np.r_[vmax, np.max(specs[kk]['hist'][ind])])

'''
E0 = 4557.591
dE = 5

ax.axvspan(
    E0-dE,
    E0+dE,
    color = 'm',
    alpha = 0.3
    )
ax.plot(
    [E0, E0],
    [0, 1.1*vmax],
    'k--',
    linewidth=lw-1
    )

ax.text(
    E0,
    vmax,
    'Ne-like Xe; 3D',
    color = 'k'
    )
'''
#ax.set_xlim(xlim)
ax.set_xlim(hc/xlim[1], hc/xlim[0])
ax.grid('on')
ax.set_ylim(0,1.1*vmax)
#ax.legend(labelcolor='linecolor', loc = 'upper left')
#ax.set_xlabel('Photon energy [eV]')
leg = ax.legend(labelcolor='linecolor', loc = 'upper right')
leg.set_draggable('on')
ax.set_xlabel('Photon wavelength [AA]')
ax.set_ylabel('[counts/s]')
ax.yaxis.set_major_formatter(formatter)





specs, scram = rh.plot_histos(
    specs=specs,
    scram=scram,
    #x_lim = [4300,5000],
    x_line = [2179, 2459, 2595],
    #Te_init = 10.0,
    cs_init = 14,
    spec_init = 11,
    #plt_scram = True,
    scram_method = 'cs',
    num_ax = 4,
    )

specs, scram = rh.plot_histos(
    specs=specs,
    scram=scram,
    #x_units = 'AA',
    #x_lim = [2.7, 2.78],
    #x_line = [2.72, 2.729, 2.736, 2.747],
    #x_highlight = [[2.71, 2.755]],
    x_lim = [4450, 4600],
    x_line = [4560, 4543, 4531, 4513],
    x_highlight = [[4500,4575]],
    scram_method = 'cs',
    #cs_init = 14,
    cs_init = 23,
    spec_init = 11,
    num_ax = 4,
    )

specs, scram = rh.plot_histos(
    specs=specs,
    scram=scram,
    x_lim=[8000,9250],
    #x_lim = [7600,8200],
    scram_method='cs',
    cs_init = 14,
    spec_init = 13,
    num_ax = 4
    )


specs, scram = rh.plot_histos(
    specs=specs,
    scram=scram,
    x_lim=[1000,2000],
    #x_lim = [7600,8200],
    #x_line = [4560],
    scram_method='cs',
    cs_init = 14,
    spec_init = 17,
    num_ax = 4
    )

# XRS-LR
specs, scram = rh.plot_histos(
    specs=specs,
    scram=scram,
    x_lim=[6850, 8730],
    #x_lim = [7600,8200],
    #x_line = [4560],
    scram_method='cs',
    cs_init = 14,
    spec_init = 14,
    num_ax = 4
    )



################################################
#
#           Plot side-by-side
#
################################################

# Modules
from atomic_world.read_LLNL_EBIT.main import read_histos as rh

# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2

# Feature 1
specs, scram = rh.plot_side_by_side(
    specs=specs,
    scram=scram,
    plt_keys = ['8.71 kV', '8.97 kV', '9.23 kV'],
    x_lim=[[2740, 2850], [4500, 4580]],
    scram_method='cs',
    cs_init = 14,
    x_line = [4560, 4543, 4531, 4513],
    plt_scram = [True, False],
    )

# Feature 2
specs, scram = rh.plot_side_by_side(
    specs=specs,
    scram=scram,
    plt_keys = ['8.97 kV', '9.23 kV', '9.49 kV'],
    x_lim=[[2740, 2850], [4500, 4580]],
    scram_method='cs',
    cs_init = 14,
    x_line = [4560, 4543, 4531, 4513],
    plt_scram = [True, False],
    )

# Feature 3
specs, scram = rh.plot_side_by_side(
    specs=specs,
    scram=scram,
    plt_keys = ['9.23 kV', '9.49 kV', '9.75 kV'],
    x_lim=[[2740, 2850], [4500, 4580]],
    scram_method='cs',
    cs_init = 14,
    x_line = [4560, 4543, 4531, 4513],
    plt_scram = [True, False],
    )

# Feature 4
specs, scram = rh.plot_side_by_side(
    specs=specs,
    scram=scram,
    plt_keys = ['8.71 kV', '8.97 kV', '9.23 kV'],
    x_lim=[[2450, 2550], [4500, 4580]],
    scram_method='cs',
    cs_init = 14,
    x_line = [4560, 4543, 4531, 4513],
    plt_scram = [True, False],
    )

# Feature 5
specs, scram = rh.plot_side_by_side(
    specs=specs,
    scram=scram,
    plt_keys = ['9.23 kV', '9.49 kV', '9.75 kV'],
    x_lim=[[2450, 2550], [4500, 4580]],
    scram_method='cs',
    cs_init = 14,
    x_line = [4560, 4543, 4531, 4513],
    plt_scram = [True, False],
    )

# Feature 6
specs, scram = rh.plot_side_by_side(
    specs=specs,
    scram=scram,
    plt_keys = ['8.71 kV', '8.97 kV', '9.23 kV'],
    #plt_keys = ['8.45 kV', '8.71 kV', '8.97 kV'],
    x_lim=[[2600, 2700], [4500, 4580]],
    scram_method='cs',
    cs_init = 14,
    x_line = [4560, 4543, 4531, 4513],
    plt_scram = [True, False],
    )

# Feature 7
specs, scram = rh.plot_side_by_side(
    specs=specs,
    scram=scram,
    plt_keys = ['9.23 kV', '9.49 kV', '9.75 kV'],
    x_lim=[[2600, 2700], [4500, 4580]],
    scram_method='cs',
    cs_init = 14,
    x_line = [4560, 4543, 4531, 4513],
    plt_scram = [True, False],
    )

# Feature 8
specs, scram = rh.plot_side_by_side(
    specs=specs,
    scram=scram,
    plt_keys = ['9.23 kV', '9.49 kV', '9.75 kV'],
    x_lim=[[8900, 9150], [4500, 4580]],
    scram_method='cs',
    cs_init = 14,
    x_line = [4560, 4543, 4531, 4513],
    plt_scram = [True, False],
    )

# Feature 9
specs, scram = rh.plot_side_by_side(
    specs=specs,
    scram=scram,
    plt_keys = ['9.75 kV', '15(lowCB) kV', '15(highCB) kV'],
    x_lim=[[8900, 9150], [4500, 4580]],
    scram_method='cs',
    cs_init = 14,
    x_line = [4560, 4543, 4531, 4513],
    plt_scram = [True, False],
    )


# Feature 10
specs, scram = rh.plot_side_by_side(
    specs=specs,
    scram=scram,
    plt_keys = ['8.71 kV', '8.97 kV', '9.23 kV'],
    x_lim=[[3250, 3350], [4500, 4580]],
    scram_method='cs',
    cs_init = 14,
    x_line = [4560, 4543, 4531, 4513],
    plt_scram = [True, False],
    )

# Feature 11
specs, scram = rh.plot_side_by_side(
    specs=specs,
    scram=scram,
    plt_keys = ['9.23 kV', '9.49 kV', '9.75 kV'],
    x_lim=[[3250, 3350], [4500, 4580]],
    scram_method='cs',
    cs_init = 14,
    x_line = [4560, 4543, 4531, 4513],
    plt_scram = [True, False],
    )

# Feature 12
specs, scram = rh.plot_side_by_side(
    specs=specs,
    scram=scram,
    plt_keys = ['9.75 kV', '15(lowCB) kV', '15(highCB) kV'],
    x_lim=[[3250, 3350], [4500, 4580]],
    scram_method='cs',
    cs_init = 14,
    x_line = [4560, 4543, 4531, 4513],
    plt_scram = [True, False],
    )

# Feature 13
specs, scram = rh.plot_side_by_side(
    specs=specs,
    scram=scram,
    #plt_keys = ['8.71 kV', '8.97 kV', '9.23 kV'],
    #plt_keys = ['9.23 kV', '9.49 kV', '9.75 kV'],
    plt_keys = ['9.75 kV', '15(lowCB) kV', '15(highCB) kV'],
    x_lim=[[4100, 4200], [4500, 4580]],
    scram_method='cs',
    cs_init = 14,
    x_line = [4560, 4543, 4531, 4513],
    plt_scram = [True, False],
    )

# Feature 14
specs, scram = rh.plot_side_by_side(
    specs=specs,
    scram=scram,
    #plt_keys = ['8.71 kV', '8.97 kV', '9.23 kV'],
    #plt_keys = ['9.23 kV', '9.49 kV', '9.75 kV'],
    plt_keys = ['9.75 kV', '15(lowCB) kV', '15(highCB) kV'],
    x_lim=[[4550, 4650], [4500, 4580]],
    scram_method='cs',
    cs_init = 14,
    x_line = [4560, 4543, 4531, 4513],
    plt_scram = [True, False],
    )

# Feature 15
specs, scram = rh.plot_side_by_side(
    specs=specs,
    scram=scram,
    #plt_keys = ['8.71 kV', '8.97 kV', '9.23 kV'],
    #plt_keys = ['9.23 kV', '9.49 kV', '9.75 kV'],
    plt_keys = ['9.75 kV', '15(lowCB) kV', '15(highCB) kV'],
    x_lim=[[4450, 4550], [4500, 4580]],
    scram_method='cs',
    cs_init = 14,
    x_line = [4560, 4543, 4531, 4513],
    plt_scram = [True, False],
    )