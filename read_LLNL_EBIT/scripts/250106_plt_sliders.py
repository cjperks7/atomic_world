'''

Script to generate slider plots of LLNL EBIT 
microcalorimeter spectra

cjperks
Jan 6th, 2025

'''

# Modules
from atomic_world.read_LLNL_EBIT.main import read_histos as rh
from atomic_world.run_SCRAM.scripts.tests import _read_sHANSEN82 as rh82
import sys,os

# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2

# Make histogram plots
specs, scram = rh.plot_histos(
    x_line = [2179],
    #dE = 5, # [eV]
    scram_method = 'Gauss',
    num_ax = 4,
    spec_init = 15,
    )

# Make histogram plots
specs, scram = rh.plot_histos(
    specs = specs,
    scram = scram,
    x_line = [2179],
    #dE = 5, # [eV]
    x_lim = [4450, 4650],
    scram_method = 'Gauss',
    num_ax = 4,
    spec_init = 15,
    Eb_init = 9.75
    )

# Make histogram plots
specs, scram = rh.plot_histos(
    specs = specs,
    scram = scram,
    x_line = [2179],
    #dE = 5, # [eV]
    x_lim = [2740, 2850],
    scram_method = 'Gauss',
    num_ax = 4,
    spec_init = 15,
    Eb_init = 9.75
    )



################################################
#
#           Plot side-by-side
#
################################################

# Feature 9
specs, scram = rh.plot_side_by_side(
    specs=specs,
    scram=scram,
    plt_keys = ['9.75 kV', '15(lowCB) kV', '15(highCB) kV'],
    #x_lim=[[8900, 9150], [4500, 4580]],
    x_lim=[[7700, 8350], [4400, 4700]],
    scram_method='Gauss',
    Eb_init = 9.49,
    x_line = [4560, 4543, 4531, 4513],
    #plt_scram = [True, False],
    plt_scram = [True, True]
    )

################################################
#
#           Line list
#
################################################

from atomic_world.run_SCRAM.scripts.tests import _read_sHANSEN82 as rh82
import sys,os

# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2

com = os.path.join(
    '/home/cjperks/work',
    '2008_SPARCxray/SCRAM',
    'from_sHANSEN/241220'
    )
name = 'W_ebit_8.71.linelist'
#name = 'XeL_9.23.linelist'

dline = rh82.read_linelist82(
    file = os.path.join(com, name),
    x_units = 'eV',
    x_lim = [7800, 8400]
    )


default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'H']
lines = ['-', '--', ':', '-.', '-', '--', ':', '-.', '-', '--', ':', '-.']

fig, ax = plt.subplots()

for kk in dline.keys():
    cs = rh82._get_nele(sym=kk)
    cs_ind = (cs-10)%10

    ax.plot(
        dline[kk]['dE']['data'],
        dline[kk]['ems']['data'],
        markers[cs_ind],
        color = default_colors[cs_ind],
        label = kk
        )

ax.grid('on')
leg = ax.legend(labelcolor='linecolor')
leg.set_draggable('on')



### ---- Finds significant lines ---- ####

vmax = 1
for kk in dline.keys():
    vmax = max(np.r_[vmax, max(dline[kk]['ems']['data'])])

dfilt = {}
tol = 0.01
for kk in dline.keys():
    ind = np.where(dline[kk]['ems']['data'] > tol*vmax)
    dfilt[kk] = {}

    for dd in dline[kk].keys():
        dfilt[kk][dd] = dline[kk][dd]['data'][ind]






####### --- L-shell Xe ---- #########


file_path = os.path.join(
    '/home/cjperks/work',
    '2008_SPARCxray/SCRAM',
    'from_sHANSEN/241220'
    )


dXe = rh82.main82(
    filepath = file_path,
    name = 'XeL_',
    quants = ['spec', 'indspec']
    )



###### ---- Plots of whole spectra ------- #########

Eb = '9.23'
key = Eb+' kV'
sp = 'spec_%skV'%(Eb)
fig, ax = plt.subplots()

ax.step(
    specs[key]['binlow'],
    specs[key]['hist'],
    where='post',
    color = 'k'
    )

ax2 = ax.twinx()

ax2.plot(
    scram['spec'][sp]['E_eV']['data'],
    scram['spec'][sp]['j(ph/s/cm^3/eV)'],
    #scram['spec'][sp]['j(W/cm^3/eV)'],
    color = 'b'
    )

ax2.plot(
    dXe['spec'][sp]['E_eV']['data'],
    dXe['spec'][sp]['j(ph/s/cm^3/eV)'],
    #scram['spec'][sp]['j(W/cm^3/eV)'],
    color = 'r'
    )

ax.set_xlim(1e3, 1e4)
ax.grid('on')