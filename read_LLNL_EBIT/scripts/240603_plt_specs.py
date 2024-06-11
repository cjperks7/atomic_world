'''

Script to plot EBIT spectra

cjperks
June 3rd, 2024

'''

# Modules
from atomic_world.read_LLNL_EBIT.main import load_events as le
import os

# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2

# Files
file_path = os.path.join(
    '/home/cjperks/work',
    '2008_SPARCxray/LLNL_EBIT',
    '240531_W_data'
    )

'''
cases = [
    '240520_4', # 15kV
    '240521_1', # 15kV
    '240521_2', # 15kV, Be window
    '240521_3', # 15kV, Be window
    '240521_5', # 7.15kV, Be window
    '240522b_1', # 15kV
    '240522b_2', # 15kV
    '240522b_3', # 7.15kV
    '240522b_4', # 7.28kV
    '240523b_1', # 15kV, Be window
    '240523b_2', # 7.28kV, Be window
    '240523b_3', # 7.41kV, Be window
    '240523b_4', # 4kV, Be window
    '240523b_5', # 4kV, Be window
    '240523b_6', # 4.04kV, Be window
    '240523b_7', # 4.04kV, Be window
    '240523b_9', # 4.04kV, Be window
    '240524a_1', # 15kV, Be window
    '240524b_0', # 7.54kV, Be window
    '240524b_1', # 7.67kV, Be window
    ]

'''

case = '240521_5'

# Loads events list
df, time = le.load(
    file_path=file_path,
    case = case,
    )

le.plt_all(
    df=df,
    title = case
    )

# Filter by UTime
df_ut = le.fil_utime(
    df=df,
    UTime = None,
    case = case
    )


le.plt_all(
    df=df_ut,
    title =case,
    )


'''
le.plt_ebit_phase(
    df = df_ut,
    )

'''

# Filter by EBIT Phase
df_eb = le.fil_ebph(
    df=df_ut,
    case=case,
    )


le.plt_ebit_phase(
    df = df_eb,
    )



# Makes spectrum histogram
spec = le.make_hist(
    df=df_eb,
    plot=True,
    title =case,
    )







################# v. Beam energy

cases = [
    '240521_5', # 7.15kV, Be window
    '240522b_4', # 7.28kV
    '240523b_3', # 7.41kV, Be window
    '240524b_0', # 7.54kV, Be window
    '240524b_1', # 7.67kV, Be window
    ]
lab = [
    '7.15kV',
    '7.28kV',
    '7.41kV',
    '7.54kV',
    '7.67kV'
    ]


spec_tot = {}



for ii, cc in enumerate(cases):
    spec_tot[cc], _ = le.main(
        file_path=file_path,
        case = cc
        )

import matplotlib.gridspec as gridspec

plt.rcParams.update({'font.size': 16})
#fig, ax = plt.subplots(len(cases),1)
fig = plt.figure(figsize = (16,8))
gs = gridspec.GridSpec(len(cases), 1, height_ratios= [1]*len(cases))

axs = []
axs.append(fig.add_subplot(gs[0,0]))

xlim = [4200,4800]
#xlim = [2000,3000]

colormap = plt.get_cmap('winter')
colors = colormap(np.linspace(0,1,len(cases)))
    
for ii, cc in enumerate(cases):
    #ind = np.where(
    #    (hc/spec_tot[cc]['binlow'] >= xlim[0])
    #    & (hc/spec_tot[cc]['binlow'] <= xlim[1])
    #    )[0]
    ind = np.where(
        (spec_tot[cc]['binlow'] >= xlim[0])
        & (spec_tot[cc]['binlow'] <= xlim[1])
        )[0]

    if ii >0:
        axs.append(fig.add_subplot(gs[ii,0], sharex=axs[0]))
   

    axs[ii].step(
        spec_tot[cc]['binlow'],
        #hc/spec_tot[cc]['binlow'],
        spec_tot[cc]['hist'],
        #(
        #    spec_tot[cc]['hist']
        #    /np.max(spec_tot[cc]['hist'][ind])
        #    ),
        label = lab[ii],
        color = colors[ii]
        )

    if ii < len(cases)-1:
        plt.setp(axs[ii].get_xticklabels(), visible=False)
        axs[ii].set_xlabel('')
        plt.subplots_adjust(hspace=0)

    axs[ii].grid('on')
    #ax.set_xlabel(r'Wavelength, $\lambda$ [$\AA]')
    axs[ii].set_ylabel('Counts')

    axs[ii].set_xlim(xlim[0], xlim[1])
    #ax.set_xlim(hc/xlim[1], hc/xlim[0])
    axs[ii].set_ylim(0,1.2*np.max(spec_tot[cc]['hist'][ind]))

    leg = axs[ii].legend()
    leg.set_draggable('on')

    axs[ii].spines['left'].set_color(colors[ii])
    axs[ii].yaxis.label.set_color(colors[ii])
    axs[ii].tick_params(axis='y', color=colors[ii])
    for ll in axs[ii].get_yticklabels():
        ll.set_color(colors[ii])

    aa = 0.3
    axs[ii].axvspan(
        hc/2.76,
        hc/2.71,
        color = 'm',
        alpha = aa
        )

    axs[ii].plot(
        [hc/2.72,hc/2.72],
        [0,1.2*np.max(spec_tot[cc]['hist'][ind])],
        'k--'
        )

    ax1 = axs[ii].twinx()

    ax1.plot(
        out['indspec']['E_eV']['data'],
        (
            out['indspec'][num]['j(ph/ion/s/AA)'][ion,:]
            /np.max(out['indspec'][num]['j(ph/ion/s/AA)'][ion,:])
            ),
        color = 'red'
        )
    ax1.set_ylim(0,1.5)

    ax1.set_ylabel('Ni-like')

    ax1.spines['right'].set_color('red')
    ax1.yaxis.label.set_color('red')
    ax1.tick_params(axis='y', color='red')
    for ll in ax1.get_yticklabels():
        ll.set_color('red')

    ax2 = axs[ii].twinx()
    ax2.plot(
        out['indspec']['E_eV']['data'],
        (
            out['indspec'][num2]['j(ph/ion/s/AA)'][ion2,:]
            /np.max(out['indspec'][num2]['j(ph/ion/s/AA)'][ion2,:])
            ),
        color = 'm'
        )
    ax2.spines['right'].set_position(('outward', 60))

    ax2.set_ylim(0,1.5)

    ax2.set_ylabel('Ti-like')

    ax2.spines['right'].set_color('m')
    ax2.yaxis.label.set_color('m')
    ax2.tick_params(axis='y', color='m')
    for ll in ax2.get_yticklabels():
        ll.set_color('m')

axs[-1].set_xlabel('Photon Energy [eV]')



num = 'spec_10'
ion = 11

num2 = 'spec_11'
ion2 = 6

fig1, ax1 = plt.subplots()

ax1.plot(
    out['indspec']['E_eV']['data'],
    out['indspec']['spec_11']['j(ph/ion/s/AA)'][6,:]
)
ax1.set_xlim(2000,3000)










############ SCRAM

from atomic_world.run_SCRAM.scripts.tests import _read_sHANSEN as rh

# Species of interest
sp = 'W'
Znuc = 74
fi = 1.5e-5

ne_cm3 = 1e14

# Common path to files
filep = os.path.join(
    '/home/cjperks/work',
    '2008_SPARCxray/SCRAM',
    'from_sHANSEN/240529'
    )
name = sp+'_h_scram_'

out = rh.main(
    filepath=filep,
    name=name,
    #quants = ['spec', 'table']
    )










'''
# Modules
from atomic_world.read_LLNL_EBIT.main import load_events as le
import os

# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2

# Files
file_path = os.path.join(
    '/home/cjperks/work',
    '2008_SPARCxray/LLNL_EBIT',
    '240531_W_data'
    )

# case = '240520_4'
case = '240521_3' # 15kV, Be window

spec, _ = le.main(
    file_path=file_path,
    case = case
    )

plt.rcParams.update({'font.size': 16})
fig, ax = plt.subplots()
fig.set_size_inches(16,7)

ax.step(
    spec['binlow'],
    spec['hist'],
    color = 'blue'
    )

ax.grid('on')
ax.set_xlabel('Photon Energy [eV]')
ax.set_ylabel('Counts')

ax.set_xlim(1e3,1e4)
ax.set_ylim(0,350)

ax.set_title('240521; E_beam =15 kV')

ax.spines['left'].set_color('blue')
ax.yaxis.label.set_color('blue')
ax.tick_params(axis='y', color='blue')
for ll in ax.get_yticklabels():
    ll.set_color('blue')

import scipy.constants as cnt
hc = cnt.h*cnt.c/cnt.e*1e10 # [eV*AA]

aa = 0.3
ax.axvspan(
    hc/2.76,
    hc/2.71,
    color = 'm',
    alpha = aa
    )
ax.plot(
    [hc/2.72, hc/2.72],
    [0,350],
    color = 'red'
    )

ax.text(
    2000,
    325,
    'M-shell (n=4->3)',
    color = 'red'
    )
ax.text(
    8000,
    #100,
    200,
    'L-shell (n=3->2)',
    color = 'red'
    )

num = 'spec_26'
ind = 7

ax1 = ax.twinx()

ax1.plot(
    out['indspec']['E_eV']['data'],
    (
        out['indspec'][num]['j(ph/ion/s/AA)'][ind,:]
        /np.max(out['indspec'][num]['j(ph/ion/s/AA)'][ind,:])
        ),
    color = 'green'
    )
ax1.set_ylim(0,5)

ax1.set_ylabel('SCRAM simulation (Ne-like W)')

ax1.spines['right'].set_color('green')
ax1.yaxis.label.set_color('green')
ax1.tick_params(axis='y', color='green')
for ll in ax1.get_yticklabels():
    ll.set_color('green')
'''


'''

ne_cm3 = 1e14
fi = 1.5e-5

case = '240520_4'
#case = '240521_3' # 15kV, Be window

spec, _ = le.main(
    file_path=file_path,
    case = case
    )


plt.rcParams.update({'font.size': 16})
fig, ax = plt.subplots()
fig.set_size_inches(16,7)

ax.plot(
    [hc/2.72,hc/2.72],
    [0,60],
    'm--',
    linewidth = 3
)

ax.step(
    spec['binlow'],
    spec['hist'],
    color = 'k'
    )

ax.grid('on')
ax.set_xlabel('Photon Energy [eV]')
ax.set_ylabel('Counts')

ax.set_xlim(1e3,1e4)
ax.set_ylim(0,350)

ax.set_title('240520; E_beam =15 kV')

#ax.spines['left'].set_color('blue')
#ax.yaxis.label.set_color('blue')
#ax.tick_params(axis='y', color='blue')
#for ll in ax.get_yticklabels():
#    ll.set_color('blue')

#aa = 0.3
#ax.axvspan(
#    hc/2.76,
#    hc/2.71,
#    color = 'm',
#    alpha = aa
#    )
#ax.plot(
#    [hc/2.72, hc/2.72],
#    [0,350],
#    color = 'red'
#    )
lw = 3
ax.plot(
    [hc/2.76, hc/2.76],
    [0,60],
    'b--',
    linewidth = lw
    )
ax.plot(
    [hc/2.71, hc/2.71],
    [0,60],
    'b--',
    linewidth = lw
    )


ax.set_xlim(4400, 4700)
ax.set_ylim(0,60)

ax1 = ax.twinx()

num = 'spec_19'
ions = [5,6,7,8]
lab = ['Na-like (11)', 'Mg-like (12)', 'Al-like (13)', 'Si-like (14)']

for ii, ion in enumerate(ions):
    ax1.plot(
    out['indspec']['E_eV']['data'],
    (
        out['indspec'][num]['j(ph/ion/s/AA)'][ion,:]
        * out['indspec'][num]['abund'][ion]
        /out['table']['ne_cm3']['data'][0]
        * ne_cm3
        * fi * ne_cm3
        ),
    label = lab[ii],
    linewidth = 2
    )

ax1.set_ylim(0,8e12)

ax1.set_ylabel('SCRAM simulation')

leg = ax1.legend()
leg.set_draggable('on')

ax1.spines['right'].set_color('green')
ax1.yaxis.label.set_color('green')
ax1.tick_params(axis='y', color='green')
for ll in ax1.get_yticklabels():
    ll.set_color('green')