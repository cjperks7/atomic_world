'''

Module to handle reading and plotting histograms from the
LLNL EBIT recieved from N. Hell on June 1st, 2024

cjperks
June 7th, 2024

'''

# Modules
import numpy as np
import matplotlib.pyplot as plt
import os, sys
from atomic_world.run_SCRAM.scripts.tests import _read_sHANSEN as rh
import matplotlib.gridspec as gridspec
import scipy.constants as cnt
from matplotlib.widgets import Slider

# Useful constants
hc = cnt.h*cnt.c/cnt.e*1e10 # [eV*AA]

# Plasma
ne_cm3 = 1e14
fi = 1.5e-5

__all__ = [
    'read_histos',
    'read_SCRAM',
    'plt_histos'
    ]


############################################
#
#           Main
#
############################################

# Read histogram data
def read_histos(
    file_path = None,
    histos = None,
    common = ['W_Escan_', '-keV.hist']
    ):

    # Defaults
    if file_path is None:
        file_path = os.path.join(
            '/home/cjperks/work',
            '2008_SPARCxray/LLNL_EBIT',
            '240601_W_histos'
            )

    if histos is None:
        histos = [
            '3.93',
            '4.04',
            '7.15',
            '7.28',
            '7.41',
            '7.54',
            '7.67(highCB)',
            '7.67(lowCB)',
            '7.93',
            '8.19',
            '8.27',
            '8.45',
            '8.71',
            '8.97',
            '15(highCB)',
            '15(lowCB)'
            ]

    # Init
    specs = {}

    # Loop over beam energies
    for histo in histos:
        if common is not None:
            if len(common) == 2:
                name = '/'+common[0]+histo+common[1]
                key = histo[0:] + ' kV'
        else:
            name = '/'+histo
            key = name[0:] + ' kV'

        specs[key] = _load_histo_ascii(
            file_path = file_path,
            name = name
            )

    # Output
    return specs

# Loads histogram ascii files
def _load_histo_ascii(
    file_path = None,
    name = None
    ):

    # Init
    out = {}
    out['binlow'] = []
    out['binhigh'] = []
    out['hist'] = []
    out['hist_all'] = []

    # Reads file
    with open(file_path+name, 'r') as ff:
        for line in ff:
            out['binlow'].append(
                float(line.split()[0])
                )
            out['binhigh'].append(
                float(line.split()[1])
                )
            out['hist_all'].append(
                float(line.split()[2])
                )
            out['hist'].append(
                float(line.split()[3])
                )

    ff.close()

    # Reformat
    out['binlow'] = np.asarray(out['binlow'])
    out['binhigh'] = np.asarray(out['binhigh'])
    out['hist_all'] = np.asarray(out['hist_all'])
    out['hist'] = np.asarray(out['hist'])

    # Output
    return out


############################################
#
#           SCRAM data
#
############################################

def read_SCRAM(
    file_path = None,
    name = 'W_h_scram_'
    ):

    # Default
    if file_path is None:
        file_path = os.path.join(
            '/home/cjperks/work',
            '2008_SPARCxray/SCRAM',
            'from_sHANSEN/240529'
            )

    # Reads SCRAM data files
    out = rh.main(
        filepath=file_path,
        name=name,
        #quants = ['spec', 'table']
        )

    # Output
    return out

############################################
#
#           Plotting
#
############################################

def _cs_label(
    cs = None
    ):

    labs = {
        '01': 'H-like',
        '02': 'He-like',
        '03': 'Li-like',
        '04': 'Be-like',
        '05': 'B-like',
        '06': 'C-like',
        '07': 'N-like',
        '08': 'O-like',
        '09': 'F-like',
        '10': 'Ne-like',
        '11': 'Na-like',
        '12': 'Mg-like',
        '13': 'Al-like',
        '14': 'Si-like',
        '15': 'P-like',
        '16': 'S-like',
        '17': 'Cl-like',
        '18': 'Ar-like',
        '19': 'K-like',
        '20': 'Ca-like',
        '21': 'Sc-like',
        '22': 'Ti-like',
        '23': 'V-like',
        '24': 'Cr-like',
        '25': 'Mn-like',
        '26': 'Fe-like',
        '27': 'Co-like',
        '28': 'Ni-like',
        '29': 'Cu-like',
        '30': 'Zn-like',
        '31': 'Ga-like',
        '32': 'Ge-like',
        '33': 'As-like',
        '34': 'Se-like',
        '35': 'Br-like',
        '36': 'Kr-like'
        }

    # Output
    try:
        return labs[cs]
    except:
        return 'missing '+cs

# Updating spectrum plotting
def _update_specs(
    # Figure
    fig = None,
    axs_l = None,
    axs_r = None,
    spec_ind = None,
    # ECS Data
    specs = None,
    keys = None,
    # SCRAM data
    scram = None,
    # X labels
    x_units = None,
    x_lim = None,
    # Color controls
    colors = None,
    # Extra
    x_line = None,
    dE = None,
    plt_scram = None,
    # Slider
    scale1 = None,
    scale2 = None, # [keV]
    ):

    # Finds spectra to plot
    kks = [keys[ii] for ii in spec_ind+int(scale1)]

    # Finds SCRAM spectra to plot
    scr_ind = np.argmin(
        abs(scram['table']['Te_eV']['data']/1e3 - scale2)
        )
    scr_kk = 'spec_%i'%(scr_ind)

    # Plot over axes
    for ii in np.arange(len(axs_l)):
        axs_l[ii].clear()
        # Prepares x-data
        if x_units == 'eV':
            bins = np.where(
                (specs[kks[ii]]['binlow'] >= x_lim[0])
                & (specs[kks[ii]]['binlow'] <= x_lim[1])
                )[0]
            xbin = specs[kks[ii]]['binlow']

            xscram = scram['indspec']['E_eV']['data'] -dE
            tscram = np.where(
                (scram['indspec']['E_eV']['data'] >= x_lim[0])
                & (scram['indspec']['E_eV']['data'] <= x_lim[1])
                )[0]
        elif x_units == 'AA':
            bins = np.where(
                (specs[kks[ii]]['binlow'] >= hc/x_lim[1])
                & (specs[kks[ii]]['binlow'] <= hc/x_lim[0])
                )[0]
            xbin = hc/specs[kks[ii]]['binlow']

            xscram = hc/(scram['indspec']['E_eV']['data']-dE)
            tscram = np.where(
                (scram['indspec']['lambda_A']['data'] >= hc/x_lim[1])
                & (scram['indspec']['lambda_A']['data'] <= hc/x_lim[0])
                )[0]

        # Plots
        axs_l[ii].step(
            xbin,
            specs[kks[ii]]['hist'],
            color = colors[ii]
            )

        # Removes x-axis
        if ii < len(axs_l)-1:
            plt.setp(axs_l[ii].get_xticklabels(), visible=False)
            axs_l[ii].set_xlabel('')
            plt.subplots_adjust(hspace=0)

        # Prepares y-axis
        axs_l[ii].grid('on')

        axs_l[ii].spines['left'].set_color(colors[ii])
        axs_l[ii].yaxis.label.set_color(colors[ii])
        axs_l[ii].tick_params(axis='y', color=colors[ii])
        for ll in axs_l[ii].get_yticklabels():
            ll.set_color(colors[ii])

        # Axis limits
        axs_l[ii].set_xlim(x_lim[0], x_lim[1])
        axs_l[ii].set_ylim(0,1.2*np.max(specs[kks[ii]]['hist'][bins]))

        # Y-label
        if len(kks[ii].split('('))>1:
            axs_l[ii].set_ylabel(
                kks[ii].split('(')[0]
                +' kV\n('
                +kks[ii].split('(')[1][:-3]
                )

        else:
            axs_l[ii].set_ylabel(kks[ii])

        ### --- SCRAM data --- ###
        # Init axis
        axs_r[ii].clear()

        scr_data = (
            scram['indspec'][scr_kk]['j(ph/ion/s/AA)']
            * scram['indspec'][scr_kk]['abund'][:,None]
            /scram['table']['ne_cm3']['data'][0]
            * ne_cm3
            * fi * ne_cm3
            ) # dim(ncs, nlambda)

        mm1 = np.max(scr_data[:,tscram], axis = 1)
        mm2 = np.argsort(mm1)[-10:][::-1]

        scr_lim = np.max(scr_data[:,tscram].flatten())

        # Loop over charge states
        for jj in np.arange(scr_data.shape[0]):
            #if scram['indspec'][scr_kk]['ions'][jj].split('-')[0] == 'd':
            if jj in mm2 and plt_scram:
                axs_r[ii].plot(
                    xscram,
                    scr_data[jj,:],
                    label = _cs_label(
                        cs = scram['indspec'][scr_kk]['ions'][jj].split('-')[1]
                        )
                    )

        # Color axis
        axs_r[ii].spines['right'].set_color('red')
        axs_r[ii].set_ylabel('SCRAM')
        axs_r[ii].yaxis.label.set_color('red')
        axs_r[ii].tick_params(axis='y', color='red')
        for ll in axs_r[ii].get_yticklabels():
            ll.set_color('red')
        axs_r[ii].set_ylim(0, 1.2*scr_lim)


        ### --- Extras --- ###
        if x_line is not None:
            for xx in x_line:
                axs_l[ii].plot(
                    [xx,xx],
                    [0, 1.2*np.max(specs[kks[ii]]['hist'][bins])],
                    'k--'
                    )

    if plt_scram:
        axs_r[0].legend(
            loc = 'upper left',
            bbox_to_anchor = (1,0.5),
            labelcolor='linecolor'
            )

    fig.canvas.draw_idle()  # Redraw the plot

# Imports data to update plots
def _update_func(
    fig = None,
    axs_l = None,
    axs_r = None,
    spec_ind = None,
    specs = None,
    keys = None,
    scram = None,
    x_units = None,
    x_lim = None,
    colors = None,
    x_line = None,
    dE = None,
    plt_scram = None,
    ):

    def update(scale1,scale2):
        _update_specs(
            fig = fig,
            axs_l = axs_l,
            axs_r = axs_r,
            spec_ind = spec_ind,
            specs = specs,
            keys = keys,
            scram = scram,
            x_units = x_units,
            x_lim = x_lim,
            colors = colors,
            x_line = x_line,
            dE = dE,
            plt_scram = plt_scram,
            scale1 = scale1,
            scale2 = scale2, # [eV]
            )
    return update

# Plots histograms
def plot_histos(
    specs = None,
    scram = None,
    # Figure controls
    font_size = 16,
    fig_size = (16,8),
    title = 'LLNL EBIT ECS spectra (N.Hell, 06/01/24)',
    # Axis controls
    num_ax = 5,
    x_units = 'eV',
    x_lim = [2000,3000],
    # Color contols
    spec_map = 'winter',
    # Slider controls
    Te_init = 5, # [keV]
    spec_init = 2,
    # Extra
    x_line = None,
    dE = 0.0,
    plt_scram = True,
    ):

    # Default values
    if specs is None:
        specs = read_histos()

    # Loads SCRAM data files
    if scram is None:
        scram = read_SCRAM()

    # Loaded spectra
    keys = list(specs.keys())

    # Init figure
    plt.rcParams.update({'font.size': font_size})
    fig = plt.figure(figsize = fig_size)
    gs = gridspec.GridSpec(
        num_ax, 
        1, 
        height_ratios= [1]*num_ax
        )

    # Init colors
    colormap = plt.get_cmap(spec_map)
    colors = colormap(np.linspace(0,1,num_ax))

    # Init axis labels
    if x_units == 'eV':
        xlabel = 'Photon energy [eV]'
    elif x_units == 'AA':
        xlabel = r'Photon wavelength, $\lambda$ [$\AA$]'

    # Init axes
    axs_l = []
    axs_r = []
    for ii in np.arange(num_ax):
        if ii >0:
            axs_l.append(fig.add_subplot(gs[ii,0], sharex=axs_l[0]))
        else:
            axs_l.append(fig.add_subplot(gs[0,0]))

        axs_r.append(axs_l[-1].twinx())

    axs_l[-1].set_xlabel(xlabel)
    fig.suptitle(title)

    # Init slider indexing
    spec_ind = np.asarray(np.arange(
        num_ax,
        ))
    spec_ind_1 = int(np.floor(num_ax/2))
    spec_ind_2 = len(keys) - int(np.ceil(num_ax/2))

    # Init update data wrapper
    interact_func = _update_func( 
        fig = fig,
        axs_l = axs_l,
        axs_r = axs_r,
        spec_ind = spec_ind-spec_ind_1,
        specs = specs,
        keys = keys,
        scram = scram,
        x_units = x_units,
        x_lim = x_lim,
        colors = colors,
        x_line = x_line,
        dE = dE,
        plt_scram = plt_scram,
        )
    interact_func(spec_init, Te_init)
    
    # Add sliders below the subplots
    ax_slider1 = plt.axes(
        [0.2, 0.01, 0.65, 0.03], 
        facecolor='lightgoldenrodyellow'
        )
    ax_slider2 = plt.axes(
        [0.2, 0.06, 0.65, 0.03], 
        facecolor='lightgoldenrodyellow'
        )

    # Sliders
    slider1 = Slider(ax_slider1, 'Spectra',
        spec_ind_1, spec_ind_2,
        valinit=spec_init, valstep=1
        )
    slider2 = Slider(ax_slider2, 'Te [keV]',
        np.min(scram['table']['Te_eV']['data']/1e3), 
        np.max(scram['table']['Te_eV']['data']/1e3),
        valinit=Te_init,
        )

    # Function to be called anytime a slider's value changes
    def sliders_on_changed(val):
        interact_func(slider1.val, slider2.val)

    # Attach the update function to sliders
    slider1.on_changed(sliders_on_changed)
    slider2.on_changed(sliders_on_changed)

    # Output
    return specs, scram
