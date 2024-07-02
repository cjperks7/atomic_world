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
from Pinhole_Optimization import _get_transmission

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

# Be window thickness
def _get_Be(key=None):

    mil2meter = 25.4e-6
    dwindows = {
        '3.93 kV': 0.5*mil2meter, # [m]
        '4.04 kV': 5.5*mil2meter, # [m]
        '7.15 kV': 0.5*mil2meter, # [m]
        '7.28 kV': 0.5*mil2meter, # [m]
        '7.41 kV': 5.5*mil2meter, # [m]
        '7.54 kV': 5.5*mil2meter, # [m]
        '7.67(highCB) kV': 5.5*mil2meter, # [m]
        '7.67(lowCB) kV': 5.5*mil2meter, # [m]
        '7.93 kV': 5.5*mil2meter, # [m]
        '8.19 kV': 5.5*mil2meter, # [m]
        '8.27 kV': 5.5*mil2meter, # [m]
        '8.45 kV': 5.5*mil2meter, # [m]
        '8.71 kV': 5.5*mil2meter, # [m]
        '8.97 kV': 5.5*mil2meter, # [m]
        '15(highCB) kV': 5.5*mil2meter, # [m]
        '15(lowCB) kV': 0, # [m]
        }

    # Output
    return dwindows[key]

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
        return labs[cs.rjust(2,'0')] +' ('+ cs+')'
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
    x_highlight = None,
    dE = None,
    plt_scram = None,
    scram_method = None,
    css = None,
    # Slider
    scale1 = None,
    scale2 = None, # [keV]
    ):

    # Finds spectra to plot
    kks = [keys[ii] for ii in spec_ind+int(scale1)]

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
                (scram['indspec']['lambda_A'] >= x_lim[0])
                & (scram['indspec']['lambda_A'] <= x_lim[1])
                )[0]

        # Get Be window filtering
        fil_Be = _get_transmission(
            key = 'Be',
            E = scram['indspec']['E_eV']['data'], # [eV]
            thick = _get_Be(key=kks[ii]), # [m]
            ) # dim(nlambda,)

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

        # Just plot the Maxwellian envolope
        if scram_method == 'Max':
            # Finds SCRAM spectra to plot
            scr_ind = np.argmin(
                abs(scram['table']['Te_eV']['data']/1e3 - scale2)
                )
            scr_kk = 'spec_%i'%(scr_ind)

            scr_data = (
                scram['indspec'][scr_kk]['j(ph/ion/s/AA)']
                * scram['indspec'][scr_kk]['abund'][:,None]
                /scram['table']['ne_cm3']['data'][0]
                * ne_cm3
                * fi * ne_cm3
                ) # dim(ncs, nlambda)
            scr_data *= fil_Be[None,:]

            mm1 = np.max(scr_data[:,tscram], axis = 1)
            mm2 = np.argsort(mm1)[-10:][::-1]

            scr_lim = 1

            # Loop over charge states
            for jj in np.arange(scr_data.shape[0]):
                if scram['indspec'][scr_kk]['ions'][jj].split('-')[0] == 's':
                    continue
                if jj in mm2 and plt_scram:
                    scr_lim = np.max((scr_lim, np.max(scr_data[jj,tscram])))
                    axs_r[ii].plot(
                        xscram,
                        scr_data[jj,:],
                        label = _cs_label(
                            cs = scram['indspec'][scr_kk]['ions'][jj].split('-')[1]
                            )
                        )
        
        # Just plot the most abundant charge states
        elif scram_method == 'cs':
            plt_cs = list(range(scale2-4, scale2+6))

            for cs in plt_cs:
                if cs not in css.keys():
                    continue

                scr_data = (
                    scram['indspec']['spec_%i'%(css[cs]['spec'])]['j(ph/ion/s/AA)'][css[cs]['ind'],:]
                    ).copy()
                #scr_data *= fil_Be
                scr_data /= np.max(scr_data[tscram])
                scr_lim = 2

                axs_r[ii].plot(
                    xscram,
                    scr_data,
                    label = _cs_label(cs=str(cs))
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

        if x_highlight is not None:
            for xxs in x_highlight:
                axs_l[ii].axvspan(
                    xxs[0],
                    xxs[1],
                    color = 'm',
                    alpha = 0.3
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
    x_highlight = None,
    dE = None,
    plt_scram = None,
    scram_method = None,
    css = None,
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
            x_highlight = x_highlight,
            dE = dE,
            plt_scram = plt_scram,
            scram_method = scram_method,
            css = css,
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
    scram_method = 'Max', # 'Max or 'cs'
    cs_init = 28,
    Te_init = 5, # [keV]
    spec_init = 2,
    # Extra
    x_line = None,
    x_highlight = None,
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

    # Finds SCRAM run where each charge state was maximized
    if scram_method == 'cs':
        css = {}
        for cs in range(32,5,-1):
            css[cs] = {
                'abund': 0,
                'spec': -1,
                'ind': -1,
                }
            for tt in range(scram['table']['Te_eV']['data'].shape[0]):
                if 'd-%02i'%(cs) in scram['indspec']['spec_%i'%(tt)]['ions']:
                    ind = scram['indspec']['spec_%i'%(tt)]['ions'].index('d-%02i'%(cs))

                    if scram['indspec']['spec_%i'%(tt)]['abund'][ind] > css[cs]['abund']:
                        css[cs]['abund'] = scram['indspec']['spec_%i'%(tt)]['abund'][ind]
                        css[cs]['spec'] = tt
                        css[cs]['ind'] = ind
    elif scram_method == 'Max':
        css = None

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
        x_highlight = x_highlight,
        dE = dE,
        plt_scram = plt_scram,
        scram_method = scram_method,
        css = css,
        )
    if scram_method == 'Max':
        interact_func(spec_init, Te_init)
    elif scram_method == 'cs':
        interact_func(spec_init, cs_init)
    
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

    if scram_method == 'Max':
        slider2 = Slider(ax_slider2, 'Te [keV]',
            np.min(scram['table']['Te_eV']['data']/1e3), 
            np.max(scram['table']['Te_eV']['data']/1e3),
            valinit=Te_init,
            )
        interact_func(spec_init, Te_init)
    elif scram_method == 'cs':
        slider2 = Slider(ax_slider2, 'num. elec.',
            #32, 6,
            6, 32,
            valinit=cs_init,
            valstep=1
            )
        interact_func(spec_init, cs_init)

    # Function to be called anytime a slider's value changes
    def sliders_on_changed(val):
        interact_func(slider1.val, slider2.val)

    # Attach the update function to sliders
    slider1.on_changed(sliders_on_changed)
    slider2.on_changed(sliders_on_changed)

    # Output
    return specs, scram
