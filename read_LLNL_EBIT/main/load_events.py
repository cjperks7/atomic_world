'''

Script to prepare Event list from
LLNL's EBIT ECS (microcalorimeter)

cjperks
June 3rd, 2024

'''


# Modules
import numpy as np
import matplotlib.pyplot as plt
import os
import ebit_ecs as eb

__all__ = [
    'load'
    ]

#########################################
#
#           Main
#
#########################################

#Preloads main workflow
def main(
    file_path=None,
    case = None,
    ):

    # Loads events list
    df, time = load(
        file_path=file_path,
        case = case,
        )

    # Filter by UTime
    df_ut = fil_utime(
        df=df,
        case = case
        )

    # Filter by EBIT Phase
    df_eb = fil_ebph(
        df=df_ut,
        case=case,
        )

    # Makes spectrum histogram
    spec = make_hist(
        df=df_eb,
        )

    # Output
    return spec, df_eb

# Loads, calibrates, and filters event list
def load(
    file_path=None,
    events = None,
    calib = None,
    case = None,
    flag_list = [64], # list of event quality
    ):

    if case is not None:
        tmp = _def_240531(case=case)
        events = tmp['events']
        calib = tmp['calib']

    # Loads raw events list
    df, time = eb.ecs_import(
        os.path.join(
            file_path,
            events
            )
        )
    
    # Calibrates with Gain file
    df_cal = eb.ecs_apply_gain(
        os.path.join(
            file_path,
            calib
            ),
        df
        )
    
    # List of lower energy pixels
    lp = eb.low_energy_pix()

    # Filters by pixel and event quality
    df_cal_pix = df_cal.loc[
        df_cal['Flags'].isin(flag_list) 
        & df_cal['Pixel'].isin(lp)
        ]

    # Ouput events
    return df_cal_pix, time

# Filter by UTime
def fil_utime(
    df = None,
    UTime = None, # [ks]
    case = None,
    ):

    # Loads default UTime
    if UTime is None:
        UTime = _def_240531(case=case)['UTime']

    df_ut = df.loc[
        (df['UTime'] >= UTime[0]) 
        & (df['UTime'] <= UTime[1])
        ]

    # Ouput
    return df_ut

# Filter by EBIT phase
def fil_ebph(
    df = None,
    ebph = None,
    case = None,
    ):

    # Loads default EBIT Phase
    if ebph is None:
        ebph = _def_240531(case=case)['EBIT_Phase']

    df_eb = df.loc[
        (df['EbitPhase'] >= ebph[0]) 
        & (df['EbitPhase'] <= ebph[1])
        ]

    # Ouput
    return df_eb

# Histograms the event list
def make_hist(
    df = None,   # --> assumes filtered by UTime
    plot = False,
    title = None,
    ):

    # Makes histogram
    spec = eb.histo(df['PH'], 1,-100,10000)

    # Plots
    if plot:
        plt.step(spec['binlow'],spec['hist'])
        plt.xlabel('Photon Energy [eV]')
        plt.ylabel('Counts [arb]')
        plt.grid('on')

        if title is not None:
            plt.title(title)

    # Ouput
    return spec
    

#########################################
#
#           Plots
#
#########################################

# Plots all pulses
def plt_all(
    df=None,
    s = 0.1,
    title =None,
    ):

    plt.scatter(
        df['UTime'], 
        df['PH'],
        s=s
        )

    plt.xlabel(
        'Utime [ks]'
        )
    plt.ylabel(
        'Photon Energy [eV]'
        )

    plt.ylim(0,17500)

    if title is not None:
        plt.title(title)
    
# Plots EBIT Phase
def plt_ebit_phase(
    df = None,          # --> assumes filtered by UTime
    # Binning controls
    utime_stepsize = 1, # [ks]
    ebitphase_stepsize = 0.1,
    ph_range = np.arange(0, 10001, 1), # [eV]
    # Plot controls
    plt_utime  = False,
    plt_ebitphase = True,
    vmin = 0,
    vmax = 1,
    ):

    # UTime binning
    utime_range = np.arange(
        df['UTime'].min(), 
        df['UTime'].max() + utime_stepsize, 
        utime_stepsize
        )
    
    # EBIT Phase binning
    ebitphase_range = np.arange(
        df['EbitPhase'].min(), 
        df['EbitPhase'].max() + ebitphase_stepsize, 
        ebitphase_stepsize
        )
    
    # Makes UTime histogram
    hist, x_edges, y_edges = np.histogram2d(
        df['UTime'], 
        df['PH'], 
        bins=[utime_range, ph_range]
        )
    
    # Makes EBIT phase histogram
    hist_ebitphase, x_edges_ebitphase, y_edges_ebitphase = np.histogram2d(
        df['PH'], 
        df['EbitPhase'],
        bins=[ph_range, ebitphase_range]
        )

    # Plot pulse height v. UTime
    if plt_utime:
        plt.imshow(
            hist.T, 
            aspect='auto', 
            cmap='plasma', 
            origin='lower', 
            extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]], 
            vmin=vmin, vmax=vmax
            )
        plt.xlabel('UTime [ks]')
        plt.ylabel('Photon Energy [eV]')
        plt.colorbar(label='Counts [arb]')
        plt.title('Intensity Plot - Pulse Height')
        plt.show()

    # Plot pusle height v. EBIT Phase
    if plt_ebitphase:
        # Plot the 2D histogram with adjusted colorbar range
        plt.imshow(
            hist_ebitphase.T, 
            aspect='auto', 
            cmap='plasma', 
            origin='lower',
            extent=[x_edges_ebitphase[0], x_edges_ebitphase[-1], y_edges_ebitphase[0], y_edges_ebitphase[-1]], 
            vmin=vmin, 
            vmax=vmax
            )
        plt.xlabel('Photon Energy [eV]')
        plt.ylabel('EBIT Phase [s]')
        plt.colorbar(label='Counts [arb]')
        plt.title('Intensity Plot - EBIT Phase')
        plt.show()



#########################################
#
#           Defaults
#
#########################################

# Default settings for exp. over 240520-240531
def _def_240531(
    case = None,
    ):

    cases = {

        #### --- May 20th, 2024 --- ####
        '240520_0': {
            'UTime': [0,5.448], # [ks]
            'description': 'background',
            'beam_energy': '17kV' ,
            'events': 'Events.20240520T183043.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0,10]
            },
        '240520_1': {
            'UTime': [5.56,6.52], # [ks]
            'description': 'start W',
            'beam_energy': '17kV',
            'events': 'Events.20240520T183043.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0,10]
            },
        '240520_2': {
            'UTime': [6.58,10.23], # [ks]
            'description': 'W',
            'beam_energy': '17kV',
            'events': 'Events.20240520T183043.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0,10]
            },
        '240520_3': {
            'UTime': [10.25,11.25], # [ks]
            'description': 'no W',
            'beam_energy': '15kV',
            'events': 'Events.20240520T183043.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0,10]
            },
        '240520_4': {
            'UTime': [11.27,38.51], # [ks]
            'description': 'W',
            'beam_energy': '15kV',
            'events': 'Events.20240520T183043.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.7,10]
            },
        '240520_5': {
            'UTime': [38.67,47.5], # [ks]
            'description': 'calib',
            'beam_energy': '8kV',
            'events': 'Events.20240520T183043.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0,10]
            },

        #### --- May 21st, 2024 --- ####
        '240521_0': {
            'UTime': [0,4.78], # [ks]
            'description': 'background',
            'beam_energy': '17kV' ,
            'events': 'Events.20240521T161246.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0,10]
            },
        '240521_1': {
            'UTime': [5.32,11.68], # [ks]
            'description': 'W',
            'beam_energy': '15kV' ,
            'events': 'Events.20240521T161246.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.6,10]
            },
        '240521_2': {
            'UTime': [11.89, 15.09], # [ks]
            'description': 'W, Be window',
            'beam_energy': '15kV' ,
            'events': 'Events.20240521T161246.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.7,10]
            },
        '240521_3': {
            'UTime': [15.67, 38.23], # [ks]
            'description': 'W, Be window',
            'beam_energy': '15kV' ,
            'events': 'Events.20240521T161246.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [1.0,10]
            },
        '240521_4': {
            'UTime': [38.76, 39.70], # [ks]
            'description': 'W, Be window',
            'beam_energy': '7.1kV' ,
            'events': 'Events.20240521T161246.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [1.0,10]
            },
        '240521_5': {
            'UTime': [40.59, 48.55], # [ks]
            'description': 'W, be window',
            'beam_energy': '7.15kV' ,
            'events': 'Events.20240521T161246.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.6,10]
            },
        '240521_6': {
            'UTime': [48.81, 57.61], # [ks]
            'description': 'calib',
            'beam_energy': '17kV' ,
            'events': 'Events.20240521T161246.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0,10]
            },

        #### --- May 22nd, 2024 --- ####
        '240522a_0': {
            'UTime': [0,4.78], # [ks]
            'description': 'RTS calib',
            'beam_energy': '17kV' ,
            'events': 'Events.20240522T082101-RTS.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0,10]
            },
        '240522b_0': {
            'UTime': [0,4.85], # [ks]
            'description': 'background',
            'beam_energy': '15kV' ,
            'events': 'Events.20240522T162932.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0,10]
            },
        '240522b_1': {
            'UTime': [5.62,8.89], # [ks]
            'description': 'W',
            'beam_energy': '15kV' ,
            'events': 'Events.20240522T162932.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.8,10]
            },
        '240522b_2': {
            'UTime': [9.11, 11.25], # [ks]
            'description': 'W, open EUV',
            'beam_energy': '15kV' ,
            'events': 'Events.20240522T162932.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.7,10]
            },
        '240522b_3': {
            'UTime': [11.35, 33.93], # [ks]
            'description': 'W, close EUV',
            'beam_energy': '7.15kV' ,
            'events': 'Events.20240522T162932.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.5,10]
            },
        '240522b_4': {
            'UTime': [34.02, 46.89], # [ks]
            'description': 'W',
            'beam_energy': '7.28kV' ,
            'events': 'Events.20240522T162932.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.6,10]
            },
        '240522b_5': {
            'UTime': [47.01, 53.61], # [ks]
            'description': 'calib',
            'beam_energy': '15kV' ,
            'events': 'Events.20240522T162932.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0,10]
            },

        #### --- May 23rd, 2024 --- ####
        '240523a_0': {
            'UTime': [0,4.78], # [ks]
            'description': 'RTS calib',
            'beam_energy': '17kV' ,
            'events': 'Events.20240523T073223-RTS.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0,10]
            },
        '240523b_0': {
            'UTime': [1009, 1013.43], # [ks]
            'description': 'background',
            'beam_energy': '15kV' ,
            'events': 'Events.20240523T161003.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0,10]
            },
        '240523b_1': {
            'UTime': [1013.78, 1017.22], # [ks]
            'description': 'W, Be window',
            'beam_energy': '15kV' ,
            'events': 'Events.20240523T161003.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.9,10]
            },
        '240523b_2': {
            'UTime': [1017.45, 1024.05], # [ks]
            'description': 'W, Be window',
            'beam_energy': '7.28kV' ,
            'events': 'Events.20240523T161003.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.5,10]
            },
        '240523b_3': {
            'UTime': [1024.13,1039.01], # [ks]
            'description': 'W, Be window',
            'beam_energy': '7.41kV' ,
            'events': 'Events.20240523T161003.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.6,10]
            },
        '240523b_4': {
            'UTime': [1039.62, 1040.99], # [ks]
            'description': 'W, Be window',
            'beam_energy': '4kV' ,
            'events': 'Events.20240523T161003.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.1,10]
            },
        '240523b_5': {
            'UTime': [1041.28, 1047.68], # [ks]
            'description': 'W, Be window',
            'beam_energy': '4kV' ,
            'events': 'Events.20240523T161003.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.1,10]
            },
        '240523b_6': {
            'UTime': [1047.88,1048.66], # [ks]
            'description': 'W, Be window',
            'beam_energy': '4.04kV' ,
            'events': 'Events.20240523T161003.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.1,10]
            },
        '240523b_7': {
            'UTime': [1053.10, 1053.67], # [ks]
            'description': 'W, Be window',
            'beam_energy': '4.04kV' ,
            'events': 'Events.20240523T161003.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.1,10]
            },
        '240523b_8': {
            'UTime': [1053.73, 1057.46], # [ks]
            'description': 'W, Be window',
            'beam_energy': '4.04kV' ,
            'events': 'Events.20240523T161003.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.1,10]
            },
        '240523b_9': {
            'UTime': [1057.53, 1065], # [ks]
            'description': 'calib',
            'beam_energy': '7.41kV' ,
            'events': 'Events.20240523T161003.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0,10]
            },

        #### --- May 24th, 2024 --- ####
        '240524a_0': {
            'UTime': [0,4.05], # [ks]
            'description': 'Background',
            'beam_energy': '15kV' ,
            'events': 'Events.20240524T163222.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0,10]
            },
        '240524a_1': {
            'UTime': [4.46,7.67], # [ks]
            'description': 'W, Be window',
            'beam_energy': '15kV' ,
            'events': 'Events.20240524T163222.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.9,10]
            },
        '240524b_0': {
            'UTime': [0.95,17.08], # [ks]
            'description': 'W, Be window',
            'beam_energy': '7.54kV' ,
            'events': 'Events.20240524T185245.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.5,10]
            },
        '240524b_1': {
            'UTime': [17.23, 35.13], # [ks]
            'description': 'W, Be window',
            'beam_energy': '7.67kV' ,
            'events': 'Events.20240524T185245.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0.4,10]
            },
        '240524b_2': {
            'UTime': [35.28, 44], # [ks]
            'description': 'calib',
            'beam_energy': '7.54kV' ,
            'events': 'Events.20240524T185245.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0,10]
            },

        #### --- May 25th, 2024 --- ####
        '240525_0': {
            'UTime': [0,4.78], # [ks]
            'description': 'RTS calib',
            'beam_energy': '17kV' ,
            'events': 'Events.20240525T075405-RTS.txt',
            'calib': 'Gain240515T00_gain.itx',
            'EBIT_Phase': [0,10]
            },
        } 

    # Output
    return cases[case]