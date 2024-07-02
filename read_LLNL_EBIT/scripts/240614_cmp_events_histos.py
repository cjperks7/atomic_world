'''

Script comparing LLNL EBIT ECS spectra loaded from using
the Events List analysis python workflow v. emailed histograms

Here, use the 7.15kV beam energy measured on May 21st, 2024

cjperks
June 14, 2024

'''

# Modules
import os
from atomic_world.read_LLNL_EBIT.main import load_events as le
from atomic_world.read_LLNL_EBIT.main import read_histos as rh
import ebit_ecs as eb


####################################
#
#       Event List
#
####################################

# File management
file_path = os.path.join(
    '/home/cjperks/work',
    '2008_SPARCxray/LLNL_EBIT',
    '240531_W_data'
    )
events = 'Events.20240521T161246.txt'
calib = 'Gain240515T00_gain.itx'


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
flag_list = [64]

# Filters by pixel and event quality
df_cal_pix = df_cal.loc[
    df_cal['Flags'].isin(flag_list) 
    & df_cal['Pixel'].isin(lp)
    ]

# Filter by UTime
df_ut = df_cal_pix.loc[
    (df_cal_pix['UTime'] >= 40.59) 
    & (df_cal_pix['UTime'] <= 48.55)
    ]

# Filter by EBIT Phase
df_eb = df_ut.loc[
    (df_ut['EbitPhase'] >= 0.6) 
    & (df_ut['EbitPhase'] <= 10)
    ]

'''
le.plt_all(
    df=df_cal_pix,
    title = '7.15kV'
    )
le.plt_all(
    df=df_ut,
    title ='7.15kV',
    )

'''

# Makes histogram
spec_ev = eb.histo(df_eb['PH'], 1,-100,10000)


###########################################
#
#           Histogram ASCII files
#
###########################################

# File management
file_path_hi = os.path.join(
    '/home/cjperks/work',
    '2008_SPARCxray/LLNL_EBIT',
    '240601_W_histos'
    )
name = '/W_Escan_7.15-keV.hist'

# Init
spec_hi = {}
spec_hi['binlow'] = []
spec_hi['binhigh'] = []
spec_hi['hist'] = []
spec_hi['hist_all'] = []

# Reads file
with open(file_path_hi+name, 'r') as ff:
    for line in ff:
        spec_hi['binlow'].append(
            float(line.split()[0])
            )
        spec_hi['binhigh'].append(
            float(line.split()[1])
            )
        spec_hi['hist_all'].append(
            float(line.split()[2])
            )
        spec_hi['hist'].append(
            float(line.split()[3])
            )

ff.close()

# Reformat
spec_hi['binlow'] = np.asarray(spec_hi['binlow'])
spec_hi['binhigh'] = np.asarray(spec_hi['binhigh'])
spec_hi['hist_all'] = np.asarray(spec_hi['hist_all'])
spec_hi['hist'] = np.asarray(spec_hi['hist'])

#############################################
#
#               SCRAM
#
#############################################

# Loads SCRAM spectra
scram = rh.read_SCRAM()

# Specific case to plat
Te_plot = 4e3 # [eV]
ion_plot = 22 # Ti-like

# Maxwellian spectrum to plot
Te_ind = np.argmin(
    abs(scram['table']['Te(eV)'] - Te_plot)
    )
key = 'spec_%i'%(Te_ind)

# Finds particular ion of interest
css = [int(ii.split('-')[-1]) for ii in scram['indspec'][key]['ions']]

try:
    ion_ind = css.index(ion_plot)

    E_sc = scram['indspec']['E_eV']['data']
    spec_sc = scram['indspec'][key]['j(W/ion/eV)'][ion_ind, :]


except:
    E_sc = 0
    spec_sc = 0

spec_sc_11 = scram['indspec']['spec_11']['j(W/ion/eV)'][6, :]



#############################################
#
#               Plotting
#
#############################################

fig, ax = plt.subplots(figsize = (16,8))


ax.step(
    spec_ev['binlow'],
    spec_ev['hist']/np.max(spec_ev['hist']),
    label = 'Event List',
    color = 'red'
    )
ax.step(
    spec_hi['binlow'],
    spec_hi['hist']/np.max(spec_hi['hist']),
    label = 'Histo ASCII',
    color = 'blue'
    )

ax.plot(
    E_sc,
    spec_sc /np.max(spec_sc),
    label = 'SCRAM (s-22)',
    color = 'green'
)

ax.plot(
    E_sc,
    spec_sc_11 /np.max(spec_sc_11),
    label = 'SCRAM (d-22)',
    color = 'm'
)

ax.set_xlim(2000,3000)

ax.legend(labelcolor='linecolor')
ax.set_xlabel('Photon Energy [eV]')
ax.set_ylabel('Counts [arb]')
ax.grid('on')
fig.set_title('7.15kV')