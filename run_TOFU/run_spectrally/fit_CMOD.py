'''

Workflow to use spectrally to fit XEUS and LoWEUS spectra
on C-Mod

cjperks
Aug 30, 2024

'''

# Modules
import os, sys
import numpy as np
import spectrally as spec

from atomic_world.run_TOFU.run_spectrally.CMOD import _get_VUV as gV
from atomic_world.run_TOFU.run_spectrally.CMOD import _def_pumpout as dP

# CMOD controls
shot = 1140221013

vuvs = ['XEUS', 'LoWEUS']

####################################
#
#           Functions
#
####################################

# Adds full spectrum data
def _add_data(
    coll = None,
    ddata = None,
    system = None,
    case = 'all',
    shot = None,
    alias = None,
    ):

    # Shorthand
    dd = ddata[system]

    # Slices the data if requested
    time = dd['time_s']
    lamb = dd['lamb_AA']
    sig = dd['signal']  
    if case == 'pumpout':
        lims = dP._def_slicing(
            shot = shot,
            system = system,
            )
        indt = np.where(
            (time >= lims['time_s'][0])
            & (time <= lims['time_s'][1])
            )[0]
        indy = np.where(
            (lamb >= lims['lamb_AA'][0])
            & (lamb <= lims['lamb_AA'][1])
            )[0]
        time = time[indt]
        lamb = lamb[indy]
        sig = sig[np.ix_(indt,indy)]

    if alias is not None:
        nt = len(time)
        ny = len(lamb)
        nt_rem = nt%alias

        if nt_rem != 0:
            sig = sig[:-nt_rem, :]
            time = time[:-nt_rem]

        sig = np.sum(
            sig.reshape(nt //alias, alias ,ny),
            axis = 1
            )
        time = time[::alias]

    # Adds time data
    coll.add_ref(
        'nt_'+case+'_'+system, 
        size=len(time)
        )
    coll.add_data(
        't_'+case+'_'+system,
        data = np.linspace(time[0], time[-1], time.size), # Ensures uniformity
        ref = 'nt_'+case+'_'+system,
        units =  's',
        )

    # Adds wavlength data
    coll.add_ref(
        'nlamb_'+case+'_'+system, 
        size = len(lamb)
        )
    coll.add_data(
        'lamb_'+case+'_'+system,
        data = lamb*1e-10,
        ref = 'nlamb_'+case+'_'+system,
        units = 'm'
        )

    # Adds spectrum data
    coll.add_data(
        'spectra_'+case+'_'+system,
        data = sig,
        ref = ('nt_'+case+'_'+system, 'nlamb_'+case+'_'+system),
        units = 'counts'
        )

    # Output
    return coll

# Perfroms spectral fitting
def _get_fits(
    coll = None,
    system = None,
    shot = None,
    ):

    # Gets defaults setting
    case = 'pumpout'
    dset = dP._def_dsettings(shot=shot,system=system)

    # Adds spectral model
    dmodel = {}
    dmodel['bck'] = dset['dmodel']['bck']
    for kk in dset['dmodel'].keys():
        if kk == 'bck':
            continue
        for ii in np.arange(len(dset['dmodel'][kk]['vlamb'])):
            dmodel[kk+str(ii)] = {
                'type': 'gauss',
                'lamb0': dset['dmodel'][kk]['vlamb'][ii],
                'mz': dset['dmodel'][kk]['mz'][ii] # [kg]
                }
    coll.add_spectral_model(
        key='sm_'+system, 
        dmodel=dmodel,
        dconstraints = dset['dconstraints'],
        )

    # Init spectral fit
    coll.add_spectral_fit(
        key='sf_'+system,
        key_model='sm_'+system,
        key_data='spectra_pumpout_'+system,
        key_lamb='lamb_pumpout_'+system,
        dvalid=dset['dvalid'],
        )

    # Lower bounds
    dbounds_low = {}
    if dmodel['bck'] == 'poly':
        dbounds_low['bck_a0'] = 0
    elif dmodel['bck'] == 'exp_lamb':
        dbounds_low['bck_amp'] = 50e-10
        dbounds_low['bck_rate'] = 1e-10

    for kk in dset['dmodel'].keys():
        if kk == 'bck':
            continue
        for ii in np.arange(len(dset['dmodel'][kk]['vlamb'])):
            dbounds_low[kk+str(ii)+'_amp'] = 0
            dbounds_low[kk+str(ii)+'_vccos'] = (
                -1 *dset['dbnds'][kk]['vc_bnd'] 
                / (dset['dmodel'][kk]['vlamb'][ii]*1e10)
                )
            dbounds_low[kk+str(ii)+'_sigma'] = (
                dset['dbnds'][kk]['dlamb']
                ) # [m]

    # Upper bounds
    tmp = coll.ddata['spectra_pumpout_'+system]['data']
    dbounds_up = {}
    if dmodel['bck'] == 'poly':
        dbounds_up['bck_a0'] = 0.5 if system == 'XEUS' else 2
    elif dmodel['bck'] == 'exp_lamb':
        dbounds_up['bck_amp'] = 500e-10
        dbounds_up['bck_rate'] = 100e-10

    for kk in dset['dmodel'].keys():
        if kk == 'bck':
            continue
        for ii in np.arange(len(dset['dmodel'][kk]['vlamb'])):
            dbounds_up[kk+str(ii)+'_amp'] = 1.1*np.max(tmp.flatten())
            dbounds_up[kk+str(ii)+'_vccos'] = (
                dset['dbnds'][kk]['vc_bnd'] 
                / (dset['dmodel'][kk]['vlamb'][ii]*1e10)
                )
            dbounds_up[kk+str(ii)+'_sigma'] = (
                dset['dbnds'][kk]['dlamb'] * dset['dbnds'][kk]['sig_bnd']
                ) # [m]

    # Initial guess
    dx0 = {}
    if dmodel['bck'] == 'poly':
        dx0['bck_a0'] = 0.1 if system == 'XEUS' else 0.1
    elif dmodel['bck'] == 'exp_lamb':
        dx0['bck_amp'] = 200e-10
        dx0['bck_rate'] = 5e-10

    for kk in dset['dmodel'].keys():
        if kk == 'bck':
            continue
        for ii in np.arange(len(dset['dmodel'][kk]['vlamb'])):
            dx0[kk+str(ii)+'_amp'] = 0.3*np.max(tmp.flatten())
            dx0[kk+str(ii)+'_vccos'] = 0
            dx0[kk+str(ii)+'_sigma'] = (
                dset['dbnds'][kk]['dlamb']
                ) # [m]

    # secoComputend spectral model with user-provided scales, bounds and x0
    coll.compute_spectral_fit(
        'sf_'+system,
        # user-provided
        dbounds_low=dbounds_low,
        dbounds_up=dbounds_up,
        dx0=dx0,
        # options
        chain=False, # True= how is initial guess defined -> if samples correlated in time
        verb=2,
        strict=False, # True= if time step doesn't converge stop
        binning = False, # True= by default the number of subbins (10) discretizing the model over each bin; scalar= customize subdivision (longer comp)
        solver = 'scipy.curve_fit'# default= (error bars) 'scipy.curve_fit'; others = 'scipy.least_squares'
        #solver = 'scipy.least_squares',
        )

    # Output
    return coll






####################################
#
#           Functions
#
####################################



# Loads spectrum data
ddata = {}
for vuv in vuvs:
    ddata[vuv] = gV._get_spectrum(
        shot = shot,
        system = vuv
        )

# Init
coll = spec.Collection()

# Path to line labeling data
path2lines = os.path.join(
    '/home/cjperks/work',
    '2201_Pumpout/CMOD',
    'spectrally'
    )

# Looks at full spectrum time trace
case = 'all'
for vuv in vuvs:
    # Adds line labeling data to Collection object
    files = [
        f for f in os.listdir(path2lines) if f.startswith(vuv+'_shifted')
        #f for f in os.listdir(path2lines) if f.startswith(vuv+'_FAC')
        ]
    for fil in files:
        coll.add_spectral_lines_from_file(
            pfe=os.path.join(path2lines, fil)
            )


    coll = _add_data(
        coll = coll,
        system = vuv,
        ddata = ddata,
        case = case
        )

    # Plots
    if False:
        dax = coll.plot_as_array(
            'spectra_'+case+'_'+vuv, 
            keyX='lamb_'+case+'_'+vuv, 
            keyY='t_'+case+'_'+vuv,
            aspect='auto'
            )

# Looks at spectrum time trace for pumpout
case = 'pumpout'
for vuv in vuvs:
    coll = _add_data(
        coll = coll,
        system = vuv,
        ddata = ddata,
        case = case,
        shot = shot,
        alias = 20,
        )

    # Plots
    if False:
        #dax = coll.plot_as_array(
        #    'spectra_'+case+'_'+vuv, 
        #    keyX='lamb_'+case+'_'+vuv, 
        #    keyY='t_'+case+'_'+vuv,
        #    aspect='auto'
        #    )
        dax = coll.plot_spectral_data(
            key_data='spectra_'+case+'_'+vuv, 
            key_lamb='lamb_'+case+'_'+vuv, 
            key_lines=[
                kk for kk in coll.dobj['spect_lines'].keys() if '_'+vuv[0] in kk
                ], 
            lines_labels_color='ion'
            )


colors = {}
colors['LoWEUS'] = {
    **{'Ar'+str(ii): 'b' for ii in ['0']},
    **{'Mo'+str(ii): 'm' for ii in ['0', '1']}
    }
colors['XEUS'] = {
    **{'Ar'+str(ii): 'r' for ii in ['0','2','3','4','5']},
    **{'Ar'+str(jj): 'b' for jj in ['1','6','7','10']},
    **{'Ar'+str(kk): 'g' for kk in ['8','9']},
    **{'Misc'+str(zz): 'k' for zz in [str(ii) for ii in np.arange(1,18)]},
    }

labels = {}
labels['LoWEUS'] = {
    **{'Ar'+str(ii): 'Ar14+' for ii in ['0']},
    **{'Mo'+str(ii): 'Mo30+' for ii in ['0']},
    **{'Mo'+str(ii): 'Mo31+' for ii in ['1']}
    }
labels['XEUS'] = {
    **{'Ar'+str(ii): 'Ar15+' for ii in ['0','2','3','4','5']},
    **{'Ar'+str(jj): 'Ar14+' for jj in ['1','6','7','10']},
    **{'Ar'+str(kk): 'Ar13+' for kk in ['8','9']},
    **{'Misc'+str(zz): '' for zz in [str(ii) for ii in np.arange(18)]},
    }

# Performs spectral fitting
for vuv in ['XEUS']:
    coll = _get_fits(
        coll = coll,
        system = vuv,
        shot = shot,
        )

    if False:
        dax = coll.plot_spectral_fit_input_validity(
            'sf_'+vuv
            )

    if True:
        dax = coll.plot_spectral_fit(
            'sf_'+vuv,
            lines_labels = labels[vuv],
            lines_labels_color = colors[vuv]
            )



db = {}

for vuv in ['XEUS']:#['LoWEUS']:
    moms = coll.get_spectral_model_moments('sf_'+vuv)
    db[vuv] = {}
    fig, ax = plt.subplots()

    keys = [ss for ss in moms.keys() if ss.endswith('_amp')]
    for key in keys:
        kk = key.split('_')[0]

        if 'Misc' in kk:
            continue
        elif 'bck' in kk:
            continue

        db[vuv][kk] = (
            moms[kk+'_integ']['data'][0] - moms[kk+'_integ']['data'][-1]
            )/moms[kk+'_integ']['data'][0]

        ax.plot(
            coll.ddata['t_pumpout_'+vuv]['data'],
            #moms[key]['data'],
            moms[kk+'_integ']['data'],
            #moms[kk+'_integ']['data'] * (moms['Ar0_integ']['data'][0]/moms[kk+'_integ']['data'][0]),
            '*-',
            label = r'%s, %0.2f $\AA$'%(
                kk,
                coll.dobj['spect_model']['sm_'+vuv]['dmodel'][kk]['param']['lamb0']
                *1e10
                #* (1 + moms[key.split('_')[0]+'_vccos']['data'])
                ),
            color = colors[vuv][kk]
            )

        ax.text(
            coll.ddata['t_pumpout_'+vuv]['data'][0],
            #moms[key]['data'][0],
            moms[kk+'_integ']['data'][0],
            #moms['Ar0_integ']['data'][0],
            kk,
            #fontsize=12,
            color = colors[vuv][kk]
        )
        '''
        ax.fill_between(
            coll.ddata['t_pumpout_'+vuv]['data'],
            moms[key+'_min']['data'],
            moms[key+'_max']['data'],
            #moms[key+'_min']['data'] *(
            #    moms[key.split('_')[0]+'_integ']['data']
            #    /moms[key]['data']
            #    ),
            #moms[key+'_max']['data'] *(
            #    moms[key.split('_')[0]+'_integ']['data']
            #    /moms[key]['data']
            #    ),
            alpha = 0.3
            )
        '''
    '''
    ax.plot(
        coll.ddata['t_pumpout_'+vuv]['data'],
        moms['bck_a0']['data'],
        '*-',
        label = 'bckgd'
        )
    ax.fill_between(
            coll.ddata['t_pumpout_'+vuv]['data'],
            moms['bck_a0_min']['data'],
            moms['bck_a0_max']['data'],
            alpha = 0.3
        )
    '''       

    leg = ax.legend(labelcolor='linecolor')
    leg.set_draggable('on')
    ax.grid('on')
    ax.set_xlabel('time [s]')
    #ax.set_ylabel('amp [arb]')
    ax.set_ylabel('integ [arb]')

    ax.set_title(vuv)


'''

import pandas as pd
from scipy.interpolate import interp1d
res = pd.read_csv('LoWEUSresponsivity_LepsonRSI2016.csv')
yys = np.r_[118.6, 130.7, 224.9]
vals = interp1d(res['wavelength [AA]'], res[' ratio calculated/measured'])(yys)



'''