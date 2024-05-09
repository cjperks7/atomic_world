'''

_calc_xics.py is a script meant to manage processing XICS data from WEST

cjperks
Apr 10, 2024

'''

# Modules
import numpy as np
import os
import matplotlib.pyplot as plt

from transport_world.run_AURORA import calc_imp as ci
from transport_world.run_AURORA import calc_spectra as csp
from transport_world.plot_utils import plt_slider
from atomic_world.run_TOFU.run_TOFU_WEST import get_west

__all__ = [
    'get_pumpout_est'
    ]

########################################################
#
#           Main
#
########################################################

# Calculates estimate of density pumpout from line brightness
def get_pumpout_est(
    # Machine settings
    device = None,
    shot = None,
    # Data settings
    path = None,
    # Line settings
    sp = None,
    cs = None,
    source = 'raw',
    dsource = None,
    dt = None,
    # Analysis settings
    davg = None,
    ):

    # lie brightness data
    if device == 'WEST':
        line, Te_tan, ne_tan, fz, time, dout = _get_line_WEST(
            shot = shot,
            path = path,
            sp = sp,
            cs = cs,
            source = 'raw',
            dsource = dsource,
            dt = dt,
            davg=davg,
            ) # dim(ntime,)

    elif device == 'CMOD':
        print('NOT IMPLEMENTED YET!!!')
        sys.exit(1)

    print(dout)

    if davg is not None:
        H_HD = 0
        cnt = 0
        for key in dout['H/H+D'].keys():
            for hh in dout['H/H+D'][key]:
                cnt += 1
                H_HD += hh
        H_HD /= cnt
        print(H_HD*100)

        dB_Ar = 1-dout['line'][1]/dout['line'][0]
        print(dB_Ar*100)

        PEC1 = csp.get_PEC(
            sp='Ar',
            cs=cs,
            Te_eV = np.r_[dout['Te0'][0]],
            ne_cm3 = np.r_[dout['ne0'][0]/1e6],
            lamb0 = 3.949075,
            dlamb = 0.001
            )

        PEC2 = csp.get_PEC(
            sp='Ar',
            cs=cs,
            Te_eV = np.r_[dout['Te0'][1]],
            ne_cm3 = np.r_[dout['ne0'][1]/1e6],
            lamb0 = 3.949075,
            dlamb = 0.001
            )

        nAr_2o1 = (
            (1-dB_Ar)
            * dout['fz0'][0]*dout['ne0'][0]*PEC1[0,1]
            / (dout['fz0'][1]*dout['ne0'][1]*PEC2[0,1])
            )
        print(nAr_2o1*100)
        print((1-nAr_2o1)*100)


        Te_arr = np.logspace(np.log10(4e2), np.log10(1.5e3), 500)
        ne_arr = np.ones_like(Te_arr)*1e14

        PEC_arr = csp.get_PEC(
            sp='Ar',
            cs=cs,
            Te_eV = Te_arr,
            ne_cm3 = ne_arr,
            lamb0 = 3.949075,
            dlamb = 0.001
            )
        fz_arr = ci.get_ion_bal(
            sp=sp,
            Te_eV = Te_arr,
            ne_cm3 = ne_arr
            )

        fig, ax = plt.subplots()

        ax.plot(
            Te_arr,
            PEC_arr[:,0]*fz_arr[:,cs-1],
            label = 'ioniz'
            )
        ax.plot(
            Te_arr,
            PEC_arr[:,1]*fz_arr[:,cs],
            label = 'excit'
            )
        ax.plot(
            Te_arr,
            PEC_arr[:,2]*fz_arr[:,cs+1],
            label = 'recomb'
            )
        ax.set_xlabel('Te [eV]')
        ax.set_ylabel(r'$\epsilon\sim f^{Z+}_{V=D=0}*PEC^{Z+}$')
        ax.legend()
        ax.grid('on')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(1e-17,1e-11)

    

########################################################
#
#           Utilities  (WEST)
#
########################################################

# Analyze line data from WEST
def _get_line_WEST(
    shot = None,
    path = None,
    sp = None,
    cs = None,
    source = None,
    dsource = None,
    dt = None,
    davg = None,
    ):

    if path is None:
        path = '/home/cjperks/2201_Pumpout/WEST/C9_results'

    # gets summary data
    dsum = _get_WEST_summary(shot=shot, path=path)

    if davg is not None:
        get_west(
            dout = dsum,
            quants = [],
            plt = True,
            t_plt = davg['t2'],
            save = False,
            )


    # Data source
    if source == 'raw':
        line, dsource, time, Te_xics = _get_line_WEST_raw(
            shot=shot,
            path=path,
            dsource = dsource,
            dsum = dsum,
            ) # dim(ntime,)
    else:
        print("NOT IMPLEMENTED YET!!!!")
        sys.exit(1)

    # Time shift if necessary
    if dt is not None:
        dsource, time = _shift_time(
            line=line,
            Te_xics=Te_xics,
            dsource=dsource,
            time=time,
            dsum=dsum,
            dt=dt
            )

    # Tangency density/temperature
    if 'rho_p' in dsum['exp']['profs'].keys():
        ne_tan = dsum['exp']['profs']['ne_m3'][:,0]
        Te_tan = dsum['exp']['profs']['Te_eV'][:,0]

        # Charge state distribution from ionization balance
        fz = ci.get_ion_bal(
            sp=sp,
            Te_eV = Te_tan,
            ne_cm3 = ne_tan/1e6
            )
    else:
        ne_tan = Te_tan = fz = 0

    # Plotting
    if dsource['plt']:
        fig, ax = plt.subplots()

        Prf = (
            np.sum(dsum['exp']['icrf']['Prf_kW'],axis=1)
            /np.max(np.sum(dsum['exp']['icrf']['Prf_kW'],axis=1))
            *np.max(line)*0.5
            )
        Te_norm = (
            1
            /np.max(Te_tan)
            *np.max(line)*0.5
            )
        Te = Te_tan*Te_norm
        ne_norm = (
            1
            /np.max(ne_tan)
            *np.max(line)*0.5
            )
        ne = ne_tan*ne_norm

        ax.plot(
            time,
            line,
            label = 'w-line'
            )
        ax.plot(
            dsum['exp']['icrf']['time_s'][:,0],
            Prf,
            label = 'ICRF [arb.]'
            )

        if 'rho_p' in dsum['exp']['profs'].keys():
            ax.plot(
                dsum['exp']['profs']['time_s'],
                Te,
                label = 'Te0 [arb.]'
                )
            ax.plot(
                dsum['exp']['profs']['time_s'],
                ne,
                label = 'ne0 [arb.]'
                )

        ax.set_xlim(0,10)
        ax.set_ylim(0, np.max(line)*1.1)
        ax.grid('on')
        ax.set_title(shot)
        ax.legend()

    if davg is not None:
        dout = {}
        # Perfroms averaging
        ind1 = np.where(
            (time >= davg['t1'][0])
            & (time <= davg['t1'][1])
            )[0]
        ind2 = np.where(
            (time >= davg['t2'][0])
            & (time <= davg['t2'][1])
            )[0]

        dout['line'] = [
            np.mean(line[ind1]),
            np.mean(line[ind2])
            ]

        # Perfroms averaging
        ind1 = np.where(
            (dsum['exp']['profs']['time_s'] >= davg['t1'][0])
            & (dsum['exp']['profs']['time_s'] <= davg['t1'][1])
            )[0]
        ind2 = np.where(
            (dsum['exp']['profs']['time_s'] >= davg['t2'][0])
            & (dsum['exp']['profs']['time_s'] <= davg['t2'][1])
            )[0]

        dout['ne0'] = [
            np.mean(ne_tan[ind1]),
            np.mean(ne_tan[ind2])
            ]
        dout['Te0'] = [
            np.mean(Te_tan[ind1]),
            np.mean(Te_tan[ind2])
            ]
        dout['fz0'] = [
            np.mean(fz[ind1, cs]),
            np.mean(fz[ind2, cs])
            ]

        dout['H/H+D'] = {}
        for key in dsum['exp']['isotope'].keys():
            ind1 = np.where(
                (dsum['exp']['isotope'][key]['time_s'] >= davg['t1'][0])
                & (dsum['exp']['isotope'][key]['time_s'] <= davg['t1'][1])
                )[0]
            ind2 = np.where(
                (dsum['exp']['isotope'][key]['time_s'] >= davg['t2'][0])
                & (dsum['exp']['isotope'][key]['time_s'] <= davg['t2'][1])
                )[0]
            dout['H/H+D'][key] = [
                np.mean(dsum['exp']['isotope'][key]['H/H+D'][ind1]),
                np.mean(dsum['exp']['isotope'][key]['H/H+D'][ind2])
                ]

        ax.plot(
            davg['t1'],
            [dout['line'][0], dout['line'][0]],
            'r--'
            )
        ax.plot(
            davg['t2'],
            [dout['line'][1], dout['line'][1]],
            'r--'
            )
        ax.plot(
            davg['t1'],
            [dout['ne0'][0]*ne_norm, dout['ne0'][0]*ne_norm],
            'b--'
            )
        ax.plot(
            davg['t2'],
            [dout['ne0'][1]*ne_norm, dout['ne0'][1]*ne_norm],
            'b--'
            )
        ax.plot(
            davg['t1'],
            [dout['Te0'][0]*Te_norm, dout['Te0'][0]*Te_norm],
            'g--'
            )
        ax.plot(
            davg['t2'],
            [dout['Te0'][1]*Te_norm, dout['Te0'][1]*Te_norm],
            'g--'
            )
    
    else:
        dout = None

    # Output
    return line, Te_tan, ne_tan, fz, time, dout


# Analyze raw data from WEST
def _get_line_WEST_raw(
    shot = None,
    path = None,        # Common path to data
    dsource = None,     # Controls to handle source data
    dsum = None,        # Shot summary data
    ):

    # Loads raw data
    draw = np.load(
        os.path.join(path,'XICS',str(shot)+'_raw.npz'),
        allow_pickle=True
        )['dout'][()]

    # Image data
    cnts = draw['exp']['xics']['raw'] # dim(ntime, nvert, nhorz)

    # Line controls
    if dsource is None:
        dsource = {}
        dsource['plt'] = True
        dsource['pix_hor'] = [385,405]
        dsource['time_s'] = 4.0 # [s]
        dsource['cnt_max'] = 150
        dsource['pix_vert'] = 800

    # Integrates over line profile
    line = np.sum(
        cnts[:,  dsource['pix_vert'], dsource['pix_hor'][0]:dsource['pix_hor'][1]],
        axis =1
        ) # dim(ntime,)

    if dsource['plt']:
        plt_slider(
            xxx=np.arange(0,draw['exp']['xics']['nxi']),
            yyy=draw['exp']['xics']['time_s'],
            dzzz = cnts[:,dsource['pix_vert'],:],
            xxlabel='horz',
            yylabel='time',
            zzlabel='counts',
            xlim = [dsource['pix_hor'][0]-5, dsource['pix_hor'][1]+5],
            ylim = [0,dsource['cnt_max']]
            )

        plt_slider(
            xxx=np.arange(0,draw['exp']['xics']['nxj']),
            yyy=draw['exp']['xics']['time_s'],
            dzzz = cnts[:,:, int(np.mean(dsource['pix_hor']))],
            xxlabel='vert',
            yylabel='time',
            zzlabel='counts',
            ylim = [0,dsource['cnt_max']]
            )

        ind = np.argmin(abs(draw['exp']['xics']['time_s']-dsource['time_s']))
        fig, ax = plt.subplots()
        aa = ax.imshow(cnts[ind,:,:], vmax = 50,vmin=0)
        fig.colorbar(aa,ax=ax)

    # Te from midplane view
    nt = draw['exp']['xics']['rho_t'].shape[0]
    ind = np.argmin(abs(
        draw['exp']['xics']['rho_t'][int(2/3*ind), :]
        ))

    # Output, [counts], dim(ntime,)
    return line, dsource, draw['exp']['xics']['time_s'], draw['exp']['xics']['Te_line_eV'][:,ind]
    

# Gets summary data for this shot
def _get_WEST_summary(
    shot=None,
    path=None,
    ):

    # output
    return np.load(
        os.path.join(path,'shots',str(shot)+'_summary.npz'),
        allow_pickle=True
        )['dout'][()]

# Shifts time basis if needed
def _shift_time(
    line=None,
    Te_xics=None,
    dsource=None,
    time=None,
    dsum=None,
    dt=None
    ):

    # Shift
    time -= dt

    # Documents
    dsource['shift_s'] = dt

    # Plots
    fig,ax=plt.subplots()

    ax.plot(
        time,
        Te_xics/np.mean(Te_xics),
        label = 'XICS'
        )
    if 'rho_p' in dsum['exp']['profs'].keys():
        ax.plot(
            dsum['exp']['profs']['time_s'],
            dsum['exp']['profs']['Te_eV'][:,0]/np.mean(dsum['exp']['profs']['Te_eV'][:,0]),
            label = 'ECE'
            )

    ax.grid('on')
    ax.legend()
    ax.set_xlabel('time [s]')
    ax.set_ylabel('Te [arb.]')

    # Output
    return dsource, time