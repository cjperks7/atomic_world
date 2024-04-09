'''

Module to handle loading WEST data

cjperks
April 1st, 2024

'''

# Modules
try:
    import tofu_west as tfw
    import imas
    import imas_west
    import pywed as pw
except:
    print('No access to WEST data')

import numpy as np
import matplotlib.pyplot as plt
import sys, os

plt.rcParams.update({'font.size': 14})

plt.switch_backend('Qt5Agg')
plt.ion()

__all__ = [
    'get_west',
    'plt_west',
    'save_west',
    ]

###############################################
#
#               Main
#
###############################################

def get_west(
    dout = None,
    shot = None,
    quants = None,
    # Plotting controls
    plt = False,
    t_plt=None,
    # Saving controls
    save = True,
    save_path = None,
    save_name = None,
    ):

    # If empty, load in everything
    if quants is None:
        quants = [
            'eq',
            'profs',
            'interf',
            'ece',
            'icrf',
            'rad',
            'ohm',
            'xics',
            'isotope_ratio',
            'Zeff',
            ]

    if dout is None:
        dout = {}
        dout['exp'] = {}
        dout['shot'] = int(shot)
        dout['device'] = 'WEST'

        # Start time
        dout['t_ignitron'] = pw.tsmat(dout['shot'], 'IGNITRON|1')

    # Loads equilibrium data
    if 'eq' in quants:
        dout = _get_eq(
            dout=dout
            )

    # Loads interferometry data
    if 'interf' in quants:
        dout = _get_interf(
            dout = dout
            )

    # Loads ece data
    if 'ece' in quants:
        dout = _get_ece(
            dout=dout
            )

    # Loads profile data
    if 'profs' in quants:
        dout = _get_profs(
            dout=dout
            )

    # Loads ICRF data
    if 'icrf' in quants:
        dout = _get_icrf(
            dout=dout
            )

    # Loads bolometry data
    if 'rad' in quants:
        dout = _get_rad(
            dout=dout
            )

    # Loads ohmic power
    if 'ohm' in quants:
        dout = _get_ohm(
            dout=dout
            )

    # Loads isotopic ratio
    if 'isotope_ratio' in quants:
        dout = _get_isotope(
            dout=dout
            )

    # Loads Zeff from visible bremsstrahlung
    if 'Zeff' in quants:
        dout = _get_zeff(
            dout=dout
            )

    # Plotting
    if plt:
        plt_west(
            dout=dout,
            t_plt=t_plt
            )

    # Saving
    if save:
        save_west(
            dout=dout,
            save_path=save_path,
            save_name=save_name,
            )

    # Output
    return dout


# Saves data to .npz
def save_west(
    dout=None,
    save_path=None,
    save_name=None,
    ):

    # Error check
    if save_path is None:
        save_path = os.getcwd()

    # Saves data
    path = os.path.join(save_path, save_name+'.npz')
    np.savez(path, dout=dout)


##############################################
#
#               Plotting
#
##############################################

def plt_west(
    dout=None,
    t_plt=None, # [s]
    nlvls = 30,
    ):

    fig, ax = plt.subplots(2,4)
    fig.tight_layout(pad=0.25)

    # Plot interferomtry
    if 'interf' in dout['exp'].keys():
        ax[0,0].plot(
            dout['exp']['interf']['time_s'],
            dout['exp']['interf']['ne_m2']
            )

        ax[0,0].grid('on')
        ax[0,0].set_xlabel('time [s]')
        ax[0,0].set_ylabel('line-integrated ne [m^-2]')

        ax[0,0].set_ylim(
            0,
            np.nanmax(dout['exp']['interf']['ne_m2'].flatten())*1.1
            )
        ax[0,0].set_xlim(0,10)

    # Plots ECE
    if 'ece' in dout['exp'].keys() and len(dout['exp']['ece']['Te_eV']) > 1:
        ax[0,1].plot(
            dout['exp']['ece']['time_s'],
            dout['exp']['ece']['Te_eV']
            )

        ax[0,1].grid('on')
        ax[0,1].set_xlabel('time [s]')
        ax[0,1].set_ylabel('Te [eV]')

        ax[0,1].set_ylim(
            0,
            np.nanmax(dout['exp']['ece']['Te_eV'].flatten())*1.1
            )
        ax[0,1].set_xlim(0,10)

    # PLots H/H+D
    if 'isotope' in dout['exp'].keys():
        for spec in dout['exp']['isotope'].keys():
            ax[0,2].plot(
                dout['exp']['isotope'][spec]['time_s'],
                dout['exp']['isotope'][spec]['H/H+D']*100,
                label = spec
                )

            if isinstance(t_plt, list):
                t_ind = np.where(
                    (dout['exp']['isotope'][spec]['time_s'] >= t_plt[0])
                    & (dout['exp']['isotope'][spec]['time_s'] <= t_plt[1])
                    )[0]

                H_HD_mean = np.mean(dout['exp']['isotope'][spec]['H/H+D'][t_ind])*100
                ax[0,2].plot([t_plt[0],t_plt[1]], [H_HD_mean, H_HD_mean], '--', linewidth=5)



        ax[0,2].set_xlabel('time [s]')
        ax[0,2].set_ylabel('H/H+D [%]')
        leg = ax[0,2].legend()
        leg.set_draggable('on')
        ax[0,2].grid('on')
        ax[0,2].set_xlim(0,10)

    # Plots powers
    if 'icrf' in dout['exp'].keys():
        for ch in np.arange(dout['exp']['icrf']['time_s'].shape[1]):
            ax[1,1].plot(
                dout['exp']['icrf']['time_s'][:,ch],
                dout['exp']['icrf']['Prf_kW'][:,ch],
                label = 'IC antenna %i, %0.1f MHz'%(
                    ch, dout['exp']['icrf']['freq_MHz'][ch]
                    )
                )

    if 'rad' in dout['exp'].keys():
        ax[1,1].plot(
            dout['exp']['rad']['time_s'],
            dout['exp']['rad']['Prad_tot_kW'],
            label = 'Radiated power'
            )

    if 'ohm' in dout['exp'].keys():
        ax[1,1].plot(
            dout['exp']['ohm']['time_s'],
            dout['exp']['ohm']['Pohm_tot_kW'],
            label = 'Ohmic power'
            )

    ax[1,1].grid('on')
    ax[1,1].set_xlabel('time [s]')
    ax[1,1].set_ylabel('power [kW]')
    leg1 = ax[1,1].legend()
    leg1.set_draggable('on')
    ax[1,1].set_xlim(0,10)

    # Plots equilibrium
    if 'eq' in dout['exp'].keys():
        if t_plt is None:
            t_ind = int(len(dout['exp']['eq']['time_s'])/2)

        elif isinstance(t_plt,float):
            t_ind = np.argmin(abs(
                t_plt - dout['exp']['eq']['time_s']
                ))

        elif isinstance(t_plt,list):
            t_ind = np.argmin(abs(
                np.mean(t_plt) - dout['exp']['eq']['time_s']
                ))

        # Plots stored energy
        ax[1,2].plot(
            dout['exp']['eq']['time_s'],
            dout['exp']['eq']['Wmhd_kJ'],
            color = 'blue'
            )
        ax[1,2].set_xlabel('time [s]')
        ax[1,2].set_ylabel('Wmhd [kJ]', color='blue')
        ax[1,2].grid('on')
        ax[1,2].set_xlim(0,10)

        ax2 = ax[1,2].twinx()
        ax2.plot(
            dout['exp']['eq']['time_s'],
            abs(dout['exp']['eq']['Ip_kA']),
            color = 'red'
            )
        ax2.set_ylabel('Ip [kA]', color='red')

        lvls = np.linspace(
            np.nanmin(dout['exp']['eq']['psi_vals'][t_ind, ...]),
            np.nanmax(dout['exp']['eq']['psi_vals'][t_ind, ...]),
            nlvls
            )

        ax[1,0].contour(
            dout['exp']['eq']['psi_R'],
            dout['exp']['eq']['psi_Z'],
            np.squeeze(dout['exp']['eq']['psi_vals'][t_ind, ...]),
            levels = lvls,
            colors = 'blue',
            linestyles ='-',
            linewidths = 0.5
            )

        ax[1,0].contour(
            dout['exp']['eq']['psi_R'],
            dout['exp']['eq']['psi_Z'],
            np.squeeze(dout['exp']['eq']['psi_vals'][t_ind, ...]),
            levels = (dout['exp']['eq']['psi_LCFS'][t_ind],),
            colors = 'tab:red',
            linestyles ='-',
            )

        ax[1,0].plot(
            dout['exp']['eq']['wall_R'],
            dout['exp']['eq']['wall_Z'],
            color='k'
            )

        ax[1,0].grid('on')
        ax[1,0].set_xlabel('R [m]')
        ax[1,0].set_ylabel('Z [m] (t=%1.2f s)'%(dout['exp']['eq']['time_s'][t_ind]))

    # Plots fitted profiles
    if 'profs' in dout['exp'].keys() and 'rho_p' in dout['exp']['profs'].keys():
        if t_plt is None:
            t_ind = int(len(dout['exp']['profs']['time_s'])/2)

        elif isinstance(t_plt, float):
            t_ind = np.argmin(abs(
                t_plt - dout['exp']['profs']['time_s']
                ))

        elif isinstance(t_plt,list):
            t_ind = np.where(
                (dout['exp']['profs']['time_s'] >= t_plt[0])
                & (dout['exp']['profs']['time_s'] <= t_plt[1])
                )[0]
            

        ax[0,3].plot(
            dout['exp']['profs']['rho_p'][t_ind,:].T,
            dout['exp']['profs']['ne_m3'][t_ind,:].T
            )

        ax[0,3].grid('on')
        ax[0,3].set_xlabel('rho_pol (sq. pol. flux) [norm]')
        

        ax[1,3].plot(
            dout['exp']['profs']['rho_p'][t_ind,:].T,
            dout['exp']['profs']['Te_eV'][t_ind,:].T
            )

        ax[1,3].grid('on')
        ax[1,3].set_xlabel('rho_pol (sq. pol. flux) [norm]')

        if isinstance(t_plt, float):
            ax[0,3].set_ylabel('ne [m^-3] (t=%1.2f s)'%(dout['exp']['profs']['time_s'][t_ind]))
            ax[1,3].set_ylabel('Te [eV] (t=%1.2f s)'%(dout['exp']['profs']['time_s'][t_ind]))
        elif isinstance(t_plt, list):
            ax[0,3].set_ylabel('ne [m^-3] (t=%1.2f -%1.2f s)'%(t_plt[0], t_plt[1]))
            ax[1,3].set_ylabel('Te [eV] (t=%1.2f -%1.2f s)'%(t_plt[0], t_plt[1]))

    fig.suptitle(dout['shot'])

    fig.show()
        


##############################################
#
#               Utilities
#
##############################################

# Loads isotopic ratio
def _get_isotope(
    dout=None
    ):

    # Init
    dout['exp']['isotope'] = {}
    data = dout['exp']['isotope'] 

    # Interferometry IDS
    vis = imas_west.get(dout['shot'], 'spectrometer_visible') 

    # Name of spectrometers
    specs = ['LODIVOU15', 'INBUM04', 'LODIVIN19']

    for ch in np.arange(len(vis.channel)):
        if vis.channel[ch].name in specs:
            data[vis.channel[ch].name] = {}
            data[vis.channel[ch].name]['H/H+D'] = vis.channel[ch].isotope_ratios.isotope[0].density_ratio
            data[vis.channel[ch].name]['time_s'] = vis.channel[ch].isotope_ratios.time - dout['t_ignitron']
        else:
            continue

    # Output
    return dout

# Loads interferometry data
def _get_interf(
    dout=None
    ):

    # Init
    dout['exp']['interf'] = {}
    data = dout['exp']['interf'] 

    # Interferometry IDS
    interf = imas_west.get(dout['shot'], 'interferometer')

    # Time basis
    data['time_s'] = interf.time - dout['t_ignitron'] 

    # Number of channels
    nch_int = len(interf.channel)

    # Loads line-integrated interferometry data
    data['ne_m2'] = np.zeros((len(data['time_s']), nch_int)) # dim(ntime, nch)
    for ch in np.arange(nch_int):
        if len(interf.channel[ch].n_e_line.data) > 0:
            data['ne_m2'][:,ch] = interf.channel[ch].n_e_line.data # [m^-2]

    # Output
    return dout

# Loads ECE data
def _get_ece(
    dout=None
    ):

    # Init
    dout['exp']['ece'] = {}
    data = dout['exp']['ece']

    # ECE IDS
    ece = imas_west.get(dout['shot'], 'ece')

    # Time basis
    data['time_s'] = ece.time - dout['t_ignitron']    

    # Number of channels
    nch_ece = len(ece.channel)

    # Loads ECE data
    data['Te_eV'] = np.zeros((len(data['time_s']), nch_ece)) # dim(ntime, nch)
    for ch in np.arange(nch_ece):
        data['Te_eV'][:,ch] = ece.channel[ch].t_e.data

    # Output
    return dout

# Loads fitted core profiles
def _get_profs(
    dout=None,
    ):

    # Init
    dout['exp']['profs'] = {}
    data = dout['exp']['profs']

    # Core profiles IDS
    corep = imas_west.get(dout['shot'], 'core_profiles')

    # Time basis
    data['time_s'] = corep.time - dout['t_ignitron'] 

    # Error check
    if len(corep.profiles_1d) == 0:
        return dout

    # Radial basis
    nrho = len(corep.profiles_1d[-1].grid.rho_pol_norm)

    # Loads fitted data
    data['rho_p'] = np.zeros((len(data['time_s']), nrho)) # dim(ntime,nrho)
    data['ne_m3'] = np.zeros((len(data['time_s']), nrho)) # dim(ntime,nrho)
    data['Te_eV'] = np.zeros((len(data['time_s']), nrho)) # dim(ntime,nrho)

    for tt in np.arange(len(data['time_s'])):
        data['rho_p'][tt,:] = corep.profiles_1d[tt].grid.rho_pol_norm

        data['ne_m3'][tt,:] = corep.profiles_1d[tt].electrons.density

        data['Te_eV'][tt,:] = corep.profiles_1d[tt].electrons.temperature

    # Output
    return dout

# Loads ICRF data
def _get_icrf(
    dout=None
    ):

    # Init
    dout['exp']['icrf'] = {}
    data = dout['exp']['icrf']

    # ICRF IDS
    icrf = imas_west.get(dout['shot'], 'ic_antennas')

    # Number of antennas
    nant = len(icrf.antenna)
    ntime = len(icrf.antenna[0].power_launched.time)

    # Loads launched RF power data
    data['time_s'] = np.zeros((ntime, nant)) # dim(ntime, nantenna)
    data['Prf_kW'] = np.zeros((ntime, nant)) # dim(ntime, nantenna)
    data['freq_MHz'] = np.zeros(nant) # dim(nantenna,)

    for aa in np.arange(nant):
        data['time_s'][:,aa] = icrf.antenna[aa].power_launched.time - dout['t_ignitron']

        data['Prf_kW'][:,aa] = icrf.antenna[aa].power_launched.data/1e3

        data['freq_MHz'][aa] = np.mean(icrf.antenna[aa].frequency.data)/1e6

    # Output
    return dout

# Loads bolometry data
def _get_rad(
    dout=None
    ):

    # Init
    dout['exp']['rad'] = {}
    data = dout['exp']['rad']

    # ICRF IDS
    bolo = imas_west.get(dout['shot'], 'bolometer')

    # Time basis
    data['time_s'] = bolo.time - dout['t_ignitron'] # dim(ntime,)

    # Total radiated power
    data['Prad_tot_kW'] = bolo.power_radiated_total/1e3 # dim(ntime,)

    # Output
    return dout

# Loads ohmic power
def _get_ohm(
    dout=None
    ):

    # Init
    dout['exp']['ohm'] = {}
    data = dout['exp']['ohm']

    # summary IDS
    summ = imas_west.get(dout['shot'], 'summary')

    # Time basis
    data['time_s'] = summ.time - dout['t_ignitron'] # dim(ntime,)

    # Total radiated power
    data['Pohm_tot_kW'] = summ.global_quantities.power_ohm.value/1e3 # dim(ntime,)

    # Output
    return dout

# Get 2D equilibrium data
def _get_eq(
    dout=None
    ):

    # Init
    dout['exp']['eq'] = {}
    data = dout['exp']['eq']

    # Equilibrium IDS
    equi = imas_west.get(dout['shot'], 'equilibrium')
    wall = imas_west.get(dout['shot'], 'wall')

    # Time basis
    data['time_s'] = equi.time - dout['t_ignitron'] # dim(ntime,)

    # Limiter contour
    data['wall_R'] = wall.description_2d[0].limiter.unit[0].outline.r
    data['wall_Z'] = wall.description_2d[0].limiter.unit[0].outline.z

    # Closes contour
    data['wall_R'] = np.append(data['wall_R'], data['wall_R'][0])
    data['wall_Z'] = np.append(data['wall_Z'], data['wall_Z'][0])

    # LCFS
    data['psi_LCFS'] = equi.boundary.psi # dim(ntime,)

    # Poloidal flux contours
    data['psi_vals'] = equi.interp2D.psi # dim(ntime, nR, nZ)
    data['psi_R'] = equi.interp2D.r # dim(nR, nZ)
    data['psi_Z'] = equi.interp2D.z # dim(nR, nZ)

    # Total B-field
    data['B_tot_T'] = np.sqrt(
        equi.interp2D.b_field_r**2
        + equi.interp2D.b_field_z**2
        + equi.interp2D.b_field_tor**2
        ) # dim(ntime, nR, nZ)

    # Stored energy
    data['Wmhd_kJ'] = equi.global_quantities.w_mhd/1e3

    # Plasma current
    data['Ip_kA'] = equi.global_quantities.ip/1e3

    # Output
    return dout

# Gets Zeff from visible bremsstrahlung
def _get_zeff(
    dout = None
    ):

     # Init
    dout['exp']['vis_brem'] = {}
    data = dout['exp']['vis_brem']

    # Equilibrium IDS
    vis = imas_west.get(dout['shot'], 'bremsstrahlung_visible', 0,1)

    # Error check
    if len(vis.channel) == 0:
        return dout

    # Number of antennas
    nch = len(vis.channel)
    ntime = len(vis.channel[0].zeff_line_average.time)

    # Loads launched RF power data
    data['time_s'] = np.zeros((ntime, nch)) # dim(ntime, nchannels)
    data['Zeff'] = np.zeros((ntime, nch)) # dim(ntime, nchannels)

    for ch in np.arange(nch):
        data['time_s'][:,ch] = vis.channel[ch].zeff_line_average.time - dout['t_ignitron']

        data['Zeff'][:,ch] = vis.channel[ch].zeff_line_average.data

    # Output
    return dout
