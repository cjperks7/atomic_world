'''

Module to handle loading WEST data

cjperks
April 1st, 2024

'''

# Modules
import tofu_west as tfw
import imas
import imas_west
import pywed as pw

__all__ = [
    'get_west',
    'plt_west',
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
    plt = False,
    t_plt=None,
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
        dout = _get_intef(
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

    # Plotting
    if plt:
        plt_west(
            dout=dout
            t_plt=t_plt
            )

    # Output
    return dout

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

    fig, ax = plt.subplot(2,2)

    # Plot interferomtry
    if 'interf' in dout['exp'].keys():
        ax[0,0].plot(
            dout['exp']['interf']['time_s'],
            dout['exp']['interf']['ne_m2']
            )

        ax[0,0].grid('on')
        ax[0,0].set_xlabel('time [s]')
        ax[0,0].set_ylabel('line-integrated ne [m^-2]')

    # Plots ECE
    if 'ece' in dout['exp'].keys():
        ax[0,1].plot(
            dout['exp']['ece']['time_s'],
            dout['exp']['ece']['Te_eV']
            )

        ax[0,1].grid('on')
        ax[0,1].set_xlabel('time [s]')
        ax[0,1].set_ylabel('Te [eV]')

    # Plots powers
    if 'icrf' in dout['exp'].keys():
        for ch in np.arange(dout['exp']['icrf']['time_s'].shape[1]):
            ax[1,0].plot(
                dout['exp']['icrf']['time_s'][:,ch],
                dout['exp']['icrf']['Prf_kW'][:,ch],
                label = 'IC antenna %i'%(ch)
                )

    if 'rad' in dout['exp'].keys():
        ax[1,0].plot(
            dout['exp']['rad']['time_s'],
            dout['exp']['rad']['Prad_tot_kW'],
            label = 'Radiated power'
            )

    ax[1,0].grid('on')
    ax[1,0].set_xlabel('time [s]')
    ax[1,0].set_ylabel('power [kW]')
    ax[1,0].legend()

    # Plots equilibrium
    if 'eq' in dout['exp'].keys():
        if t_plt is None:
            t_plt = int(len(dout['exp']['eq']['time_s'])/2)

        else:
            t_plt = np.argmin(abs(
                t_plt - dout['exp']['eq']['time_s']
                ))

        lvls = np.linspace(
            np.nanmin(dout['exp']['eq']['psi_vals'][t_plt, ...]),
            np.nanmax(dout['exp']['eq']['psi_vals'][t_plt, ...]),
            nlvls
            )

        ax[1,1].contour(
            dout['exp']['eq']['psi_R'],
            dout['exp']['eq']['psi_Z'],
            np.squeeze(dout['exp']['eq']['psi_vals'][t_plt, ...]),
            levels = lvls,
            colors = 'blue',
            linestyles ='-',
            linewidths = 0.5
            )

        ax[1,1].contour(
            dout['exp']['eq']['psi_R'],
            dout['exp']['eq']['psi_Z'],
            np.squeeze(dout['exp']['eq']['psi_vals'][t_plt, ...]),
            levels = (dout['exp']['eq']['psi_LCFS'][t_plt],),
            colors = 'tab:red',
            linestyles ='-',
            )

        ax[1,1].plot(
            dout['exp']['eq']['wall_R'],
            dout['exp']['eq']['wall_Z'],
            color='k'
            )

        ax[1,1].grid('on')
        ax[1,1].set_xlabel('R [m]')
        ax[1,1].set_ylabel('Z [m]')
        


##############################################
#
#               Utilities
#
##############################################

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
    data['Prf_W'] = np.zeros((ntime, nant)) # dim(ntime, nantenna)

    for aa in np.arange(nant):
        data['time_s'][:,aa] = icrf.antenna[aa].power_launched.time - dout['t_ignitron']

        data['Prf_kW'][:,aa] = icrf.antenna[aa].power_launched.data/1e3

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

    # Output
    return dout