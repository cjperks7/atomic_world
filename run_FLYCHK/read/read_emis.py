'''

Script to read emissivity files simulated by FLYCHK

cjperks
May 20, 2024

'''

# Modules
import numpy as np
import os, sys
from aurora import line_broaden as lb
from scipy.interpolate import interp1d
import scipy.constants as cnt
import matplotlib.pyplot as plt

__all__ = [
    'read',
    'plot_spec'
    ]

##########################################
#
#               Main
#
##########################################

# Main script to read data files
def read(
    path = '/home/cjperks/work/2008_SPARCxray/FLYCHK_NLTE/emis',
    suffix = None,
    sp = None,
    Te_eV = None,
    ne_cm3 = None,
    Ti_eV = None,
    nimp_cm3 = None,
    convolve = True,
    ):

    # File name
    file = os.path.join(
        path,
        sp,
        sp+'_Te%0.2fkeV_ne%0.2ecm3'%(Te_eV/1e3, ne_cm3)
        )
    if suffix is not None:
        file += suffix+'.dat'
    else:
        file += '.dat'

    # Loads spectral emissivity data
    if os.path.isfile(file):
        out = _read_data(
            file = file,
            sp = sp,
            nimp_cm3 = nimp_cm3
            )
    else:
        print('FILE NOT FOUND')
        print(file)
        sys.exit(1)

    # If convolve with a Voigt profile
    if convolve:
        if Ti_eV is None:
            Ti_eV = Te_eV
        out['Ti_eV'] = Ti_eV

        out = _convolve(
            out = out,
            sp = sp,
            Ti_eV = Ti_eV
            )

    # Output
    return out


##########################################
#
#               Utilities
#
##########################################

# Convolves emissivity data with Voigt profile
def _convolve(
    out = None,
    sp = None,
    Ti_eV = None,
    ):

    # Ion mass
    if sp == 'W':
        ion_A = 183.8


    # Loop over charge states
    for cs in out[sp].keys():
        # Broadening mechanism
        dbroad = {
            'Voigt':{
                'Ti_eV':Ti_eV,
                'ion_A':ion_A,
                'key_options':'isel'
                }
            }

        for ii in range(len(out[sp][cs]['Arates_1/s'])):
            dbroad['Voigt'][str(ii+1)] = out[sp][cs]['Arates_1/s'][ii]

        # Calculates line shapes, dim(ntrs,nlambda), [AA], [1/AA]
        lams_profs_A, theta = lb.get_line_broaden(
            dbroad=dbroad,
            wave_A = out[sp][cs]['lambda_A']
            )

        # non-equally spaced wavelenght
        wave_final_A = np.unique(lams_profs_A)

        # Calculates spectrum
        out[sp][cs]['spec_ph/cm3/s/AA'] = np.zeros_like(wave_final_A)

        # Adds line shape
        for ii in np.arange(lams_profs_A.shape[0]):
            out[sp][cs]['spec_ph/cm3/s/AA'] += (
                out[sp][cs]['emis_ph/cm3/s'][ii]
                * interp1d(
                    lams_profs_A[ii,:],
                    theta[ii,:],
                    bounds_error=False,
                    fill_value = 0.0
                    )(wave_final_A)
                )
        out[sp][cs]['spec_lambda_A'] = wave_final_A

    # Output
    return out

# Reads emissivity data from FLYCHK
def _read_data(
    file = None,
    sp = None,
    nimp_cm3 = None
    ):

    # Init
    out = {}

    # Opens data file
    f = open(file, 'r')

    # Reads header data
    xx = f.readline()
    out['version'] = xx.split(' ')[-2]

    xx = f.readline()
    out['Znuc'] = int(xx.split(' ')[-1])

    xx = f.readline()
    out['Te_eV'] = float(xx.split(' ')[-1])

    xx = f.readline()
    out['ne_cm3'] = float(xx.split(' ')[-1])

    xx = f.readline()
    if nimp_cm3 is None:
        out['nimp_cm3'] = float(xx.split(' ')[-1])
        rescale_nimp = 1
    else:
        out['nimp_cm3'] = nimp_cm3
        rescale_nimp = nimp_cm3/float(xx.split(' ')[-1])

    # Flags
    take_data = False

    out[sp] = {}

    # Loop over file
    for line in f:
        if take_data:
            # Makes entry for charge state
            cs = int(line.split()[0])
            if cs not in out[sp].keys():
                out[sp][cs] = {}
                out[sp][cs]['Energy_eV'] = []
                out[sp][cs]['lambda_A'] = []
                out[sp][cs]['emis_ph/cm3/s'] = []
                out[sp][cs]['Arates_1/s'] = [] 
                out[sp][cs]['upper'] = []
                out[sp][cs]['lower'] = []

            # Save wavelength data
            out[sp][cs]['Energy_eV'].append(
                float(line.split()[1])
                )
            out[sp][cs]['lambda_A'].append(
                float(line.split()[2])
                )

            # Saves emissivity data
            out[sp][cs]['emis_ph/cm3/s'].append(
                float(line.split()[4])
                * 4 * np.pi
                * rescale_nimp
                *6.242e11/out[sp][cs]['Energy_eV'][-1]
                )

            # Saves A-value
            out[sp][cs]['Arates_1/s'].append(
                float(line.split()[6])
                )

            # Saves level configs
            out[sp][cs]['lower'].append(
                line.split()[-2]
                )
            out[sp][cs]['upper'].append(
                line.split()[-1]
                )

        elif line.split(' ')[0] == 'Iso':
            take_data = True

    # Reformatting
    for cs in out[sp].keys():
        out[sp][cs]['Energy_eV'] = np.asarray(out[sp][cs]['Energy_eV'])
        out[sp][cs]['lambda_A'] = np.asarray(out[sp][cs]['lambda_A'])
        out[sp][cs]['emis_ph/cm3/s'] = np.asarray(out[sp][cs]['emis_ph/cm3/s'])
        out[sp][cs]['Arates_1/s'] = np.asarray(out[sp][cs]['Arates_1/s'])

    # Output
    return out

##########################################
#
#               Plotting
#
##########################################

def plot_spec(
    out=None,
    sp = None,
    x_units ='AA',
    y_units = 'ph/cm3/s/AA',
    xlim = np.r_[2.7,2.8], # [AA]
    incl_cs = None,
    tol = 1e-2,
    lw = 2,
    calib_vals = None,
    calib_wdth = 0.0005,
    ):

    # Useful constants
    hc = cnt.h*cnt.c/cnt.e*1e10 # [eV*AA]

    # Total spectrum
    xx_tot = np.linspace(xlim[0], xlim[1], 1000)
    spec_tot = np.zeros_like(xx_tot)

    for cs in out[sp].keys():
        spec_tot += interp1d(
            out[sp][cs]['spec_lambda_A'],
            out[sp][cs]['spec_ph/cm3/s/AA'],
            bounds_error = False,
            fill_value = 0.0
            )(xx_tot)


    plt.rcParams.update({'font.size': 20})
    fig,ax = plt.subplots()

    if x_units == 'AA':
        xx_tmp = xx_tot
        xxlabel = r'$\lambda$, wavelength [$\AA$]'
        xlim_og = xlim.copy()
    elif x_units == 'eV':
        xx_tmp = hc/xx_tot
        xlim_og = xlim.copy()
        xlim = hc/xlim
        xxlabel = 'E_photon [eV]'

    if y_units == 'ph/cm3/s/AA':
        yy_tmp = spec_tot
        yy_max = np.max(spec_tot)
        yylabel = r'$\epsilon$, emissivity [$ph/s/cm^3/\AA$]'

    ax.plot(
        xx_tmp,
        yy_tmp,
        'k-',
        linewidth = lw,
        label = 'total'
        )

    # Plots charge-state resolved spectra
    for cs in out[sp].keys():
        if incl_cs is not None:
            if cs not in incl_cs:
                continue

        if x_units == 'AA':
            xx_tmp = out[sp][cs]['spec_lambda_A']
        elif x_units == 'eV':
            xx_tmp = hc/out[sp][cs]['spec_lambda_A']

        if y_units == 'ph/cm3/s/AA':
            yy_tmp = out[sp][cs]['spec_ph/cm3/s/AA']

        ind = np.where(
            (out[sp][cs]['spec_lambda_A'] >= xlim_og[0])
            & (out[sp][cs]['spec_lambda_A'] <= xlim_og[1])
            )[0]
                
        if np.max(yy_tmp[ind]) >= tol*yy_max:
            ax.plot(
                xx_tmp,
                yy_tmp,
                label = sp + '%i+, %i'%(out['Znuc']-cs,cs),
                linewidth = lw,
                linestyle = 'solid',
                )

    # Plots extra for calibration sources
    if calib_vals is not None:
        for ii, cc in enumerate(calib_vals):

            if x_units == 'AA':
                xx1 = cc - calib_wdth
                xx2 = cc + calib_wdth
            elif x_units == 'eV':
                xx1 = hc/(cc - calib_wdth)
                xx2 = hc/(cc + calib_wdth)

            ax.axvspan(
                xx1,
                xx2,
                color = 'm',
                alpha = 0.3
                )

    ax.set_xlim(np.sort(xlim))
    ax.set_ylim(0, 1.1*yy_max)

    leg = ax.legend()
    leg.set_draggable('on')
    ax.grid('on')

    ax.set_ylabel(yylabel)
    ax.set_xlabel(xxlabel)

    ax.set_title(
        r'FLYCHK spectrum, $T_e$ = %01.2f keV, $T_i$ = %01.2f keV'%(
            out['Te_eV']/1e3, out['Ti_eV']/1e3
            )
        )


    

