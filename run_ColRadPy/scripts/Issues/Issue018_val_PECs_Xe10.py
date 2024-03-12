'''

Script to compare My FAC-made PEC data files to Adam Foster's

cjperks
Feb 20, 2024


'''

# Module
import os, sys
import aurora

# Ion
sp = 'Xe'
nele = 10
Znuc = 54
ion_A = 131.3
fi = 1e-3

# Plasma
Te_eV = 9e3
Ti_eV = 7e3
ne_cm3 = 1e14

# Files
file_me = os.path.join(
    '/home/cjperks',
    '2008_SPARCxray/ColRadPy/FAC',
    'test_FAC_'+sp+'%02d.dat'%(nele)
    )

file_ADAS = os.path.join(
    '/home/cjperks',
    'tofu_sparc/atomic_data',
    'ADAS_PEC_files/Xenon',
    'fs#0.50A_8.50A#'+sp.lower()+'%02d.dat'%(Znuc-nele)
    )

# Broadening mechanisms
dbroad_ADAS_G = {
    'Voigt':{
        'Ti_eV': Ti_eV,
        'ion_A': ion_A,
        'key_options': 'wavelength',
        }
    }
dbroad_ADAS_V = {
    'Voigt':{
        'Ti_eV': Ti_eV,
        'ion_A': ion_A,
        'key_options': 'wavelength',
        '2.72040': 4.619864e14,       # [1/s], 3D
        '2.72900': 2.314308e14,       # [1/s], 3F
        }
    }

# Charge states per ionization balance
atom_data = aurora.atomic.get_atom_data(sp,['scd','acd'])
_, fz = aurora.atomic.get_frac_abundances(
    atom_data, 
    ne_cm3, Te_eV, 
    plot=False,
    )

# Aurora spectrum modeling
out_ADAS_G = aurora.get_local_spectrum(
    file_ADAS,
    ne_cm3,
    Te_eV,
    ion_exc_rec_dens = [
        fi*ne_cm3*fz[0,Znuc-nele-1],
        fi*ne_cm3*fz[0,Znuc-nele],
        fi*ne_cm3*fz[0,Znuc-nele+1]
        ],
    Ti_eV=Ti_eV,
    n0_cm3 = 0.0,
    dlam_A=0.0,
    plot_spec_tot=False,
    # New Material
    dbroad=dbroad_ADAS_G
    )

out_ADAS_V = aurora.get_local_spectrum(
    file_ADAS,
    ne_cm3,
    Te_eV,
    ion_exc_rec_dens = [
        fi*ne_cm3*fz[0,Znuc-nele-1],
        fi*ne_cm3*fz[0,Znuc-nele],
        fi*ne_cm3*fz[0,Znuc-nele+1]
        ],
    Ti_eV=Ti_eV,
    n0_cm3 = 0.0,
    dlam_A=0.0,
    plot_spec_tot=False,
    # New Material
    dbroad=dbroad_ADAS_V
    )

fig,ax = plt.subplots()


ax.plot(
    out_ADAS_G[0],
    out_ADAS_G[2],
    label = 'ADAS, exc, Gauss'
    )

ax.plot(
    out_ADAS_V[0],
    out_ADAS_V[2],
    label = 'ADAS, exc, Voigt'
    )


ax.set_xlim(2.71, 2.74)
leg = ax.legend()
leg.set_draggable('on')
ax.grid('on')
ax.set_ylabel('emissivity [ph/s/cm3/AA]')
ax.set_xlabel('wavelength [AA]')
ax.set_title(sp+'%02d+, Te = %01.1f keV, Ti = %01.1f keV'%(Znuc-nele, Te_eV/1e3, Ti_eV/1e3))
