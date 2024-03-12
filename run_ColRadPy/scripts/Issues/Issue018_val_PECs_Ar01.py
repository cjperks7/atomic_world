'''

Script to compare My FAC-made PEC data files to Adam Foster's

cjperks
Feb 20, 2024


'''

# Module
import os, sys
import aurora

# Ion
sp = 'Ar'
nele = 1
Znuc = 18
ion_A = 40
fi = 1e-3

# Plasma
Te_eV = 2e3
Ti_eV = 2e3
ne_cm3 = 1e14

# Files
file_me = os.path.join(
    '/home/cjperks',
    '2008_SPARCxray/ColRadPy/FAC',
    'test_FAC_'+sp+'%02d.dat'%(nele)
    )

file_AF = os.path.join(
    '/home/cjperks',
    'tofu_sparc/atomic_data',
    'ADAS_PEC_files/fsciortino_share/atomdb',
    'pec#'+sp.lower()+'%02d.dat'%(Znuc-nele)
    )

# Broadening mechanisms
dbroad_AF_G = {
    'Voigt':{
        'Ti_eV': Ti_eV,
        'ion_A': ion_A,
        'key_options': 'wavelength',
        }
    }
dbroad_AF_V = {
    'Voigt':{
        'Ti_eV': Ti_eV,
        'ion_A': ion_A,
        'key_options': 'wavelength',
        '3.736532': 6.572170e13,       # [1/s], 2P(0.5)
        '3.731117': 6.596815e13,       # [1/s], 2P(1.5)
        }
    }

dbroad_me= {
    'Voigt':{
        'Ti_eV': Ti_eV,
        'ion_A': ion_A,
        'key_options': 'wavelength',
        '3.73653': 6.572170e13,       # [1/s], 2P(0.5)
        '3.73112': 6.596815e13,       # [1/s], 2P(1.5)
        },
    '2-photon':{
        'wavelength': [3.73594],
        'Znuc': 18,
        'nele': 1,
        },
    }
'''
wave_A = np.r_[3.736532, 3.731117]

out_me = lb.get_line_broaden(
    dbroad = dbroad_me,
    wave_A = wave_A,
    use_scipy=True,
    use_pseudo=False,
    )

out_adf = lb.get_line_broaden(
    dbroad = dbroad_AF_G,
    wave_A = wave_A,
    use_scipy=True,
    use_pseudo=False,
    )
'''
# Charge states per ionization balance
atom_data = aurora.atomic.get_atom_data(sp,['scd','acd'])
_, fz = aurora.atomic.get_frac_abundances(
    atom_data, 
    ne_cm3, Te_eV, 
    plot=False,
    )

# Aurora spectrum modeling
out_AF_G = aurora.get_local_spectrum(
    file_AF,
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
    dbroad=dbroad_AF_G
    )

out_AF_V = aurora.get_local_spectrum(
    file_AF,
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
    dbroad=dbroad_AF_V
    )

out_me = aurora.get_local_spectrum(
    file_me,
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
    dbroad=dbroad_me
    )

fig,ax = plt.subplots()

ax.plot(
    out_AF_G[0],
    out_AF_G[1],
    label = 'AF, ion, Gauss'
    )
ax.plot(
    out_AF_G[0],
    out_AF_G[2],
    label = 'AF, exc, Gauss'
    )
ax.plot(
    out_AF_G[0],
    out_AF_G[3],
    label = 'AF, RR, Gauss'
    )

ax.plot(
    out_AF_V[0],
    out_AF_V[1],
    label = 'AF, ion, Voigt'
    )
ax.plot(
    out_AF_V[0],
    out_AF_V[2],
    label = 'AF, exc, Voigt'
    )
ax.plot(
    out_AF_V[0],
    out_AF_V[3],
    label = 'AF, RR, Voigt'
    )

ax.plot(
    out_me[0],
    out_me[2],
    label = 'me, exc, Voigt'
    )

if out_me[3] is not None:
    ax.plot(
        out_me[0],
        out_me[3],
        label = 'me, rec, Voigt'
        )

ax.set_xlim(3.725, 3.745)
leg = ax.legend()
leg.set_draggable('on')
ax.grid('on')
ax.set_ylabel('emissivity [ph/s/cm3/AA]')
ax.set_xlabel('wavelength [AA]')
ax.set_title(sp+'%02d+, Te = %01.1f keV, Ti = %01.1f keV'%(Znuc-nele, Te_eV/1e3, Ti_eV/1e3))

'''
ind = np.where(
    (out_me[0] >= 3.725)
    & (out_me[0] <=3.7425)
    )[0]

np.trapz(out_me[2][ind], out_me[0][ind])

ind = np.where(
    (out_AF_V[0] >= 3.725)
    & (out_AF_V[0] <=3.7425)
    )[0]

np.trapz(out_AF_V[2][ind], out_AF_V[0][ind])

aurora.adf15_line_identification(file_me, Te_eV=3e3, ne_cm3=1e14)
aurora.adf15_line_identification(file_AF, Te_eV=3e3, ne_cm3=1e14)
'''
##############
temp_arr = np.logspace(np.log10(8.62), np.log10(8.5e4), 41) # [eV]

from colradpy import colradpy

# Ion
sp = 'Ar'
nele = 1
Znuc = 18

# Common FAC atomic data files name
file = os.path.join(
    '/home/cjperks/2008_SPARCxray',
    'FAC/ATOM',
    sp,
    sp+'%02d'%(nele)
    )

# Simulation grids
dens_arr = np.logspace(np.log10(1e10), np.log10(1e20), 21) # [cm^-3]
meta_arr = np.array([0])

# Load atomic data
crm = colradpy(
    file,
    meta_arr,
    temp_arr,
    dens_arr,
    use_recombination=True,
    use_recombination_three_body = False,
    use_ionization=True,
    suppliment_with_ecip=False,
    # New material
    atomic_data_type = 'FAC',
    ele = sp,
    nele = nele,
    Zele = Znuc,
    EEDF = 'Maxwellian',
    atomic_physics = 'incl_all',
    )


# Reads ADF15
trs_AF = aurora.read_adf15(file_AF)
trs_me = aurora.read_adf15(file_me)

# Line of interest
#lamb_AF = 3.736532
lamb_AF = 3.731117
inds = np.where(trs_AF['lambda [A]'] == lamb_AF)

ind_AF_i = inds[0][0]
ind_AF_e = inds[0][1]
ind_AF_r = inds[0][2]

#lamb_me = 3.73612
lamb_me = 3.73112
inds = np.where(trs_me['lambda [A]'] == lamb_me)

ind_me_e = inds[0][0]
ind_me_r = inds[0][1]
#ind_me_r = None

# PECs
val_AF_i = 10**trs_AF['log10 PEC fun'][ind_AF_i].ev(np.log10(1e14), np.log10(temp_arr))
val_AF_e = 10**trs_AF['log10 PEC fun'][ind_AF_e].ev(np.log10(1e14), np.log10(temp_arr))
val_AF_r = 10**trs_AF['log10 PEC fun'][ind_AF_r].ev(np.log10(1e14), np.log10(temp_arr))

val_me_e1 = 10**trs_me['log10 PEC fun'][ind_me_e].ev(np.log10(1e14), np.log10(temp_arr))
if ind_me_r is not None:
    val_me_r1 = 10**trs_me['log10 PEC fun'][ind_me_r].ev(np.log10(1e14), np.log10(temp_arr))






fig1,ax1 = plt.subplots(2,2)
fig1.tight_layout(pad=1.5)

ax1[0,0].plot(
    temp_arr/1e3,
    val_AF_i,
    label = 'AF, ioniz'
    )
ax1[0,0].plot(
    temp_arr/1e3,
    val_AF_e,
    label = 'AF, excit'
    )
ax1[0,0].plot(
    temp_arr/1e3,
    val_AF_r,
    label = 'AF, recomb'
    )

ax1[0,0].plot(
    temp_arr/1e3,
    val_me_e,
    label = 'me, excit'
    )
if ind_me_r is not None:
    ax1[0,0].plot(
        temp_arr/1e3,
        val_me_r,
        label = 'me, recomb'
        )

ax1[0,0].plot(
    temp_arr/1e3,
    crm.data['input_file']['rates']['excit']['ratec_cm3/s'][2,:],
    label = 'me, excit, rate'
    )

ax1[0,0].plot(
    temp_arr/1e3,
    crm.data['input_file']['rates']['recomb']['recomb_excit'][3,:],
    label = 'me, recomb, rate'
    )

ax1[0,0].grid('on')
ax1[0,0].set_xlabel(r'$T_e$ [$keV$]')
ax1[0,0].set_ylabel(r'PEC [$cm^3/s$]')
ax1[0,0].set_title('2P_1.5 -> 2S_0.5')
leg = ax1[0,0].legend()
leg.set_draggable('on')
ax1[0,0].set_ylim(0,5e-12)


ax1[0,1].plot(
    temp_arr/1e3,
    val_AF_i,
    label = 'AF, ioniz'
    )
ax1[0,1].plot(
    temp_arr/1e3,
    val_AF_e,
    label = 'AF, excit'
    )
ax1[0,1].plot(
    temp_arr/1e3,
    val_AF_r,
    label = 'AF, recomb'
    )

ax1[0,1].plot(
    temp_arr/1e3,
    val_me_e,
    label = 'me, excit'
    )
if ind_me_r is not None:
    ax1[0,1].plot(
        temp_arr/1e3,
        val_me_r,
        label = 'me, recomb'
        )

ax1[0,1].plot(
    temp_arr/1e3,
    crm.data['input_file']['rates']['excit']['ratec_cm3/s'][2,:],
    label = 'me, excit, rate'
    )

ax1[0,1].plot(
    temp_arr/1e3,
    crm.data['input_file']['rates']['recomb']['recomb_excit'][3,:],
    label = 'me, recomb, rate'
    )


ax1[0,1].grid('on')
ax1[0,1].set_xlabel(r'$T_e$ [$keV$]')
ax1[0,1].set_ylabel(r'PEC [$cm^3/s$]')
ax1[0,1].set_title('2P_1.5 -> 2S_0.5')

ax1[0,1].set_yscale('log')
ax1[0,1].set_xscale('log')
ax1[0,1].set_ylim(1e-20,1e-10)

print(val_me_e/crm.data['input_file']['rates']['excit']['ratec_cm3/s'][2,:])

# Line of interest
lamb_AF = 3.736532
#lamb_AF = 3.731117
inds = np.where(trs_AF['lambda [A]'] == lamb_AF)

ind_AF_i = inds[0][0]
ind_AF_e = inds[0][1]
ind_AF_r = inds[0][2]

lamb_me = 3.73653
#lamb_me = 3.73070
inds = np.where(trs_me['lambda [A]'] == lamb_me)

ind_me_e = inds[0][0]
ind_me_r = inds[0][1]
#ind_me_r = None

# PECs
val_AF_i = 10**trs_AF['log10 PEC fun'][ind_AF_i].ev(np.log10(1e14), np.log10(temp_arr))
val_AF_e = 10**trs_AF['log10 PEC fun'][ind_AF_e].ev(np.log10(1e14), np.log10(temp_arr))
val_AF_r = 10**trs_AF['log10 PEC fun'][ind_AF_r].ev(np.log10(1e14), np.log10(temp_arr))

val_me_e2 = 10**trs_me['log10 PEC fun'][ind_me_e].ev(np.log10(1e14), np.log10(temp_arr))
if ind_me_r is not None:
    val_me_r2 = 10**trs_me['log10 PEC fun'][ind_me_r].ev(np.log10(1e14), np.log10(temp_arr))

ax1[1,0].plot(
    temp_arr/1e3,
    val_AF_i,
    label = 'AF, ioniz'
    )
ax1[1,0].plot(
    temp_arr/1e3,
    val_AF_e,
    label = 'AF, excit'
    )
ax1[1,0].plot(
    temp_arr/1e3,
    val_AF_r,
    label = 'AF, recomb'
    )

ax1[1,0].plot(
    temp_arr/1e3,
    val_me_e,
    label = 'me, excit'
    )
if ind_me_r is not None:
    ax1[1,0].plot(
        temp_arr/1e3,
        val_me_r,
        label = 'me, recomb'
        )

ax1[1,0].plot(
    temp_arr/1e3,
    crm.data['input_file']['rates']['excit']['ratec_cm3/s'][0,:],
    label = 'me, excit, rate'
    )

ax1[1,0].plot(
    temp_arr/1e3,
    crm.data['input_file']['rates']['recomb']['recomb_excit'][1,:],
    label = 'me, recomb, rate'
    )



ax1[1,0].grid('on')
ax1[1,0].set_xlabel(r'$T_e$ [$keV$]')
ax1[1,0].set_ylabel(r'PEC [$cm^3/s$]')
ax1[1,0].set_title('2P_0.5 -> 2S_0.5')
leg = ax1[1,0].legend()
leg.set_draggable('on')
ax1[1,0].set_ylim(0,5e-12)


ax1[1,1].plot(
    temp_arr/1e3,
    val_AF_i,
    label = 'AF, ioniz'
    )
ax1[1,1].plot(
    temp_arr/1e3,
    val_AF_e,
    label = 'AF, excit'
    )
ax1[1,1].plot(
    temp_arr/1e3,
    val_AF_r,
    label = 'AF, recomb'
    )

ax1[1,1].plot(
    temp_arr/1e3,
    val_me_e,
    label = 'me, excit'
    )
if ind_me_r is not None:
    ax1[1,1].plot(
        temp_arr/1e3,
        val_me_r,
        label = 'me, recomb'
        )

ax1[1,1].plot(
    temp_arr/1e3,
    crm.data['input_file']['rates']['excit']['ratec_cm3/s'][0,:],
    label = 'me, excit, rate'
    )

ax1[1,1].plot(
    temp_arr/1e3,
    crm.data['input_file']['rates']['recomb']['recomb_excit'][1,:],
    label = 'me, recomb, rate'
    )



ax1[1,1].grid('on')
ax1[1,1].set_xlabel(r'$T_e$ [$keV$]')
ax1[1,1].set_ylabel(r'PEC [$cm^3/s$]')
ax1[1,1].set_title('2P_0.5 -> 2S_0.5')

ax1[1,1].set_yscale('log')
ax1[1,1].set_xscale('log')
ax1[1,1].set_ylim(1e-20,1e-10)

print(val_me_e/crm.data['input_file']['rates']['excit']['ratec_cm3/s'][0,:])