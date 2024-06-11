'''

Scripts to read spectra files shared to me by S. Hansen

cjperks
May 14, 2024

'''

# Modules
import numpy as np
import os
from atomic_world.run_SCRAM.scripts.tests import _read_sHANSEN as rh
import pdb

# Debugging
#%reload_ext autoreload
#%autoreload 2

# Species of interest
#sp = 'Xe'
#Znuc = 54
#fi = 1e-3
sp = 'W'
Znuc = 74
fi = 1.5e-5

ne_cm3 = 1e14

# Common path to files
filep = os.path.join(
    '/home/cjperks/work',
    '2008_SPARCxray/SCRAM',
    'from_sHANSEN/240510'
    )
name = sp+'-sparc_scram_'

out = rh.main(
    filepath=filep,
    name=name,
    #quants = ['spec', 'table']
    )





#Te_plot = 4.5e3
#Te_plot = 9.5e3
#Te_plot = 2.1e3
#Te_plot = 2e4
#Te_plot = 1e3
Te_plots = [1e3, 2.1e3, 4.5e3, 9.5e3, 2e4]
xlim = [2.70, 2.80]

for Te_plot in Te_plots:

    ind = np.argmin(abs(out['table']['Te_eV']['data'] - Te_plot))

    plt.rcParams.update({'font.size': 20})
    fig, ax = plt.subplots()

    tot = (
        out['spec']['spec_%i'%(ind)]['j(ph/ion/s/AA)']
        /out['table']['ne_cm3']['data'][0]
        * ne_cm3
        * fi * ne_cm3
        )
    tol = 1e-3

    ss = np.where(
        (out['spec']['lambda_A'] >= xlim[0])
        & (out['spec']['lambda_A'] <= xlim[1])
        )[0]

    ax.plot(
        out['spec']['lambda_A'],
        tot,
        'k--',
        label = 'total'
        )

    for ii, ion in enumerate(out['indspec']['spec_%i'%(ind)]['ions']):
        tt = np.where(
            (out['indspec']['lambda_A'] >= xlim[0])
            & (out['indspec']['lambda_A'] <= xlim[1])
            )[0]
        spec_ind = (
            out['indspec']['spec_%i'%(ind)]['j(ph/ion/s/AA)'][ii,:]
            * out['indspec']['spec_%i'%(ind)]['abund'][ii]
            /out['table']['ne_cm3']['data'][0]
            * ne_cm3
            * fi * ne_cm3
            )

        if np.max(spec_ind[tt]) >= tol*np.max(tot[ss]):
            ax.plot(
                out['indspec']['lambda_A'],
                spec_ind,
                '-',
                label = sp+'%i+; '%(Znuc-int(ion.split('-')[-1])) +ion
                )

    ax.set_xlim(xlim[0], xlim[1])


    ax.set_ylim(0, 1.1*np.max(tot[ss]))

    leg = ax.legend()
    leg.set_draggable('on')
    ax.grid('on')

    ax.set_xlabel(r'$\lambda$ [$\AA$]')
    ax.set_ylabel(r'emis [$ph/s/cm^3/\AA$]')

    fig.suptitle(
        'SCRAM (S.Hansen, May 10 2024); '
        + sp
        + '; Te= %00.2f keV'%(out['table']['Te_eV']['data'][ind]/1e3)
        )






'''

diag = 'XRSHRXe'

Te_plot = 10e3
ind = np.argmin(
    abs(spec['plasma']['Te']['data']-Te_plot)
    )

fig1, ax1 = plt.subplots()

for ion in spec[diag]['emis_line'].keys():
    for cs in spec[diag]['emis_line'][ion].keys():
        ax1.plot(
            spec[diag]['lambda']['data'],
            spec[diag]['emis_line'][ion][cs]['data'][ind,0,:],
            '-',
            label = ion+cs+'+'
            )

ax1.plot(
    spec_tot[diag]['lambda']['data'],
    spec_tot[diag]['emis']['data'][ind,0,:],
    'k--',
    label = 'total'
    )

fig1.suptitle(
    'ADAS'
    + '; Te= %00.2f keV'%(spec_tot['plasma']['Te']['data'][ind]/1e3)
    )

xlim = [2.70, 2.80]
ax1.set_xlim(xlim[0], xlim[1])

ax1.set_ylabel('emis [ph/cm3/s/AA]')
ax1.set_xlabel('lambda [AA]')

ax1.grid('on')
leg1 = ax1.legend()
leg1.set_draggable('on')

'''