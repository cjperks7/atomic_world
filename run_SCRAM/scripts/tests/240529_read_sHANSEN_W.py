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
    'from_sHANSEN/240529'
    )
name = sp+'_h_scram_'

out = rh.main(
    filepath=filep,
    name=name,
    #quants = ['spec', 'table']
    )

import scipy.constants as cnt
hc = cnt.h*cnt.c/cnt.e*1e10 # [eV*AA]


#Te_plot = 4.5e3
#Te_plot = 9.5e3
#Te_plot = 2.1e3
#Te_plot = 2e4
#Te_plot = 1e3
Te_plots = [8.83e3, 9.637e3, 10.52e3, 11.48e3, 12.52e3]
#xlim = [2.70, 2.80]

#Te_plots = [
#    13.67e3,14.91e3, 16.28e3,17.76e3, 19.38e3,
 #   21.15e3, 23.08e3,25.19e3,27.49e3,30e3
 #   ]
#xlim = [0.94,0.96]
xlim = [2.71, 2.77]

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
        hc/out['spec']['lambda_A'],
        #out['spec']['lambda_A'],
        tot,
        'k--',
        label = 'total'
        )

    for ii, ion in enumerate(out['indspec']['spec_%i'%(ind)]['ions']):
        #if int(ion.split('-')[-1]) not in [13,14]:
        #    continue

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
                hc/out['indspec']['lambda_A'],
                #out['indspec']['lambda_A'],
                spec_ind,
                '-',
                label = sp+'%i+; '%(Znuc-int(ion.split('-')[-1])) +ion
                )

    Ka1 = 2.7490
    Ka2 = 2.7520

    dy = 0.0005
    aa = 0.3

    ax.axvspan(
        hc/(Ka1-dy),
        hc/(Ka1+dy),
        #(Ka1-dy),
        #(Ka1+dy),
        color = 'm',
        alpha = aa
        )
    ax.axvspan(
        hc/(Ka2-dy),
        hc/(Ka2+dy),
        #(Ka2-dy),
        #(Ka2+dy),
        color = 'm',
        alpha = aa
        )

    ax.set_xlim(hc/xlim[1], hc/xlim[0])
    #ax.set_xlim(xlim[0], xlim[1])

    ax.set_ylim(0, 1.1*np.max(tot[ss]))

    leg = ax.legend()
    leg.set_draggable('on')
    ax.grid('on')

    #ax.set_xlabel(r'$\lambda$ [$\AA$]')
    ax.set_xlabel('E_photon [eV]')
    ax.set_ylabel(r'emis [$ph/s/cm^3/\AA$]')

    fig.suptitle(
        'SCRAM (S.Hansen, May 29 2024); '
        + sp
        + '; Te= %00.2f keV'%(out['table']['Te_eV']['data'][ind]/1e3)
        )


'''
for ii, Te in enumerate(out['table']['Te_eV']['data']):
    for jj, ion in enumerate(out['indspec']['spec_%i'%(ii)]['ions']):

        with open(
            '/home/cjperks/toEBIT_240529/W_cjperks_SCRAM_spec_Te%0.1fkeV_%i.txt'%(Te/1e3,int(ion.split('-')[-1])),
            'w') as file:
            # Zip the arrays and write each pair to the file
            for item1, item2 in zip(
                hc/out['indspec']['lambda_A'], 
                (
                    out['indspec']['spec_%i'%(ii)]['j(ph/ion/s/AA)'][jj,:]
                    * out['indspec']['spec_%i'%(ii)]['abund'][jj]
                    /out['table']['ne_cm3']['data'][0]
                    * ne_cm3
                    * fi * ne_cm3
                    )
                ):
                file.write(f"{item1}\t{item2}\n")
'''