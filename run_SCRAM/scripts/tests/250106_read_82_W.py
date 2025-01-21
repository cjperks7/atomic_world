'''

Scripts to read spectra files shared to me by S. Hansen

cjperks
jan 6, 2025

'''

# Modules
import numpy as np
import os
from atomic_world.run_SCRAM.scripts.tests import _read_sHANSEN82 as rh
import pdb
import matplotlib.pyplot as plt

# Debugging
%reload_ext autoreload
%autoreload 2

ddata = rh.main82(
    quants = ['spec', 'indspec'],
    name = 'W_ebit_',
    convert = True
    )


spec = 'spec_15.00kV'
fig, ax = plt.subplots(1,2)

ax[0].plot(
    ddata['spec']['E_eV']['data'],
    ddata['spec'][spec]['j(ph/s/cm^3/eV)'],
    color = 'k',
    label = 'total'
    )

for zz in np.arange(len(ddata['indspec'][spec]['ions'])):
    ax[0].plot(
        ddata['indspec']['E_eV']['data'],
        ddata['indspec'][spec]['j(ph/s/cm^3/eV)'][zz,:],
        label = ddata['indspec'][spec]['ions'][zz]
        )

ax[0].grid('on')
ax[0].set_xlabel('Photon energy [eV]')
ax[0].set_ylabel(r'$\varepsilon$ [$ph/s/cm^3/eV$]')
leg = ax[0].legend(labelcolor='linecolor')
leg.set_draggable('on')


for kk in ddata['spec'][spec]['cs_bal'].keys():
    if kk < ddata['indspec'][spec]['neli'] or kk > ddata['indspec'][spec]['nelf']:
        continue

    ax[1].step([kk-0.5, kk+0.5], [ddata['spec'][spec]['cs_bal'][kk]['Xz'], 0], where='post')
    ax[1].fill_between([kk-0.5,kk+0.5], [ddata['spec'][spec]['cs_bal'][kk]['Xz'], ddata['spec'][spec]['cs_bal'][kk]['Xz']])

ax[1].set_xlabel('nele')
ax[1].set_ylabel('Abundance')
ax[1].grid('on')