'''

This script calculates the electron impact excitation 
rate coefficient from formulas in
R. Mewe, Astronomy and Astrophysics, Vol. 20, p. 215 (1972)

Applies for H-, He-, Li-, and Ne-like

cjperks
Jan 19th, 2024

'''

# Modules
import numpy as np
import scipy.constants as cnt
from scipy.special import exp1 as E1

##################################################
#
#            Main
#
##################################################

def get_EIE(
    sp=None,
    Te=None,
    nele=None,
    ):

    # If H-like
    if nele == 1:
        return _get_H(
            sp=sp,
            Te=Te,
            )
    else:
        print('NOT IMPLEMENTED YET')



##################################################
#
#            Utilities
#
##################################################

# Get H-like formulas
def _get_H(
    sp=None,
    Te=None, # [eV], dim(ntemp,)
    ):

    # Init output
    out = {}

    # Ion information
    if sp == 'Ar':
        ZZ = 18
        dE = 3.318e3 # [eV], transition energy
        f_ij = 2.751248e-1/(2*0.5+1) # oscillator strength (from fac.TRTable)

    # Energy grid for corss-section
    # NOTE: Incident electron energy
    out['Egrid'] = np.logspace(
        np.log10(dE),
        np.log10(dE*1e3),
        int(1e3)
        )

    # States considered
    out['states'] = ['1s-2p']
    out['PEC'] = np.zeros((len(out['states']), len(Te)))
    out['cas'] = np.zeros(len(out['states']))
    out['XS'] = np.zeros((len(out['states']), len(out['Egrid'])))

    # Loop over states
    for ii, st in enumerate(out['states']):
        # Fit coefficients
        if st == '1s-2p':
            if ZZ > 2:
                AA = 0.04
                BB = 0.21
                CC = -0.04
                DD = 0.28
                cas = 1.06
            elif ZZ == 2:
                AA = 0.04
                BB = 0.04
                CC = 0.02
                DD = 0.28
                cas = 1.06

        # Stores cascade correction
        out['cas'][ii] = cas

        # Calculates Gaunt factor, Eqn 5
        UU = out['Egrid']/dE
        gg = (
            AA
            + BB *UU**-1
            + CC *UU**-2
            + DD *np.log(UU)
            )

        # Calculates cross-section, Eqn 4
        out['XS'][ii,:] = (
            4 *np.pi *cnt.physical_constants['Bohr radius'][0]**2 *1e4
            * 2 *np.pi /np.sqrt(3)
            * (13.6/dE)**2
            * f_ij
            * gg
            * UU**-1
            ) # [cm2]


        # Calculates Maxwell-averaged Gaunt factor, Eqn 8
        yy = dE/Te
        gbar = (
            AA
            + (BB*yy - CC*yy**2 +DD)
            #* (np.log((yy+1)/yy) - 0.4/(yy+1)**2)
            * np.exp(yy) *E1(yy)
            + CC*yy
            )

        # Calculates PEC, [cm3/s], Eqn 7
        out['PEC'][ii,:] = (
            1.70e-3
            * (Te*11606)**(-1/2)
            * dE**-1
            * gbar
            * 10**(-5040*dE/Te/11606)
            * f_ij
            )

    # Output
    return out



