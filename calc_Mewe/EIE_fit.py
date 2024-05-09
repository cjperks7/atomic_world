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
    elif nele == 2:
        return _get_He(
            sp=sp,
            Te=Te
            )
    else:
        print('NOT IMPLEMENTED YET')



##################################################
#
#            Utilities
#
##################################################

# Get He-like formulas
def _get_He(
    sp=None,
    Te=None, # [eV], dim(ntemp,)
    ):

    # Init output
    out = {}

    # Ion information
    if sp == 'Kr':
        ZZ = 36
        dE = [
            1,          # 1S0
            1.311369e4, # 1P1
            1.297883e4, # 3S1
            1.309014e4  # 3P2
            ]
        f_ij = [
            1,                   # 1S0
            6.291404e-1/(2*1+1), # 1P1, take the allowed transition
            6.291404e-1/(2*1+1), # 1P1, take the allowed transition
            6.291404e-1/(2*1+1), # 1P1, take the allowed transition
            ]

    # States considered
    out['states'] = [
        '1s^1S-2s^1S',
        '1s^1S-2p^1P',
        '1s^1S-2p^3S',
        '1s^1S-2p^3P'
        ]
    out['Egrid'] =  np.zeros((len(out['states']), int(1e3)))
    out['PEC'] = np.zeros((len(out['states']), len(Te)))
    out['cas'] = np.zeros(len(out['states']))
    out['XS'] = np.zeros((len(out['states']), int(1e3)))

    # Loop over states
    for ii, st in enumerate(out['states']):
        # Energy grid for cross-section
        # NOTE: Incident electron energy
        out['Egrid'][ii,:] = np.logspace(
            np.log10(dE[ii]),
            np.log10(dE[ii]*1e3),
            int(1e3)
            )

        # Fit coefficients
        if st == '1s^1S-2s^1S':
            if ZZ == 3:
                AA = 0.05
                BB = -0.04
                CC = 0.0
                DD = 0.0
                cas = 1.07
            elif ZZ > 3:
                AA = 0.05
                BB = 0.0
                CC = 0.0
                DD = 0.0
                cas = 1.0
        elif st == '1s^1S-2p^1P':
            if ZZ == 3:
                AA = 0.02
                BB = 0.02
                CC = 0.0
                DD = 0.28
                cas = 1.06
            elif ZZ > 3:
                AA = 0.04
                BB = 0.2
                CC = -0.04
                DD = 0.28
                cas = 1.06
        elif st == '1s^1S-2s^3S':
            if ZZ >= 3:
                AA = 0.0
                BB = 0.0
                CC = 0.07
                DD = 0.0
                cas = 1.30
        elif st == '1s^1S-2s^3P':
            if ZZ >= 3:
                AA = 0.0
                BB = 0.0
                CC = 0.2
                DD = 0.0
                cas = 1.20
        
        # Stores cascade correction
        out['cas'][ii] = cas

        (
            out['XS'][ii,:], 
            out['PEC'][ii,:], 
            ) = _calc(
            Egrid = out['Egrid'][ii,:],
            dE = dE[ii],
            AA = AA,
            BB = BB,
            CC = CC,
            DD = DD,
            f_ij = f_ij[ii],
            Te = Te
            )

    # Output
    return out

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
    elif sp == 'Kr':
        ZZ = 36
        dE = 1.342933e4
        f_ij = 2.680270e-1/(2*0.5+1)

    # Energy grid for cross-section
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

        (
            out['XS'][ii,:], 
            out['PEC'][ii,:], 
            ) = _calc(
            Egrid = out['Egrid'],
            dE = dE,
            AA = AA,
            BB = BB,
            CC = CC,
            DD = DD,
            f_ij = f_ij,
            Te = Te
            )

    # Output
    return out



def _calc(
    Egrid = None,
    dE = None,
    AA = None,
    BB = None,
    CC = None,
    DD = None,
    f_ij = None,
    Te = None,
    ):

    # Calculates Gaunt factor, Eqn 5
    UU = Egrid/dE
    gg = (
        AA
        + BB *UU**-1
        + CC *UU**-2
        + DD *np.log(UU)
        )

    # Calculates cross-section, Eqn 4
    XS = (
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
    PEC = (
        1.70e-3
        * (Te*11606)**(-1/2)
        * dE**-1
        * gbar
        * 10**(-5040*dE/Te/11606)
        * f_ij
        )

    # Output
    return XS, PEC



