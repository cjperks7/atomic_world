'''

Functions to manage utilizing XEUS and LoWEUS data from C-Mod

cjperks
Aug 30, 2024

'''

# Modules
import MDSplus
import scipy as sp


__all__ = [
    '_get_XEUS',
    '_get_LoWEUS',
    '_get_spectrum'
    ]

####################################
#
#           Main
#
####################################

# Loads XEUS data
def _get_XEUS(
    shot = None,
    ):

    # Output
    return _get_spectrum(
        system = 'XEUS',
        shot = shot
        )

# Loads LoWEUS data
def _get_LoWEUS(
    shot = None,
    ):

    # Output
    return _get_spectrum(
        system = 'LoWEUS',
        shot = shot
        )

####################################
#
#           Utilities
#
####################################

# Loads spectrum data
def _get_spectrum(
    system = None,
    shot = None
    ):

    # Get Tree
    tt = MDSplus.Tree("spectroscopy", shot)

    # Gets Node
    NN = tt.getNode(system + '.spec')

    # Init
    dout = {}

    # Get signal
    dout['signal'] = sp.asarray(
        NN.data(), dtype=float
        ) # dim(nt, ny)

    # Gets time
    dout['time_s'] = sp.asarray(
        NN.dim_of(idx=1).data(), dtype=float
        ) # [s], dim(nt,)

    # Gets lambda
    dout['lamb_AA'] = sp.asarray(
        NN.dim_of(idx=0).data(), dtype=float
        ) # [AA], dim(ny,)

    # Output
    return dout 
