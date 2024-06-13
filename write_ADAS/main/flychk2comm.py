'''

Script to write FLYCHK data to a common data format

'''

# Modules
from atomic_world.run_FLYCHK import read as FLY
from .comm2adas import _get_atom
import numpy as np

#######################################
#
#             Main
#
#######################################

# Converts FLYCHK data
def convert(
    filepath=None,
    filetype=None,
    sp=None,
    ):

    # If converting factional abundance data
    if filetype == 'ion_rates':
        return _convert_ion_rates(
            filepath=filepath,
            sp=sp
            )

#######################################
#
#             Utilities
#
#######################################

# Converts FLYCHK abundance rates to common format
def _convert_ion_rates(
    filepath=None,
    sp=None,
    ):

    # Atom information
    datom = _get_atom(sym=sp)

    # Obtains rate data
    #### NOTE: FLYCHK data posted on:
    #### https://www-amdis.iaea.org/FLYCHK/
    #### is calculated at Ne=1[#/cm3]
    data = FLY.main(
        filepath=filepath,
        quants = 'ion_rates',
        sp=sp,
        Znuc = datom['Zimp']
        )

    # Formats data
    Te_eV = data['ion_rates']['ci']['Te_eV'] # [eV], dim(Te,)

    scd = (
        data['ion_rates']['ci']['data']
        +data['ion_rates']['ea']['data']
        ) # [cm3/s], dim(Te,cs)
    acd = (
        data['ion_rates']['rr']['data']
        +data['ion_rates']['dr']['data']
        ) # [cm3/s] dim(Te,cs)

    # Dummy values for density
    ne_cm3 = np.r_[1e10, 1e15] # [cm^-3]

    # Accounts for density shape
    scd = np.repeat(
        scd[None,:,:], 
        axis=0,
        repeats=len(ne_cm3)
        ) # dim(ne,Te,cs)
    acd = np.repeat(
        acd[None,:,:], 
        axis=0,
        repeats=len(ne_cm3)
        ) # dim(ne,Te,cs)

    # Ouput
    return Te_eV, ne_cm3, acd, scd


