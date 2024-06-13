'''

Script to manage loading FLYCHK data

cjperks
Apr 29, 2014

'''

# Module
import numpy as np
import sys, os

from . import read_cooling_curves as rcc
from . import read_frac_abund as rfa
from . import read_ion_rates as rir

__all__ = [
    'main'
]

##########################################
#
#           Main
#
##########################################

def main(
    filepath = None,
    quants = None,
    sp = None,
    Znuc = None,
    ne = None,
    ):

    # Init
    out = {}

    if quants is None:
        quants = [
            'cooling',
            'frac_abund',
            'ion_rates',
            ]

    # Get fractional abundance at a electron density
    if 'frac_abund' in quants:
        file = os.path.join(
            filepath,
            'frac_abund',
            sp,
            sp+'_abund_Ne%0.0ecm3.dat'%(ne)
            )
        out['frac_abund'] = rfa.read(
            filepath=file
            )

    # Gets cooling curves
    if 'cooling' in quants:
        file = os.path.join(
            filepath,
            'cooling_curve',
            sp+'_cooling.dat'
            )
        out['cooling'] = rcc.read(
            filepath=file
            )

    # Gets effective ionization/recombination rates
    if 'ion_rates' in quants:
        file = os.path.join(
            filepath,
            'ion_rates'
            + '/rt.%03i.dat'%(Znuc)
            )
        out['ion_rates'] = rir.read(
            filepath=file
            )


    # Output
    return out