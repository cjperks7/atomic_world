'''

_main.py is a module meant to facilitate running FAC
to calculate the rates for various atomic processes as faciliated by the
"fac" module

NOTES:
    1) "pfac" -> python-version of FAC code base
    2) Based off of demo in FeL/fe.py, NeL/ne.py, https://www-amdis.iaea.org/FAC/
    3) This script is modeled off the provided pfac.atom.atomic_data() 
        specialized example

OUTPUTS:
    1) The "fac" module oututs files in the form SS##a.zz, SS##b.zz where
        i) SS is the atomic symbol
        ii) ## is the number of electrons
        iii) a -> ASCII format, b -> binary format
        iv) zz -> file type
            1) .en -> energy levels by fac.Structure
            2) .tr -> radiative transition rates by fac.TRTable
            3) .ce -> collisional excitation cross-sections by fac.CETable
            4) .ci -> collisional ionization cross-sections by fac.CITable
            5) .ai -> autoionization rates by fac.AITable
            6) .rr -. radiative recombination and photoionization cross-sections by fac.RRTable


cjperks
May 9, 2024

'''

# Modules
import numpy as np
from pfac import fac
import time
import os
from . import _en
from . import _ce_tr
from . import _ci_rr
from . import _ai_dr

__all__ = [
    'rates'
    ]

##########################################################
#
#                     Main
#
##########################################################

def rates(
    # Input
    ion=None,           # Name of ion species
    save_path=None,     # path/to/save
    # Configuration Settings
    nele=1,             # Number of electrons
    settings = None,    # Dictionary of modeling settings
    physics=None,       # If None -> include all
    ):

    # Sets the standard values for mass, charge, nuclear radius
    fac.SetAtom(ion)

    # Save file name
    name = os.path.join(
        save_path,
        '%s%02d'%(ion, nele)
        )

    # Default physics to include
    if physics is None:
        physics = ['en', 'tr', 'ce', 'rr', 'ci', 'ai']

    # Loads default settings
    if isinstance(settings, str):
        settings = _def(key=settings)

    # Calculates energy levels
    groups = _en._states(
        name=name,
        nele=nele,
        settings=settings,
        fac = fac,
        )

    # Calculates energy levels
    if 'en' in physics:
        _en._run_en(
            name=name,
            nele=nele,
            groups=groups,
            settings = settings,
            fac = fac,
            )

    # Calculates radiative transition rates
    if 'tr' in physics:
        _ce_tr._run_tr(
            name=name,
            groups=groups,
            settings = settings,
            fac = fac,
            )

    # Calculates collisional excitation cross-sections
    if 'ce' in physics:
        _ce_tr._run_ce(
            name=name,
            groups=groups,
            settings = settings,
            fac = fac,
            )

    # Calculates radiative recombination rates
    if 'rr' in physics:
        _ci_rr._run_rr(
            name=name,
            groups=groups,
            settings = settings,
            fac = fac,
            )

    # Calculates collisional ionization cross-section
    if 'ci' in physics:
        _ci_rr._run_ci(
            name=name,
            groups=groups,
            settings = settings,
            fac = fac,
            )

    # Calculates  autoionization rates
    if 'ai' in physics:
        _ai_dr._run_ai(
            name=name,
            nele=nele,
            groups=groups,
            fac = fac,
            )

    if len(physics)>0:
        # Builds an energy level table in memory from file
        fac.MemENTable(name+'b.en')

        # Converts binary file to ASCII
        #   Optional arg: v=1 -> verbose
        for ph in physics:
            fac.PrintTable(name+'b.'+ph, name+'a.'+ph, 1)

        # Clears all memory
        fac.Reinit(0)

##########################################################
#
#                     Utilites
#
##########################################################

# Setting defaults
def _def(key=None):

    # Custom settings
    if key.startswith('XEUS_Mo'):
        # Organized as "XEUS_Mo_nele"
        nouter = int(float(key.split('XEUS_Mo_')[-1])-18)
        outer = ''
        ioniz1 = ''
        ioniz2 = ''
        key = 'XEUS_Mo'
    elif key.startswith('Survey_W'):
        # Organized as "Survey_W_nele"
        nouter = int(float(key.split('Survey_W_')[-1])-10)
        if nouter < 3:
            outer = '3s'%(nouter)
        elif nouter < 9:
            outer = '3s2 3p'%(nouter-2)
        # Ionization state
        if nouter < 2:
            ioniz1 = '1*2 2*8'
        elif nouter < 4:
            ioniz1 = '1*2 2*8 3*%i;3p<1;3d<1'%(nouter-1)
        else:
            ioniz1 = '1*2 2*8 3*%i;3d<1'%(nouter-1)
        # Inner-shell ionization
        if nouter <3:
            ioniz2 = '1*2 2s2 2p5 3*%i;3p<1;3d<1'%(nouter)
        else:
            ioniz2 = '1*2 2s2 2p5 3*%i;3d<1'%(nouter)
        key = 'Survey_W'
    elif key.startswith('EBIT_W'):
        # Organized as "XEUS_Mo_nele"
        nouter = int(float(key.split('EBIT_W_')[-1])-12)
        outer = ''
        ioniz1 = ''
        ioniz2 = ''
        key = 'EBIT_W'
    else:
        nouter = 0
        outer = ''
        ioniz1 = ''
        ioniz2 = ''


    # Various default settings
    dsettings = {
        'debug':{
            'method': 'algo', # Define energy levels by my algorithm
            'grd':{                 # Number of electrons in each qunatum number at ground
                '1': [2],
                '2': [1]
                },
            'single':{
                'grd_ind': 0,       # Index of grd state base off of
                'n_raise': [2,1],   # List of electrons to excite
                'n_min': [3,2],     # Minimum quantum number to excite to
                'n_max': [3,2]      # Maximum quantum number to excite to
                },
            'double':{
                'n_raise':[],
                'n_min': [],
                'n_max': [],
                },
            'ion':{
                'n_remove': [2]   # List of quantum numbers to remove an electron
                },
            'spectator':{
                'ion_ind': 0,
                'rec_pw_max': 4,
                'n_frozen': 4,
                'n_spec_min': 4,
                'n_spec_max': 5,
                'n_ext': 45,
                }
            },

        '1':{
            'method': 'algo', # Define energy levels by my algorithm
            'grd':{                 # Number of electrons in each qunatum number at ground
                '1': [1]
                },
            'single':{              # --- Single excited electron --- # 
                'grd_ind': 0,       # Index of grd state base off of
                'n_raise': [1],   # List of electrons to excite
                'n_min': [2],     # Minimum quantum number to excite to
                'n_max': [10]      # Maximum quantum number to excite to
                },
            'double':{              # --- Two excited electrons --- #
                'n_raise': [],
                'n_min': [],
                'n_max': []
                },
            'ion':{                 # --- Ionized state ---- #
                'n_remove': [1]   # List of quantum numbers to remove an electron
                },
            'spectator':{           # --- Highly excited (spectator) settings --- #
                'ion_ind': 0,       # Ionized state to base spectators off of
                'rec_pw_max': 11,    # Maximum orbital angular momentum
                'n_frozen': 11,      # n>n_frozen uses frozen core approximation
                'n_spec_min': 11,    # Minimum quantum number of spectator
                'n_spec_max': 25,    # Maximum quantum number of spectator
                'n_ext': 45,      # Highly excited state if desired
                },
            'restrictions': None,
            },

        '2':{
            'method': 'algo', # Define energy levels by my algorithm
            'grd':{                 # Number of electrons in each quantum number at ground
                '1': [2]
                },
            'single':{              # --- Single excited electron --- # 
                'grd_ind': 0,       # Index of grd state base off of
                'n_raise': [1],   # List of electrons to excite
                'n_min': [2],     # Minimum quantum number to excite to
                'n_max': [10]      # Maximum quantum number to excite to
                },
            'double':{              # --- Two excited electrons --- #
                'n_raise': [[1,1]],
                'n_min': [[2,2]],
                'n_max': [[3,10]]
                },
            'ion':{                 # --- Ionized state ---- #
                'n_remove': [1]   # List of quantum numbers to remove an electron
                },
            'spectator':{           # --- Highly excited (spectator) settings --- #
                'ion_ind': 0,       # Ionized state to base spectators off of
                'rec_pw_max': 5,    # Maximum orbital angular momentum
                'n_frozen': 9,      # n>n_frozen uses frozen core approximation
                'n_spec_min': 10,    # Minimum quantum number of spectator
                'n_spec_max': 9,    # Maximum quantum number of spectator
                'n_ext': None,      # Highly excited state if desired
                },
            'restrictions': None,
            },


        '3':{
            'method': 'algo', # Define energy levels by my algorithm
            'grd':{                 # Number of electrons in each quantum number at ground
                '1': [2],
                '2': [1]
                },
            'single':{              # --- Single excited electron --- # 
                'grd_ind': 0,       # Index of grd state base off of
                'n_raise': [2,1],   # List of electrons to excite
                'n_min': [3,2],     # Minimum quantum number to excite to
                'n_max': [10,8]      # Maximum quantum number to excite to
                },
            'double':{              # --- Two excited electrons --- #
                'n_raise': [[1,2]],
                'n_min': [[3,3]],
                'n_max': [[3,4]]
                },
            'ion':{                 # --- Ionized state ---- #
                'n_remove': [2,1]   # List of quantum numbers to remove an electron
                },
            'spectator':{           # --- Highly excited (spectator) settings --- #
                'ion_ind': 0,       # Ionized state to base spectators off of
                'rec_pw_max': 5,    # Maximum orbital angular momentum
                'n_frozen': 9,      # n>n_frozen uses frozen core approximation
                'n_spec_min': 10,    # Minimum quantum number of spectator
                'n_spec_max': 9,    # Maximum quantum number of spectator
                'n_ext': None,      # Highly excited state if desired
                },
            'restrictions': None,
            },

        '4':{
            'method': 'algo', # Define energy levels by my algorithm
            'grd':{                 # Number of electrons in each quantum number at ground
                '1': [2],
                '2': [2]
                },
            'single':{              # --- Single excited electron --- # 
                'grd_ind': 0,       # Index of grd state base off of
                #'n_raise': [2,1],   # List of electrons to excite
                #'n_min': [3,2],     # Minimum quantum number to excite to
                #'n_max': [5,4]      # Maximum quantum number to excite to
                'n_raise': [2,1],
                'n_min': [3,2],
                'n_max': [8,3]
                },
            'double':{              # --- Two excited electrons --- #
                #'n_raise': [[1,2], [2,2]],
                #'n_min': [[3,3], [3,3]],
                #'n_max': [[3,3], [3,3]]
                'n_raise': [],
                'n_min': [],
                'n_max': []
                },
            'ion':{                 # --- Ionized state ---- #
                'n_remove': [2,1]   # List of quantum numbers to remove an electron
                },
            'spectator':{           # --- Highly excited (spectator) settings --- #
                'ion_ind': 0,       # Ionized state to base spectators off of
                'rec_pw_max': 5,    # Maximum orbital angular momentum
                'n_frozen': 7,      # n>n_frozen uses frozen core approximation
                'n_spec_min': 8,    # Minimum quantum number of spectator
                'n_spec_max': 7,    # Maximum quantum number of spectator
                'n_ext': None,      # Highly excited state if desired
                },
            'restrictions': {
                'l_max': 9
                },
            },
        
        '5':{
            'method': 'algo', # Define energy levels by my algorithm
            'grd':{                 # Number of electrons in each quantum number at ground
                '1': [2],
                '2': [3]
                },
            'single':{              # --- Single excited electron --- # 
                'grd_ind': 0,       # Index of grd state base off of
                #'n_raise': [2,1],   # List of electrons to excite
                #'n_min': [3,2],     # Minimum quantum number to excite to
                #'n_max': [5,4]      # Maximum quantum number to excite to
                'n_raise': [2,1],
                'n_min': [3,2],
                'n_max': [6,2]
                },
            'double':{              # --- Two excited electrons --- #
                #'n_raise': [[2,2]],
                #'n_min': [[3,3]],
                #'n_max': [[3,3]]
                'n_raise': [],
                'n_min': [],
                'n_max': []
                },
            'ion':{                 # --- Ionized state ---- #
                'n_remove': [2]   # List of quantum numbers to remove an electron
                },
            'spectator':{           # --- Highly excited (spectator) settings --- #
                'ion_ind': 0,       # Ionized state to base spectators off of
                'rec_pw_max': 4,    # Maximum orbital angular momentum
                'n_frozen': 6,      # n>n_frozen uses frozen core approximation
                'n_spec_min': 7,    # Minimum quantum number of spectator
                'n_spec_max': 6,    # Maximum quantum number of spectator
                'n_ext': None,      # Highly excited state if desired
                },
            'restrictions': {
                'l_max': 5
                },
            },

        '6':{
            'method': 'algo', # Define energy levels by my algorithm
            'grd':{                 # Number of electrons in each quantum number at ground
                '1': [2],
                '2': [4]
                },
            'single':{              # --- Single excited electron --- # 
                'grd_ind': 0,       # Index of grd state base off of
                #'n_raise': [2,1],   # List of electrons to excite
                #'n_min': [3,2],     # Minimum quantum number to excite to
                #'n_max': [5,4]      # Maximum quantum number to excite to
                'n_raise': [2],
                'n_min': [3],
                'n_max': [5]
                },
            'double':{              # --- Two excited electrons --- #
                #'n_raise': [[2,2]],
                #'n_min': [[3,3]],
                #'n_max': [[3,3]]
                'n_raise': [],
                'n_min': [],
                'n_max': []
                },
            'ion':{                 # --- Ionized state ---- #
                'n_remove': [2]   # List of quantum numbers to remove an electron
                },
            'spectator':{           # --- Highly excited (spectator) settings --- #
                'ion_ind': 0,       # Ionized state to base spectators off of
                'rec_pw_max': 4,    # Maximum orbital angular momentum
                'n_frozen': 6,      # n>n_frozen uses frozen core approximation
                'n_spec_min': 7,    # Minimum quantum number of spectator
                'n_spec_max': 6,    # Maximum quantum number of spectator
                'n_ext': None,      # Highly excited state if desired
                },
            'restrictions': {
                'l_max': 4
                },
            },

        '7':{
            'method': 'algo', # Define energy levels by my algorithm
            'grd':{                 # Number of electrons in each quantum number at ground
                '1': [2],
                '2': [5]
                },
            'single':{              # --- Single excited electron --- # 
                'grd_ind': 0,       # Index of grd state base off of
                #'n_raise': [2,1],   # List of electrons to excite
                #'n_min': [3,2],     # Minimum quantum number to excite to
                #'n_max': [5,4]      # Maximum quantum number to excite to
                'n_raise': [2],
                'n_min': [3],
                'n_max': [5]
                },
            'double':{              # --- Two excited electrons --- #
                #'n_raise': [[2,2]],
                #'n_min': [[3,3]],
                #'n_max': [[3,3]]
                'n_raise': [],
                'n_min': [],
                'n_max': []
                },
            'ion':{                 # --- Ionized state ---- #
                'n_remove': [2]   # List of quantum numbers to remove an electron
                },
            'spectator':{           # --- Highly excited (spectator) settings --- #
                'ion_ind': 0,       # Ionized state to base spectators off of
                'rec_pw_max': 3,    # Maximum orbital angular momentum
                'n_frozen': 5,      # n>n_frozen uses frozen core approximation
                'n_spec_min': 6,    # Minimum quantum number of spectator
                'n_spec_max': 5,    # Maximum quantum number of spectator
                'n_ext': None,      # Highly excited state if desired
                },
            'restrictions': {
                'l_max': 3
                },
            },

        '8':{
            'method': 'algo', # Define energy levels by my algorithm
            'grd':{                 # Number of electrons in each quantum number at ground
                '1': [2],
                '2': [6]
                },
            'single':{              # --- Single excited electron --- # 
                'grd_ind': 0,       # Index of grd state base off of
                #'n_raise': [2,1],   # List of electrons to excite
                #'n_min': [3,2],     # Minimum quantum number to excite to
                #'n_max': [5,4]      # Maximum quantum number to excite to
                'n_raise': [2],
                'n_min': [3],
                'n_max': [6]
                },
            'double':{              # --- Two excited electrons --- #
                #'n_raise': [[2,2]],
                #'n_min': [[3,3]],
                #'n_max': [[3,3]]
                'n_raise': [],
                'n_min': [],
                'n_max': []
                },
            'ion':{                 # --- Ionized state ---- #
                'n_remove': [2]   # List of quantum numbers to remove an electron
                },
            'spectator':{           # --- Highly excited (spectator) settings --- #
                'ion_ind': 0,       # Ionized state to base spectators off of
                'rec_pw_max': 3,    # Maximum orbital angular momentum
                'n_frozen': 7,      # n>n_frozen uses frozen core approximation
                'n_spec_min': 8,    # Minimum quantum number of spectator
                'n_spec_max': 7,    # Maximum quantum number of spectator
                'n_ext': None,      # Highly excited state if desired
                },
            'restrictions': {
                'l_max': 3
                },
            },

        '9':{
            'method': 'algo', # Define energy levels by my algorithm
            'grd':{                 # Number of electrons in each quantum number at ground
                '1': [2],
                '2': [7]
                },
            'single':{              # --- Single excited electron --- # 
                'grd_ind': 0,       # Index of grd state base off of
                #'n_raise': [2,1],   # List of electrons to excite
                #'n_min': [3,2],     # Minimum quantum number to excite to
                #'n_max': [5,4]      # Maximum quantum number to excite to
                'n_raise': [2,1],
                'n_min': [3,2],
                'n_max': [7,2]
                },
            'double':{              # --- Two excited electrons --- #
                #'n_raise': [[2,2]],
                #'n_min': [[3,3]],
                #'n_max': [[3,3]]
                'n_raise': [],
                'n_min': [],
                'n_max': []
                },
            'ion':{                 # --- Ionized state ---- #
                'n_remove': [2]   # List of quantum numbers to remove an electron
                },
            'spectator':{           # --- Highly excited (spectator) settings --- #
                'ion_ind': 0,       # Ionized state to base spectators off of
                'rec_pw_max': 3,    # Maximum orbital angular momentum
                'n_frozen': 7,      # n>n_frozen uses frozen core approximation
                'n_spec_min': 8,    # Minimum quantum number of spectator
                'n_spec_max': 7,    # Maximum quantum number of spectator
                'n_ext': None,      # Highly excited state if desired
                },
            'restrictions': {
                'l_max': 5
                },
            },

        '10':{
            'method': 'algo', # Define energy levels by my algorithm
            'grd':{                 # Number of electrons in each quantum number at ground
                '1': [2],
                '2': [8]
                },
            'single':{              # --- Single excited electron --- # 
                'grd_ind': 0,       # Index of grd state base off of
                #'n_raise': [2,1],   # List of electrons to excite
                #'n_min': [3,2],     # Minimum quantum number to excite to
                #'n_max': [5,4]      # Maximum quantum number to excite to
                'n_raise': [2,1],
                'n_min': [3,3],
                'n_max': [10,10]
                },
            'double':{              # --- Two excited electrons --- #
                #'n_raise': [[2,2]],
                #'n_min': [[3,3]],
                #'n_max': [[3,3]]
                'n_raise': [],
                'n_min': [],
                'n_max': []
                },
            'ion':{                 # --- Ionized state ---- #
                'n_remove': [2]   # List of quantum numbers to remove an electron
                },
            'spectator':{           # --- Highly excited (spectator) settings --- #
                'ion_ind': 0,       # Ionized state to base spectators off of
                'rec_pw_max': 3,    # Maximum orbital angular momentum
                'n_frozen': 7,      # n>n_frozen uses frozen core approximation
                'n_spec_min': 8,    # Minimum quantum number of spectator
                'n_spec_max': 7,    # Maximum quantum number of spectator
                'n_ext': None,      # Highly excited state if desired
                },
            'restrictions': {
                'l_max': 10
                },
            },

        '11':{
            'method': 'algo', # Define energy levels by my algorithm
            'grd':{                 # Number of electrons in each quantum number at ground
                '1': [2],
                '2': [8],
                '3': [1]
                },
            'single':{              # --- Single excited electron --- # 
                'grd_ind': 0,       # Index of grd state base off of
                #'n_raise': [2,1],   # List of electrons to excite
                #'n_min': [3,2],     # Minimum quantum number to excite to
                #'n_max': [5,4]      # Maximum quantum number to excite to
                'n_raise': [3,2],
                'n_min': [4,3],
                'n_max': [10,5]
                },
            'double':{              # --- Two excited electrons --- #
                #'n_raise': [[2,3]],
                #'n_min': [[4,5]],
                #'n_max': [[4,5]]
                'n_raise': [],
                'n_min': [],
                'n_max': []
                },
            'ion':{                 # --- Ionized state ---- #
                'n_remove': [3]   # List of quantum numbers to remove an electron
                },
            'spectator':{           # --- Highly excited (spectator) settings --- #
                'ion_ind': 0,       # Ionized state to base spectators off of
                'rec_pw_max': 3,    # Maximum orbital angular momentum
                'n_frozen': 7,      # n>n_frozen uses frozen core approximation
                'n_spec_min': 8,    # Minimum quantum number of spectator
                'n_spec_max': 7,    # Maximum quantum number of spectator
                'n_ext': None,      # Highly excited state if desired
                },
            'restrictions':{
                'l_max': 4
                },
            },

        '12':{
            'method': 'algo', # Define energy levels by my algorithm
            'grd':{                 # Number of electrons in each quantum number at ground
                '1': [2],
                '2': [8],
                '3': [2]
                },
            'single':{              # --- Single excited electron --- # 
                'grd_ind': 0,       # Index of grd state base off of
                #'n_raise': [3,2],   # List of electrons to excite
                #'n_min': [4,3],     # Minimum quantum number to excite to
                #'n_max': [7,4]      # Maximum quantum number to excite to
                'n_raise': [2],
                'n_min': [3],
                'n_max': [3]
                },
            'double':{              # --- Two excited electrons --- #
                #'n_raise': [[2,2]],
                #'n_min': [[3,3]],
                #'n_max': [[3,3]]
                'n_raise': [],
                'n_min': [],
                'n_max': []
                },
            'ion':{                 # --- Ionized state ---- #
                'n_remove': [3]   # List of quantum numbers to remove an electron
                },
            'spectator':{           # --- Highly excited (spectator) settings --- #
                'ion_ind': 0,       # Ionized state to base spectators off of
                'rec_pw_max': 3,    # Maximum orbital angular momentum
                'n_frozen': 7,      # n>n_frozen uses frozen core approximation
                'n_spec_min': 8,    # Minimum quantum number of spectator
                'n_spec_max': 7,    # Maximum quantum number of spectator
                'n_ext': None,      # Highly excited state if desired
                },
            'restrictions': {
                'l_max': 3,
                'skip_CI': True,    # Skip configuration interactions (ce: exc->exc)
                },
            },

        '13':{
            'method': 'algo', # Define energy levels by my algorithm
            'grd':{                 # Number of electrons in each quantum number at ground
                '1': [2],
                '2': [8],
                '3': [3]
                },
            'single':{              # --- Single excited electron --- # 
                'grd_ind': 0,       # Index of grd state base off of
                #'n_raise': [2,1],   # List of electrons to excite
                #'n_min': [3,2],     # Minimum quantum number to excite to
                #'n_max': [5,4]      # Maximum quantum number to excite to
                'n_raise': [2],
                'n_min': [3],
                'n_max': [3]
                },
            'double':{              # --- Two excited electrons --- #
                #'n_raise': [[2,2]],
                #'n_min': [[3,3]],
                #'n_max': [[3,3]]
                'n_raise': [],
                'n_min': [],
                'n_max': []
                },
            'ion':{                 # --- Ionized state ---- #
                'n_remove': [3]   # List of quantum numbers to remove an electron
                },
            'spectator':{           # --- Highly excited (spectator) settings --- #
                'ion_ind': 0,       # Ionized state to base spectators off of
                'rec_pw_max': 3,    # Maximum orbital angular momentum
                'n_frozen': 7,      # n>n_frozen uses frozen core approximation
                'n_spec_min': 8,    # Minimum quantum number of spectator
                'n_spec_max': 7,    # Maximum quantum number of spectator
                'n_ext': None,      # Highly excited state if desired
                },
            'restrictions': {
                'l_max': 3,
                'skip_CI': True,    # Skip configuration interactions (ce: exc->exc)
                },
            },

        '14':{
            'method': 'algo', # Define energy levels by my algorithm
            'grd':{                 # Number of electrons in each quantum number at ground
                '1': [2],
                '2': [8],
                '3': [4]
                },
            'single':{              # --- Single excited electron --- # 
                'grd_ind': 0,       # Index of grd state base off of
                #'n_raise': [2,1],   # List of electrons to excite
                #'n_min': [3,2],     # Minimum quantum number to excite to
                #'n_max': [5,4]      # Maximum quantum number to excite to
                'n_raise': [2],
                'n_min': [3],
                'n_max': [3]
                },
            'double':{              # --- Two excited electrons --- #
                #'n_raise': [[2,2]],
                #'n_min': [[3,3]],
                #'n_max': [[3,3]]
                'n_raise': [],
                'n_min': [],
                'n_max': []
                },
            'ion':{                 # --- Ionized state ---- #
                'n_remove': [3]   # List of quantum numbers to remove an electron
                },
            'spectator':{           # --- Highly excited (spectator) settings --- #
                'ion_ind': 0,       # Ionized state to base spectators off of
                'rec_pw_max': 3,    # Maximum orbital angular momentum
                'n_frozen': 7,      # n>n_frozen uses frozen core approximation
                'n_spec_min': 8,    # Minimum quantum number of spectator
                'n_spec_max': 7,    # Maximum quantum number of spectator
                'n_ext': None,      # Highly excited state if desired
                },
            'restrictions': {
                'l_max': 3,
                'skip_CI': True,    # Skip configuration interactions (ce: exc->exc)
                },
            },

        '15':{
            'method': 'algo', # Define energy levels by my algorithm
            'grd':{                 # Number of electrons in each quantum number at ground
                '1': [2],
                '2': [8],
                '3': [5]
                },
            'single':{              # --- Single excited electron --- # 
                'grd_ind': 0,       # Index of grd state base off of
                #'n_raise': [2,1],   # List of electrons to excite
                #'n_min': [3,2],     # Minimum quantum number to excite to
                #'n_max': [5,4]      # Maximum quantum number to excite to
                'n_raise': [2],
                'n_min': [3],
                'n_max': [3]
                },
            'double':{              # --- Two excited electrons --- #
                #'n_raise': [[2,2]],
                #'n_min': [[3,3]],
                #'n_max': [[3,3]]
                'n_raise': [],
                'n_min': [],
                'n_max': []
                },
            'ion':{                 # --- Ionized state ---- #
                'n_remove': [3]   # List of quantum numbers to remove an electron
                },
            'spectator':{           # --- Highly excited (spectator) settings --- #
                'ion_ind': 0,       # Ionized state to base spectators off of
                'rec_pw_max': 3,    # Maximum orbital angular momentum
                'n_frozen': 7,      # n>n_frozen uses frozen core approximation
                'n_spec_min': 8,    # Minimum quantum number of spectator
                'n_spec_max': 7,    # Maximum quantum number of spectator
                'n_ext': None,      # Highly excited state if desired
                },
            'restrictions': {
                'l_max': 3,
                'skip_CI': True,    # Skip configuration interactions (ce: exc->exc)
                },
            },

        'XEUS_Mo':{
            'method': 'string', # Define energy levels by a direct string
            'grd':{                 # Number of electrons in each quantum number at ground
                'grd.0': '1*2 2*8 3s2 3p6 3d%i'%(nouter),
                },
            'exc':{              # --- Excited state --- # 
                'exc.3d->4f': '1*2 2*8 3s2 3p6 3d%i 4f1'%(nouter-1),
                },
            'ion':{                 # --- Ionized state ---- #
                'ion.0': '1*2 2*8 3s2 3p6 3d%i'%(nouter-1),   # List of quantum numbers to remove an electron
                },
            'spectator':{           # --- Highly excited (spectator) settings --- #
                'ion_ind': 0,       # Ionized state to base spectators off of
                'rec_pw_max': 3,    # Maximum orbital angular momentum
                'n_frozen': 7,      # n>n_frozen uses frozen core approximation
                'n_spec_min': 8,    # Minimum quantum number of spectator
                'n_spec_max': 7,    # Maximum quantum number of spectator
                'n_ext': None,      # Highly excited state if desired
                },
            'restrictions': {
                'skip_CI': True,    # Skip configuration interactions (ce: exc->exc)
                },
            },

        'Survey_W':{
            'method': 'string', # Define energy levels by a direct string
            'grd':{                 # Number of electrons in each quantum number at ground
                'grd.0': '1*2 2*8 '+outer,
                },
            'exc':{             # --- Excited state --- # 
                'exc.3->3': '1*2 2*8 3*%i;3d<2'%(nouter),
                'exc.2p->3': '1*2 2s2 2p5 3*%i;3d<2'%(nouter+1),
                },
            'ion':{                 # --- Ionized state ---- #
                'ion.0': ioniz1,   # List of quantum numbers to remove an electron
                'ion.1': ioniz2,
                },
            'spectator':{           # --- Highly excited (spectator) settings --- #
                'ion_ind': 0,       # Ionized state to base spectators off of
                'rec_pw_max': 3,    # Maximum orbital angular momentum
                'n_frozen': 7,      # n>n_frozen uses frozen core approximation
                'n_spec_min': 8,    # Minimum quantum number of spectator
                'n_spec_max': 7,    # Maximum quantum number of spectator
                'n_ext': None,      # Highly excited state if desired
                },
            'restrictions': {
                'skip_CI': False,    # Skip configuration interactions (ce: exc->exc)
                },
            },

        'EBIT_W':{
            'method': 'string', # Define energy levels by a direct string
            'grd':{                 # Number of electrons in each quantum number at ground
                'grd.0': '1*2 2*8 3s2 3p%i'%(nouter),
                'grd.1': '1*2 2*8 3*%i;3p<%i;3d>0;3d<2'%(nouter+2, nouter),
                },
            'exc':{             # --- Excited state --- # 
                #'exc.3->4f': '1*2 2*8 3*%i 4f1;3p<%i;3d<1'%(nouter-1+2, nouter),
                #'exc.3->5f': '1*2 2*8 3*%i 5f1;3p<%i;3d<1'%(nouter-1+2, nouter),
                'exc.3->6f': '1*2 2*8 3*%i 6f1;3p<%i;3d<1'%(nouter-1+2, nouter),
                #'exc.3->3d': '1*2 2*8 3*%i;3p<%i;3d>0;3d<2'%(nouter+2, nouter),
                #'exc.3->4d': '1*2 2*8 3*%i 4d1;3p<%i;3d<1'%(nouter+1, nouter),
                'exc.3->5d': '1*2 2*8 3*%i 5d1;3p<%i;3d<1'%(nouter+1, nouter),
                },
            'ion':{                 # --- Ionized state ---- #
                'ion.0': '1*2 2*8 3*%i;3d<1'%(nouter+2-1)   # List of quantum numbers to remove an electron
                },
            'spectator':{           # --- Highly excited (spectator) settings --- #
                'ion_ind': 0,       # Ionized state to base spectators off of
                'rec_pw_max': 3,    # Maximum orbital angular momentum
                'n_frozen': 7,      # n>n_frozen uses frozen core approximation
                'n_spec_min': 8,    # Minimum quantum number of spectator
                'n_spec_max': 7,    # Maximum quantum number of spectator
                'n_ext': None,      # Highly excited state if desired
                },
            'restrictions': {
                'skip_CI': True,    # Skip configuration interactions (ce: exc->exc)
                },
            },

        }

    # Output
    return dsettings[key]

