'''

Script to read and process desired data
from SCRAM v8.6

'''

# Modules
import sys, os
import numpy as np
import copy
import scipy.constants as cnt

__all__ = [
    'main86',
    '_get_nele',
    '_organize_data',
    ]

#####################################################
#
#               Main
#
#####################################################

def main86(
    # File management
    fol = None,     # path to SCRAM output data
    elem = None,    # element of interest
    quants = None,  # Files to read, ['spec', 'indpsec']
    ):

    # Init
    dout = {}

    # Defaults
    if quants is None:
        quants = ['spec', 'indspec']

    # Reads total spectrum file
    if 'spec' in quants:
        dout['spec'] = _read_spec(
            file = os.path.join(
                fol,
                '%sh.spec'%(elem)
                ),
            )

    # Reads charge-state-resolved spectrum file
    if 'indspec' in quants:
        dout['indspec'] = _read_indspec(
            file = os.path.join(
                fol,
                '%sh.indspec'%(elem)
                ),
            )

    # Reads table of charge-state-resolved spectrum file
    if 'indspec_table' in quants:
        dout['indspec_table'] = _read_indspec_table(
            file = fol,
            )

    # Output
    return dout

#####################################################
#
#               Utilities
#
#####################################################

# Organizes batches of SCRAM table runs
def _organize_data_table(
    # Output
    ddata = None,
    # Controls
    fils = None,
    elem = None,
    ):

    # Init
    dorg = {}
    Znuc = _get_nele(sym=elem.split('_')[0])

    # Te grid
    Tes = []
    for ii, ff in enumerate(fils):
        for jj, kk in enumerate(ddata[ff]['indspec_table'].keys()):
            if kk == 'E':
                continue
            Tes.append(ddata[ff]['indspec_table'][kk]['Te']['data'])
    
    nTe = len(Tes)
    dorg['Te'] = {
        'data': np.sort(np.asarray(Tes)),
        'units': 'eV',
        'dim': 'dim(nTe,)',
        'name': r'$T_e$',
        'name_long': 'electron temperature'
        }

    # Init
    dorg['<Z>'] = {
        'data': np.zeros(nTe, dtype=float),
        'units': '',
        'dim': 'dim(nTe,)',
        'name': r'$\langle Z \rangle$',
        'name_long': 'effective plasma charge'
        }

    dorg['scd'] = {
        'data':np.zeros((nTe,Znuc+1), dtype=float),
        'units': r'$cm^3/s$',
        'dim': 'dim(nTe, nZ)',
        'name': 'scd',
        'name_long': 'effective ionization rate',
        }
    dorg['acd'] = {
        'data':np.zeros((nTe,Znuc+1), dtype=float),
        'units': r'$cm^3/s$',
        'dim': 'dim(nTe, nZ)',
        'name': 'acd',
        'name_long': 'effective recombination rate',
        }

    dorg['Xz'] = {
        'data':np.zeros((nTe,Znuc+1), dtype=float),
        'units': 'frac',
        'dim': 'dim(nTe, nZ)',
        'name': 'Xz',
        'name_long': 'fractional abundance',
        }

    # Gets grid data
    for ii, ff in enumerate(fils):
        for jj, kk in enumerate(ddata[ff]['indspec_table'].keys()):
            if kk == 'E':
                continue
            sdata = ddata[ff]['indspec_table'][kk]
            indT = np.argmin(abs(dorg['Te']['data'] - sdata['Te']['data']))
            ne = sdata['ne']['data']

            # <Z>
            dorg['<Z>']['data'][indT] = sdata['<Z>']['data']

            # Loop over bound electrons
            for zz, nn in enumerate(np.arange(sdata['neli'], sdata['nelf']+1)):
                dorg['scd']['data'][indT,Znuc-nn] = (
                    sdata['Rion']['data'][zz]
                    /ne
                    )
                dorg['acd']['data'][indT,Znuc-nn] = (
                    sdata['Rrec']['data'][zz]
                    /ne
                    )
                dorg['Xz']['data'][indT,Znuc-nn] = (
                    sdata['abund']['data'][zz]
                    )

    # Output
    return dorg


# Organizes a loop of SCRAM runs into a unified dictionary
def _organize_data(
    # Output
    ddata = None,
    # Controls
    lcases = None,
    elem = None,
    com = None,
    neles = None,   # List of number of bound electrons of interest
    ):

    # Init
    dorg = {}
    icase = len(lcases)

    dorg['Te'] = {
        'data': np.zeros(icase, dtype=float),
        'units': 'eV',
        'dim': 'dim(nTe,)',
        'name': r'$T_e$',
        'name_long': 'electron temperature'
        }
    dorg['<Z>'] = {
        'data': np.zeros(icase, dtype=float),
        'units': '',
        'dim': 'dim(nTe,)',
        'name': r'$\langle Z \rangle$',
        'name_long': 'effective plasma charge'
        }
    dorg['Prad'] = {
        'data': np.zeros(icase, dtype=float),
        'units': r'$W*m^3/atom/electron$',
        'dim': 'dim(nTe,)',
        'name': r'$P_{rad}$',
        'name_long': 'radiative cooling curve'
        }

    Znuc = _get_nele(sym=elem.split('_')[0])
    dorg['scd'] = {
        'data':np.zeros((icase,Znuc+1), dtype=float),
        'units': r'$cm^3/s$',
        'dim': 'dim(nTe, nZ)',
        'name': 'scd',
        'name_long': 'effective ionization rate',
        }
    dorg['acd'] = {
        'data':np.zeros((icase,Znuc+1), dtype=float),
        'units': r'$cm^3/s$',
        'dim': 'dim(nTe, nZ)',
        'name': 'acd',
        'name_long': 'effective recombination rate',
        }

    dorg['Xz'] = {
        'data':np.zeros((icase,Znuc+1), dtype=float),
        'units': 'ion/atom',
        'dim': 'dim(nTe, nZ)',
        'name': 'Xz',
        'name_long': 'fractional abundance',
        }

    nE = ddata[lcases[0]]['spec']['spec']['E(eV)']['num']
    dorg['E_photon'] = {
        'data': ddata[lcases[0]]['spec']['spec']['E(eV)']['data'],
        'units': 'eV',
        'dim': 'dim(nE,)',
        'name': r'$E_{photon}$',
        'name_long': 'Photon energy'
        }

    dorg['emis_ff'] = {
        'data':np.zeros((icase,nE), dtype=float),
        'units': r'$J*cm^3/s/eV/atom/electron$',
        'dim': 'dim(nTe, nE)',
        'name': 'ff',
        'name_long': 'free-free emissivity',
        }
    dorg['emis_fb'] = {
        'data':np.zeros((icase,nE), dtype=float),
        'units': r'$J*cm^3/s/eV/atom/electron$',
        'dim': 'dim(nTe, nE)',
        'name': 'fb',
        'name_long': 'free-bound emissivity',
        }
    dorg['emis_bb'] = {
        'data':np.zeros((icase,nE), dtype=float),
        'units': r'$J*cm^3/s/eV/atom/electron$',
        'dim': 'dim(nTe, nE)',
        'name': 'bb',
        'name_long': 'bound-bound emissivity',
        }
    if neles is not None:
        nion = len(neles)
        dorg['emis_ion'] = {
            'data':np.zeros((icase,nE, nion), dtype=float),
            'units': r'$J*cm^3/s/eV/ion/electron$',
            'dim': 'dim(nTe, nE, ncs)',
            'name': 'bb',
            'name_long': f"total emissivity per ionization charge within number of bound electrons={neles}",
            }
    dorg['emis_tot'] = {
        'data':np.zeros((icase,nE), dtype=float),
        'units': r'$J*cm^3/s/eV/atom/electron$',
        'dim': 'dim(nTe, nE)',
        'name': 'emis',
        'name_long': 'total emissivity',
        }

    # Loop over cases
    for ii in np.arange(icase):
        # Useful things
        sspec = ddata['%s%i'%(com, ii)]['spec']['spec']
        sset = ddata['%s%i'%(com, ii)]['spec']['settings']
        sbal = ddata['%s%i'%(com, ii)]['indspec']['cs_bal']
        sind = ddata['%s%i'%(com, ii)]['indspec']['spec']

        ne = sset['ne']['data']
        ni = sset['nion']['data']

        # Organizes data
        dorg['Te']['data'][ii] = sset['Te']['data']

        dorg['<Z>']['data'][ii] = sset['<Z>']['data']

        dorg['Prad']['data'][ii] = (
            sset['ploss']['data']
            *1e-7       # [erg] --> [J]
            *1e-6       # [cm^3] --> [m^3]
            /ne/ni
            )

        # Loop over bound electrons
        for jj, nn in enumerate(np.arange(sbal['neli'], sbal['nelf']+1)):
            dorg['scd']['data'][ii,Znuc-nn] = (
                sbal['Rion']['data'][jj]
                /ne
                )
            dorg['acd']['data'][ii,Znuc-nn] = (
                sbal['Rrec']['data'][jj]
                /ne
                )
            dorg['Xz']['data'][ii,Znuc-nn] = (
                sbal['abund']['data'][jj]
                )

        # Error check
        if sspec['E(eV)']['num'] != nE:
            print('Unequal spectrum grid lengths!!')
        if (
            np.max(abs(
                sspec['E(eV)']['data'] 
                - ddata[lcases[0]]['spec']['spec']['E(eV)']['data']
                )
                /ddata[lcases[0]]['spec']['spec']['E(eV)']['data']
                *100
                )
            ) > 1e-3:
            print('Changes in grid spacing!!!')

        # Organizes emissivity data
        dorg['emis_tot']['data'][ii,:] = (
            sspec['j_tot']['data']
            /ne/ni
            )
        dorg['emis_ff']['data'][ii,:] = (
            sspec['j_ff']['data']
            /ne/ni
            )
        dorg['emis_fb']['data'][ii,:] = (
            sspec['j_fb']['data']
            /ne/ni
            )
        dorg['emis_bb']['data'][ii,:] = (
            sspec['j_bb']['data']
            /ne/ni
            )

        # If b-b spectrum per ion
        if neles is not None:
            for nn, nele in enumerate(neles):
                # Find ion of interest
                if nele < sbal['neli'] or nele > sbal['nelf']:
                    continue
                
                nn_ind = list(range(sbal['neli'], sbal['nelf']+1)).index(nele)

                dorg['emis_ion']['data'][ii,:,nn] = (
                    sind['j_cs']['data'][:,nn_ind]
                    #*sbal['abund']['data'][nn_ind]
                    /ne/ni
                    )

    # Output
    return dorg


# Converts isoelectronic sequence to number of electrons
def _get_nele(
    sym = None
    ):

    labs = {
        'bare': 0,
        'H': 1,
        'He': 2,
        'Li': 3,
        'Be': 4,
        'B': 5,
        'C': 6,
        'N': 7,
        'O': 8,
        'F': 9,
        'Ne': 10,
        'Na': 11,
        'Mg': 12,
        'Al': 13,
        'Si': 14,
        'P': 15,
        'S': 16,
        'Cl': 17,
        'Ar': 18,
        'K': 19,
        'Ca': 20,
        'Sc': 21,
        'Ti': 22,
        'V': 23,
        'Cr': 24,
        'Mn': 25,
        'Fe': 26,
        'Co': 27,
        'Ni': 28,
        'Cu': 29,
        'Zn': 30,
        'Ga': 31,
        'Ge': 32,
        'As': 33,
        'Se': 34,
        'Br': 35,
        'Kr': 36,
        'Mo': 42,
        'Xe': 54,
        'W': 74,
        }

    return labs[sym]

#####################################################
#
#               File reading
#
#####################################################

# Reads the table file for indspec runs
def _read_indspec_table(
    file = None,
    ):

    # Init
    dtab = {}

    # Read file
    ff = open(file, 'r')

    ### --- Gets energy grid --- ###
    xx = ff.readline()

    # Number of grid points
    nE = int(xx.split()[0])     
    dtab['E'] = {
        'data': np.zeros(nE, dtype=float),
        'num': nE,
        'units': 'eV',
        'dim': 'dim(nE,)',
        'name': r'$E_{photon}$',
        'name_long': 'Photon energy'
        }

    # Gets energy grid
    for ii in np.arange(nE):
        xx = ff.readline()
        dtab['E']['data'][ii] = float(xx.split()[0])

    # Pad
    xx = ff.readline()

    ### --- Emissivity data --- ###

    # Get Te block
    xx = ff.readline()
    nTe = 0
    #while xx.split()[0] == 'Te(eV)':
    while len(xx.split()) > 0:
        # Init
        dtab[nTe] = {}
        stab = dtab[nTe]

        if nTe == 0:
            # Skip headers
            xx = ff.readline()
            xx = ff.readline()
            xx = ff.readline()
            xx = ff.readline()
            xx = ff.readline()

            # Header data
            xx = ff.readline()

        stab['Te'] = {
            'data': float(xx.split()[0]),
            'units': 'eV',
            'name': r'$T_e$',
            'name_long': 'electron temperature'
            }
        stab['ne'] = {
            'data': float(xx.split()[1]),
            'units': r'cm$^{-3}$',
            'name': r'$n_e$',
            'name_long': 'electron density'
            }
        stab['<Z>'] = {
            'data': float(xx.split()[2]),
            'units': r'',
            'name': r'$\langle Z\rangle$',
            'name_long': 'effective plasma charge'
            }

        # Charge state window
        xx = ff.readline()
        stab['neli'] = int(xx.split()[0])
        stab['nelf'] = int(xx.split()[1])

        xx = ff.readline()
        ncs = len(xx.split())

        # Abundances
        xx = ff.readline()
        dd = xx.split()
        stab['abund'] = {
            'data': np.asarray([float(dd[xz]) for xz in range(0, ncs)]),
            'units': 'frac'
            }

        # Ionization rates
        xx = ff.readline()
        dd = xx.split()
        stab['Rion'] = {
            'data': np.asarray([float(dd[xz]) for xz in range(0, ncs)]),
            'units': 's'
            }

        # Recombination rates
        xx = ff.readline()
        dd = xx.split()
        stab['Rrec'] = {
            'data': np.asarray([float(dd[xz]) for xz in range(0, ncs)]),
            'units': 's'
            }

        # !!! Skips Emissivity block
        for ii in np.arange(nE):
            xx = ff.readline()

        # Pad
        xx = ff.readline()

        # New Te block
        nTe += 1
        xx = ff.readline()

    # Output
    return dtab

        

# Reads the total spectrum file
def _read_spec(
    file = None,
    ):
    # Init
    dspec = {}

    # Read file
    ff = open(file, 'r')

    ### --- Settings --- ###
    key = 'settings'
    dspec[key] = {}
    tmp = dspec[key]

    # Version
    xx = ff.readline()
    tmp['version'] = xx.strip()

    # Model
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = dd[1]      # model

    # Supplementing super-configurations
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = dd[1]      # sup

    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = dd[1]      # msupbb

    # Time
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 's'}       # time

    # Electron density
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': '1/cm^3'}      # ne

    # Ion density
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': '1/cm^3'}      # nion

    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': '???'}     # den

    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': '???'}     # cfrac    

    # Electron temperature
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'eV'}     # Te

    # Ion temperature
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'eV'}     # Ti

    # Radiation temperature
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'eV'}     # Tr

    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'frac'}     # fTr  

    # Super-Gaussian flag
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = dd[1]    # mSG

    # Hot electron energy
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'eV'}     # Eh

    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'frac'}     # fH

    # Magnetic field
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'T'}     # B

    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = float(dd[1])     # Bdir

    # Controls reconstitution of averaged levels
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = float(dd[1])     # iros

    # Instrumental function
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'dE=E/(EoDE)'}     # EoDE

    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'dE=E/(EoDE)'}     # EoDEeff

    # Collisional broadening
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = float(dd[1])     # ncbr

    # Stark broadening
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = float(dd[1])     # nsbr

    # ????
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = float(dd[1])     # gsfac

    # Plasma size
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'um'}     # psize

    # Escape factors
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = float(dd[1])     # niop

    # Plasma geometry
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = float(dd[1])     # ngeo

    # Plasma velocity
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'cm/s'}     # vplas

    # Radiated power
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'erg/cm^3/s'}     # ploss

    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'eV'}     # <e>

    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {
        'bb': float(dd[1]), 
        'bf': float(dd[2]), 
        'ff': float(dd[3]), 
        'units': 'erg/cm^3/s'
        }     # bb,bf,ff contributions to ploss

    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'erg/cm^3/s'}     # non-E1

    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'erg/cm^3/s'}     # dn=0

    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'erg/cm^3/s'}     # dn=0(g)

    # ????
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {
        'data1': float(dd[1]), 
        'data2': float(dd[1]), 
        'units': '???'
        }     # sumems

    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': '???'}     # taumx

    # Average charge
    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'charge'}     # <Z>

    ### --- Charge balance --- ###
    key = 'spec'
    dspec[key] = {}
    tmp = dspec[key]
    tmp['cs_bal'] = {}

    xx = ff.readline()  # line break
    xx = ff.readline()  # header
    xx = ff.readline()  # first data line

    while len(xx.strip()) > 0:  # Unit next line break
        dd = xx.strip().split()

        nele = _get_nele(sym=dd[0])

        tmp['cs_bal'][nele] = {'Xz': float(dd[1]), 'XzgSC': float(dd[2])}
        xx = ff.readline()
        
    ### --- Emissivity spectrum --- ###
    EE = []
    tot = []
    bb = []
    bf = []
    frfr = []
    pp = []

    # Y-units
    xx = ff.readline()
    yunit = xx.strip()
    
    # Pad
    xx = ff.readline()

    # X-units
    xx = ff.readline()
    xlab = xx.strip().split()[0]
    xunit = xlab.split('(')[-1].split(')')[0]

    # Emissivity data
    xx = ff.readline()
    while len(xx.strip()) > 0:
        dd = xx.strip().split()

        # Error check
        for zz in np.arange(len(dd)):
            if not 'E' in dd[zz]:
                dd[zz] = dd[zz].replace('-', 'E-')

        EE.append(float(dd[0]))
        tot.append(float(dd[1]))
        bb.append(float(dd[2]))
        bf.append(float(dd[3]))
        frfr.append(float(dd[4]))
        pp.append(float(dd[5]))

        xx = ff.readline()

    tmp[xlab] = {'data': np.asarray(EE), 'units': xunit, 'num': len(EE)}

    tmp['j_tot'] = {'data': np.asarray(tot), 'units': yunit, 'num': len(tot)}

    tmp['j_bb'] = {'data': np.asarray(bb), 'units': yunit, 'num': len(tot)}
    tmp['j_fb'] = {'data': np.asarray(bf), 'units': yunit, 'num': len(tot)}
    tmp['j_ff'] = {'data': np.asarray(frfr), 'units': yunit, 'num': len(tot)}
    tmp['j_planck'] = {'data': np.asarray(pp), 'units': yunit, 'num': len(tot)}

    ff.close()

    # NOTE: I'm just skipping the opacities at the end lol
    return dspec

# Reads the charge-state-resolved spectrum file
def _read_indspec(
    file = None,
    ):
    # Init
    dindspec = {}

    # Read file
    ff = open(file, 'r')

    ### --- Charge balance header --- ###
    dindspec['cs_bal'] = {}
    tmp = dindspec['cs_bal']

    # Charge states simulated
    xx = ff.readline()
    dd = xx.strip().split()
    tmp['neli'] = int(dd[1])
    tmp['nelf'] = int(dd[2])

    # Abundances
    xx = ff.readline()
    dd = xx.strip().split()
    tmp['abund'] = {
        'data': np.asarray([float(dd[xz]) for xz in range(1, len(dd))]),
        'units': 'frac'
        }

    # Ionization rates
    xx = ff.readline()
    dd = xx.strip().split()
    tmp['Rion'] = {
        'data': np.asarray([float(dd[xz]) for xz in range(1, len(dd))]),
        'units': 's'
        }

    # Recombination rates
    xx = ff.readline()
    dd = xx.strip().split()
    tmp['Rrec'] = {
        'data': np.asarray([float(dd[xz]) for xz in range(1, len(dd))]),
        'units': 's'
        }

    ### --- Emissivity data --- ###
    dindspec['spec'] = {}
    tmp = dindspec['spec']

    # Y-units
    xx = ff.readline()
    yunit = xx.strip().split()[1]

    # X-units 
    xx = ff.readline()
    xlab = xx.strip().split()[0]
    xunit = xlab.split('(')[-1].split(')')[0]

    # Init
    EE = []
    sp = []

    xx = ff.readline()
    while len(xx.strip()) >0:   # Loop until line break
        dd = xx.strip().split()

        EE.append(float(dd[0]))
        sp.append(
            [float(dd[xz]) for xz in range(1, len(dd))]
            )

        xx = ff.readline()

    tmp[xlab] = {'data': np.asarray(EE), 'units': xunit, 'num': len(EE)}

    tmp['j_cs'] = {
        'data': np.asarray(sp),
        'units': yunit,
        'num': len(sp)
        }   # dim(nE, ncs)

    # Output
    return dindspec

    
