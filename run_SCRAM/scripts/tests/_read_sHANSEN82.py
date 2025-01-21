'''

Reads data files produced from SCRAM v8.2

cjperks
Jan 6th, 2025

'''

# Modules
import sys, os
import numpy as np
import copy
import scipy.constants as cnt

__all__ = [
    'main82',
    'read_linelist82',
    '_get_nele',
    ]

#####################################################
#
#               Main
#
#####################################################

def main82(
    filepath = None,
    name = None,
    quants = None,
    convert = True,
    ):

    # Init
    out = {}

    # Defaults
    if quants is None:
        quants = ['indspec', 'spec', 'table']
    if filepath is None:
        filepath = os.path.join(
            '/home/cjperks/work',
            '2008_SPARCxray/SCRAM',
            'from_sHANSEN/241220'
            )

    names = [f.split('.spec')[0] for f in os.listdir(filepath) if f.startswith(name) and f.endswith('.spec')]
    EEs = []

    for nn in names:
        print(nn)
        Ebeam = float(nn.split('_')[-1]) # [kV]
        EEs.append(Ebeam)

        # Reads table file
        if 'table' in quants:
            out = _read_table(
                out=out,
                file = os.path.join(
                    filepath,
                    nn+'.table'
                    ),
                Ebeam = Ebeam,
                )

        # Reads total spectrum file
        if 'spec' in quants:
            out = _read_spec(
                out=out,
                file = os.path.join(
                    filepath,
                    nn+'.spec'
                    ),
                Ebeam = Ebeam,
                )

        # Reads ion-resolved spectrum file
        if 'indspec' in quants:
            out = _read_indspec(
                out=out,
                file = os.path.join(
                    filepath,
                    nn+'.indspec'
                    ),
                Ebeam = Ebeam,
                )

    # Not checking it, but I know all the energy grids are the same
    out['spec']['E_eV'] = out['indspec']['E_eV'] = out['spec']['spec_%0.2fkV'%(Ebeam)]['E_eV']
    out['settings']['Ebeam'] = np.asarray(EEs)

    # Converts units to what I like
    if convert:
        out = _convert(out=out)

    # Output
    return out


# Reads the linelist
def read_linelist82(
    file=None,
    # Filtering
    x_units = 'eV', # 'eV' or 'AA'
    x_lim = [1,2],
    ):

    # Init
    dout = {}

    ion = ''

    upr = []
    upr_qual = []
    lwr = []
    lwr_qual = []
    g_upr = []
    g_lwr = []
    dE = []
    wl = []
    ems = []
    Arad = []
    Aug = []
    flx_ex = []
    flx_ion = []

    # Reads file
    first = True
    with open(file, 'r') as ff:
        for line in ff:
            # Skips the first line
            if first:
                first = False
                continue

            # To deal with memory, we cut up the linelist processing by charge state
            dd = line.strip().split()

            # New ion so process the last data
            if dd[0] != ion:
                ion_old = copy.copy(ion)
                ion = copy.copy(dd[0])

                # Processes data
                if len(upr) >0:
                    dout[ion_old] = {}

                    dout[ion_old]['upr'] = {'data': np.asarray(upr), 'units': ''}
                    dout[ion_old]['upr_qual'] = {'data': np.asarray(upr_qual), 'units': ''}
                    dout[ion_old]['lwr'] = {'data': np.asarray(lwr), 'units': ''}
                    dout[ion_old]['lwr_qual'] = {'data': np.asarray(lwr_qual), 'units': ''}
                    dout[ion_old]['g_upr'] = {'data': np.asarray(g_upr), 'units': ''}
                    dout[ion_old]['g_lwr'] = {'data': np.asarray(g_lwr), 'units': ''}
                    dout[ion_old]['dE'] = {'data': np.asarray(dE), 'units': 'eV'}
                    dout[ion_old]['wl'] = {'data': np.asarray(wl), 'units': 'AA'}
                    dout[ion_old]['ems'] = {'data': np.asarray(ems), 'units': 'ph/s/atom'}
                    dout[ion_old]['Arad'] = {'data': np.asarray(Arad), 'units': '1/s'}
                    dout[ion_old]['Aug'] = {'data': np.asarray(Aug), 'units': '1/s'}
                    dout[ion_old]['flx_ex'] = {'data': np.asarray(flx_ex), 'units': 'frac'}
                    dout[ion_old]['flx_ion'] = {'data': np.asarray(flx_ion), 'units': 'frac'}

                    # Clear data for next
                    upr = []
                    upr_qual = []
                    lwr = []
                    lwr_qual = []
                    g_upr = []
                    g_lwr = []
                    dE = []
                    wl = []
                    ems = []
                    Arad = []
                    Aug = []
                    flx_ex = []
                    flx_ion = []

            # Error check
            for zz in range(5, 20+1):
                if not 'E' in dd[zz]:
                    dd[zz] = dd[zz].replace('-', 'E-')


            # Filtering
            if x_units == 'eV':
                if float(dd[7]) < x_lim[0] or float(dd[7]) > x_lim[1]:
                    continue
            elif x_units == 'AA':
                if float(dd[8]) < x_lim[0] or float(dd[8]) > x_lim[1]:
                    continue

            # Gets data
            upr.append(dd[1])
            upr_qual.append(dd[2])
            lwr.append(dd[3])
            lwr_qual.append(dd[4])
            g_upr.append(float(dd[5]))
            g_lwr.append(float(dd[6]))
            dE.append(float(dd[7]))
            wl.append(float(dd[8]))
            ems.append(float(dd[9]))
            Arad.append(float(dd[10]))
            Aug.append(float(dd[11]))
            flx_ex.append(float(dd[19]))
            flx_ion.append(float(dd[20]))

    # Saves last ion
    dout[ion] = {}
    dout[ion]['upr'] = {'data': np.asarray(upr), 'units': ''}
    dout[ion]['upr_qual'] = {'data': np.asarray(upr_qual), 'units': ''}
    dout[ion]['lwr'] = {'data': np.asarray(lwr), 'units': ''}
    dout[ion]['lwr_qual'] = {'data': np.asarray(lwr_qual), 'units': ''}
    dout[ion]['g_upr'] = {'data': np.asarray(g_upr), 'units': ''}
    dout[ion]['g_lwr'] = {'data': np.asarray(g_lwr), 'units': ''}
    dout[ion]['dE'] = {'data': np.asarray(dE), 'units': 'eV'}
    dout[ion]['wl'] = {'data': np.asarray(wl), 'units': 'AA'}
    dout[ion]['ems'] = {'data': np.asarray(ems), 'units': 'ph/s/atom'}
    dout[ion]['Arad'] = {'data': np.asarray(Arad), 'units': '1/s'}
    dout[ion]['Aug'] = {'data': np.asarray(Aug), 'units': '1/s'}
    dout[ion]['flx_ex'] = {'data': np.asarray(flx_ex), 'units': 'frac'}
    dout[ion]['flx_ion'] = {'data': np.asarray(flx_ion), 'units': 'frac'}

    # Output
    return dout


#####################################################
#
#               Utilities
#
#####################################################

# Converts isoelectronic sequence to number of electrons
def _get_nele(
    sym = None
    ):

    labs = {
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
        }

    return labs[sym]

# Converts emissivty units from W/cm^3/ion/eV -> ph/ion/s/cm^3/AA
def _convert(
    out=None
    ):

    # Useful values
    hc = cnt.h*cnt.c/cnt.e*1e10 # [eV*AA]
    keys = list(out.keys())

    # Loop over spectra
    for key in keys:
        if key in ['spec', 'indspec']:
            out[key]['lambda_A'] = hc/out[key]['E_eV']['data'] # [AA]

            # Loop over spectra data points
            for spec in list(out[key].keys()):
                if spec in ['E_eV', 'lambda_A']:
                    continue

                #out[key][spec]['j(ph/ion/s/AA)'] = (
                #    out[key][spec]['j(W/ion/eV)']
                #    / cnt.e
                #    /out[key]['lambda_A']
                #    )
                out[key][spec]['j(ph/s/cm^3/eV)'] = (
                    out[key][spec]['j(W/cm^3/eV)']
                    / out[key]['E_eV']['data']
                    )

    # Output
    return out

# Reads the spectrum per charge state
def _read_indspec(
    out=None,
    file=None,
    Ebeam = None,
    ):

    # Reads file
    ff = open(file, 'r')

    # Init
    key = 'indspec'
    if key not in out.keys():
        out[key] = {}
    out[key]['spec_%0.2fkV'%(Ebeam)] = {}
    tmp = out[key]['spec_%0.2fkV'%(Ebeam)]

    # Charge states
    xx = ff.readline()
    dd = xx.strip().split()
    tmp['neli'] = float(dd[1])
    tmp['nelf'] = float(dd[2])

    # Abundances
    xx = ff.readline()
    dd = xx.strip().split()
    tmp['abund'] = np.asarray([float(dd[xz]) for xz in range(1, len(dd))])

    # Ions
    xx = ff.readline()
    xx = ff.readline()
    dd = xx.strip().split()
    tmp['ions'] = dd[1:]

    # Emissivity
    xx = ff.readline()

    EE = []
    sp = []

    while len(xx.strip()) >0:
        dd = xx.strip().split()

        EE.append(float(dd[0]))
        sp.append(
            [float(dd[xz])*tmp['abund'][xz-1] for xz in range(1, len(dd))]
            )

        xx = ff.readline()

    out[key]['E_eV'] = {'data': np.asarray(EE), 'units': 'eV', 'num': len(EE)}

    if 'settings' in out.keys():
        tmp['j(W/cm^3/eV)'] = np.asarray(sp).T/(out['settings']['psize']['data']*1e-4)
    else:
        tmp['j(W/cm^2/eV)'] = np.asarray(sp).T

    tmp['Ebeam'] = {'data': Ebeam, 'units': 'kV'}

    # Output
    return out

# Reads the simulation settings and total spectrum
def _read_spec(
    out=None,
    file=None,
    Ebeam=None,
    ):

    # Read file
    ff = open(file, 'r')

    # Settings
    key = 'settings'
    out[key] = {}
    tmp = out[key]

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
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'frac'}     # cfrac    

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
    tmp[dd[0]] = {'data': float(dd[1])*1e-7, 'units': 'W/cm^3'}     # ploss

    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1]), 'units': 'eV'}     # <e>

    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {
        'bb': float(dd[1])*1e-7, 
        'bf': float(dd[2])*1e-7, 
        'ff': float(dd[3])*1e-7, 
        'units': 'W/cm^3'
        }     # bb,bf,ff

    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1])*1e-7, 'units': 'W/cm^3'}     # non-E1

    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1])*1e-7, 'units': 'W/cm^3'}     # dn=0

    xx = ff.readline()
    dd = xx.strip().split()
    tmp[dd[0]] = {'data': float(dd[1])*1e-7, 'units': 'W/cm^3'}     # dn=0(g)

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
    tmp[dd[0]] = {'data': float(dd[1])*1e-7, 'units': 'charge'}     # <Z>

    # Charge balance
    key = 'spec'
    if key not in out.keys():
        out[key] = {}
    out[key]['spec_%0.2fkV'%(Ebeam)] = {}
    tmp = out[key]['spec_%0.2fkV'%(Ebeam)]
    tmp['cs_bal'] = {}

    xx = ff.readline()
    xx = ff.readline()
    xx = ff.readline()

    while len(xx.strip()) > 0:
        dd = xx.strip().split()

        nele = _get_nele(sym=dd[0])

        tmp['cs_bal'][nele] = {'Xz': float(dd[1]), 'XzgSC': float(dd[2])}
        xx = ff.readline()
        
    # Emissivity spectrum
    EE = []
    tot = []
    bb = []
    bf = []
    frfr = []
    pp = []

    xx = ff.readline()
    xx = ff.readline()
    xx = ff.readline()

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

    tmp['E_eV'] = {'data': np.asarray(EE), 'units': 'eV', 'num': len(EE)}

    tmp['j(W/cm^3/eV)'] = np.asarray(tot)/(out['settings']['psize']['data']*1e-4)

    tmp['bb'] = {
        'j(W/cm^3/eV)': np.asarray(bb)/(out['settings']['psize']['data']*1e-4)
        }
    tmp['bf'] = {
        'j(W/cm^3/eV)': np.asarray(bf)/(out['settings']['psize']['data']*1e-4)
        }
    tmp['ff'] = {
        'j(W/cm^3/eV)': np.asarray(frfr)/(out['settings']['psize']['data']*1e-4)
        }
    tmp['Planck'] = {
        'j(W/cm^3/eV)': np.asarray(pp)/(out['settings']['psize']['data']*1e-4)
        }
    
    tmp['Ebeam'] = {'data': Ebeam, 'units': 'kV'}

    ff.close()

    # NOTE: I'm just skipping the opacities at the end lol

    return out



        
    




