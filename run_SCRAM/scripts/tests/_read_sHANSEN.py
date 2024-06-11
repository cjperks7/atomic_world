'''

Manages reading SCRAM files passed to me by
Stephanie Hansen

cjperks
May 14, 2024

'''

# Modules
import numpy as np
import os
import copy
import scipy.constants as cnt

__all__ = [
    'main'
    ]

#####################################################
#
#               Main
#
#####################################################

def main(
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

    # Reads table file
    if 'table' in quants:
        out = _read_table(
            out=out,
            file = os.path.join(
                filepath,
                name+'table.txt'
                )
            )

    # Reads total spectrum file
    if 'spec' in quants:
        out = _read_spec(
            out=out,
            file = os.path.join(
                filepath,
                name+'spec.txt'
                )
            )

    # Reads ion-resolved spectrum file
    if 'indspec' in quants:
        out = _read_spec(
            out=out,
            file = os.path.join(
                filepath,
                name+'indspec.txt'
                ),
            ind =True
            )

    # Converts units to what I like
    if convert:
        out = _convert(out=out)

    # Output
    return out

#####################################################
#
#               Utilities
#
#####################################################

# Converts emissivty units from W/ion/eV -> ph/ion/s/AA
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

                out[key][spec]['j(ph/ion/s/AA)'] = (
                    out[key][spec]['j(W/ion/eV)']
                    / cnt.e
                    /out[key]['lambda_A']
                    )

    # Output
    return out

# Reads spectrum data file
def _read_spec(
    out=None,
    file=None,
    ind = False,
    ):

    # Init
    if ind:
        key = 'indspec'
    else:
        key = 'spec'
    out[key] = {}

    # Opens file
    f = open(file, 'r')

    # Energy grid
    out[key]['E_eV'] = {}

    xx = f.readline()
    out[key]['E_eV']['num'] = int(xx.split(' ')[-4])

    out[key]['E_eV']['data'] = np.zeros(out[key]['E_eV']['num'])

    for ii in range(out[key]['E_eV']['num']):
        xx = f.readline()
        out[key]['E_eV']['data'][ii] = float(xx)
    xx = f.readline()

    # Data
    first = True
    for ii in range(out['table']['length']):
        # Header titles only at first instance
        if first:
            out[key]['spec_%i'%(ii)] = {}

            # Data headers
            xx = f.readline()
            top = []

            # ne/Te data
            for yy in xx.split(' '):
                if len(yy)>0 and yy!='\n':
                    top.append(yy.strip('\n'))
                    out[key]['spec_%i'%(ii)][yy.strip('\n')] = 0

            # Ion type data
            if ind:
                mid = []
                xx = f.readline()
                for yy in xx.split(' '):
                    if len(yy)>0 and yy!='\n':
                        mid.append(yy.strip('\n'))
                        out[key]['spec_%i'%(ii)][yy.strip('\n')] = 0
                xx = f.readline()
                for yy in xx.split(' '):
                    if len(yy)>0 and yy!='\n':
                        mid.append(yy.strip('\n'))
                        out[key]['spec_%i'%(ii)][yy.strip('\n')] = 0
                xx = f.readline()
                for yy in xx.split(' '):
                    if len(yy)>0 and yy!='\n':
                        mid.append(yy.strip('\n'))
                        out[key]['spec_%i'%(ii)][yy.strip('\n')] = 0

            xx = f.readline()
            bot = []

            # Spectra data
            for yy in xx.split(' '):
                if len(yy)>0 and yy!='\n':
                    bot.append(yy.strip('\n'))
                    out[key]['spec_%i'%(ii)][yy.strip('\n')] = np.zeros(out[key]['E_eV']['num'])

            xx = f.readline()
            xx = f.readline()

            first = False

            out[key]['spec_init'] = copy.deepcopy(out[key]['spec_%i'%(ii)])

        else:
            out[key]['spec_%i'%(ii)] = copy.deepcopy(out[key]['spec_init'])


        # Header data
        cnt = 0
        for yy in xx.split(' '):
            if len(yy)>0 and yy!='\n':
                out[key]['spec_%i'%(ii)][top[cnt]] = float(yy)
                cnt += 1

        if ind:
            cnt = 0
            # Ions simulated
            xx = f.readline()
            for yy in xx.split(' '):
                if len(yy)>0 and yy!='\n':
                    out[key]['spec_%i'%(ii)][mid[cnt]] = int(yy)
                    cnt += 1

            # Reshapes data block, dim(nion, nE)
            for cnt in range(2):
                out[key]['spec_%i'%(ii)][bot[cnt]] = np.zeros((
                    int(out[key]['spec_%i'%(ii)]['nelf'] - out[key]['spec_%i'%(ii)]['neli']+1),
                    out[key]['E_eV']['num']
                    ))

            # Fractional abundance block
            out[key]['spec_%i'%(ii)]['abund'] = np.zeros(
                int(out[key]['spec_%i'%(ii)]['nelf'] - out[key]['spec_%i'%(ii)]['neli']+1)
                )
                
            # Ion labels
            xx = f.readline()
            out[key]['spec_%i'%(ii)]['ions'] = []
            for yy in xx.split(' '):
                if len(yy)>0 and yy!='\n':
                    out[key]['spec_%i'%(ii)]['ions'].append(yy)

        # Fractional abundance block
        if ind:
            xx = f.readline()
            ion = 0

            for yy in xx.split(' '):
                if len(yy)>0 and yy!='\n':
                    out[key]['spec_%i'%(ii)]['abund'][ion] = float(yy)
                    ion +=1

        # Emissivity data block
        for jj in range(out[key]['E_eV']['num']):
            xx = f.readline()

            if ind:
                cnt = 0
                ion = 0

                for yy in xx.split(' '):
                    if len(yy)>0 and yy!='\n':
                        out[key]['spec_%i'%(ii)][bot[cnt]][ion,jj] = float(yy)

                        ion += 1

                        if ion == len(out[key]['spec_%i'%(ii)]['ions']):
                            ion = 0
                            cnt +=1

            else:
                cnt = 0
                for yy in xx.split(' '):
                    if len(yy)>0 and yy!='\n':
                        out[key]['spec_%i'%(ii)][bot[cnt]][jj] = float(yy)
                        cnt += 1

        # Skip line break
        xx = f.readline()
        xx = f.readline()


    # Closes files
    f.close()

    # Deletes placeholder
    del out[key]['spec_init']
    
    # Output
    return out

# Reads summary data table
def _read_table(
    out=None,
    file=None
    ):

    # Init
    out['table'] = {}
    length = 1

    # Opens file
    f = open(file, 'r')

    # Stores SCRAM data
    out['table']['SCRAM'] = {}

    xx = f.readline()
    out['table']['SCRAM']['version'] = xx.split(' ')[-2]

    xx = f.readline()
    out['table']['SCRAM']['structure'] = xx.split(' ')[-1][:-1]

    # Element data
    out['table']['species'] = {}

    xx = f.readline()
    out['table']['species']['name'] = xx.split(' ')[-3]

    xx = f.readline()
    out['table']['species']['Znuclear'] = xx.split(' ')[-4]

    xx = f.readline()
    out['table']['species']['amu'] = xx.split(' ')[0]

    # Electron temperature grid
    out['table']['Te_eV'] = {}

    xx = f.readline()
    out['table']['Te_eV']['num'] = int(xx.split(' ')[-6])
    length *= out['table']['Te_eV']['num'] 

    xx = f.readline()
    out['table']['Te_eV']['data'] = np.asarray([float(yy) for yy in xx.split(' ') if len(yy) > 0])

    # Electron density grid
    out['table']['ne_cm3'] = {}

    xx = f.readline()
    out['table']['ne_cm3']['num'] = int(xx.split(' ')[-6])
    length *= out['table']['ne_cm3']['num']

    xx = f.readline()
    out['table']['ne_cm3']['data'] = np.asarray([float(yy) for yy in xx.split(' ') if len(yy) > 0])

    # Radiation field grid
    out['table']['rad field'] = {}

    xx = f.readline()
    out['table']['rad field']['num_pts'] = int(xx.split(' ')[-6])
    length *= out['table']['rad field']['num_pts']

    xx = f.readline()
    out['table']['rad field']['num_groups'] = int(xx.split(' ')[-7])
    length *= out['table']['rad field']['num_groups']

    xx = f.readline()
    out['table']['rad field']['pts_J/cc'] = float(xx.split(' ')[3])
    out['table']['rad field']['brdy_eV'] = float(xx.split(' ')[-1])

    # Table headers
    xx = f.readline()
    zz = []

    for ii in xx.split(' '):
        if len(ii) >0 and ii !='\n':
            zz.append(ii)
            out['table'][ii] = np.zeros(length)

    # Populates table data
    for ii in range(length+1):
        xx = f.readline()

        cnt = 0
        for yy in xx.split(' '):
            if len(yy) >0 and yy!='\n':
                out['table'][zz[cnt]][ii] = yy
                cnt += 1

    out['table']['length'] = length

    # Closes file
    f.close()

    # Output
    return out