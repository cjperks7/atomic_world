'''

Script to read cooling curve data files produced by FLYCHK

For NLTE, see:
https://www-amdis.iaea.org/FLYCHK/

'''

# Modules
import numpy as np
import sys, os

__all__ = [
    'read'
    ]

########################################################
#
#                   Main
#
########################################################

# Read cooling curve data file
def read(
    filepath = None,
    ):

    # Init
    ne = []
    data = []

    # Flags
    take_data = False
    take_info = False

    # Opens data file and loop over lines
    with open(filepath, 'r') as f:
        # Loop over file
        for line in f:

            # If at beginning of new data table
            if line.split(' ')[0] == '#begin' and line.split(' ')[1][0] == 'd':

                data_tmp = []
                Te_tmp = []

            # If at end of new data table
            elif line.split(' ')[0] == '#end':
                take_data = False
                if line.split(' ')[1][0] == 'd':
                    data.append(data_tmp)

            # If density line
            if line.split(':')[0] == 'fun':
                ne.append(float(line.split(' ')[-1]))

            # Loads data array
            if take_data:
                Te_tmp.append(float(
                    line.split('  ')[-2]
                    )) # [eV]
                data_tmp.append(float(
                    line.split('  ')[-1]
                    )) # [erg/s/atom]

            # If start of data line
            if line == '//\n':
                take_data = True

            # If start of info table
            if line.split(' ')[0] == '#begin' and line.split(' ')[1][0] == 'a':
                take_info = True

            if take_info:
                # Species name
                if line.split(':')[0] == 'tit':
                    sp = line.split(':')[1].split(' ')[0]

    # Closes data file
    f.close()

    # Organizes data
    out = {}

    out['Te_eV'] = np.asarray(Te_tmp) # dim(ntemp,)
    out['ne_cm3'] = np.asarray(ne) # dim(ndens,)
    out['species'] = sp

    out['cooling'] = np.zeros((len(out['ne_cm3']), len(out['Te_eV']))) # dim(ndens, ntemp)

    erg2J = 1e-7
    for nd in np.arange(len(out['ne_cm3'])):
        out['cooling'][nd,:] = np.asarray(data[nd])*erg2J

    # Output, [W/atom]
    return out