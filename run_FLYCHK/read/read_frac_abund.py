'''

Script to read fractional abundance data files produced by FLYCHK

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

# Read fractional abundance data file
def read(
    filepath = None,
    ):

    # Init
    Te = [] # [eV]
    data = []

    # Flags
    take_data = False
    take_info = False

    # Opens data file and loop over lines
    with open(filepath, 'r') as f:
        # Loop over file
        for line in f:

            # If at beginning of new data table
            if line.split(' ')[0] == '#begin' and line.split(' ')[1][0] == 'c':

                data_tmp = []
                cs_tmp = []

            # If at end of new data table
            elif line.split(' ')[0] == '#end':
                take_data = False
                if line.split(' ')[1][0] == 'c':
                    data.append(data_tmp)

            # If temperature line
            if line.split(':')[0] == 'fun':
                Te.append(float(line.split(' ')[-2]))

            # Loads data array
            if take_data:
                cs_tmp.append(int(line.split(' ')[0]))
                data_tmp.append(float(line.split(' ')[1]))

            # If start of data line
            if line == '//\n':
                take_data = True

            # If start of info table
            if line.split(' ')[0] == '#begin' and line.split(' ')[1][0] == 'a':
                take_info = True

            if take_info:
                # Species name
                if line.split(':')[0] == 'tit':
                    sp = line.split(':')[1]

                elif line.split(':')[0] == 'tit2':
                    ne = float(line.split(' ')[4])

    # Closes data file
    f.close()

    # Organizes data
    out = {}

    out['Te_eV'] = np.asarray(Te) # dim(ntemp,)
    out['cs'] = np.asarray(cs_tmp) # dim(ncs,)
    out['species'] = sp
    out['ne_cm3'] = ne

    out['abund'] = np.zeros((len(out['Te_eV']), len(out['cs']))) # dim(ntemp, ncs)

    for nt in np.arange(len(out['Te_eV'])):
        out['abund'][nt,:] = np.asarray(data[nt])

    # Output
    return out




