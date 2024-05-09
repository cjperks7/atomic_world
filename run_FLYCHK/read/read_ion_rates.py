'''

Scripts to read effective ionization/recombination rates posted
for FLYCHK NLTE

cjperks
May 7, 2024

NOTES:
    1) From some reason, there's a line break after csd that breaks
        the logic here, so I just manually delete it

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

# Read ionization rates data file
def read(
    filepath = None,
    ):

    # Init
    Te = [] # [eV]
    data = {}

    # Data types
    types = ['ci', 'ea', 'rr', 'dr', 'bb', 'bf', 'csd']

    # Flags
    take_data = False
    take_info = False

    # Opens data file and loop over lines
    with open(filepath, 'r') as f:
        # Loop over file
        for line in f:

            # If at the end of a data block
            if take_data and line == '\n':
                take_data = False

            # Reads data
            if take_data:
                tmp = [ii for ii in line.split(' ') if ii != '']

                data[key]['Te_eV'].append(float(tmp[0]))
                data[key]['data'].append([float(ii) for ii in tmp[1:]])

            # If at the beginning of a data block
            if line.split(' ')[0] in types and line.split(' ')[1] == 'temp':
                key = line.split(' ')[0]
                take_data = True

                data[key] = {}
                data[key]['Te_eV'] = []
                data[key]['data'] = []

    # Closes data file
    f.close()

    # Formats output
    out = {}

    for key in data.keys():
        out[key] = {}
        out[key]['Te_eV'] = np.asarray(data[key]['Te_eV'])

        out[key]['data'] = np.zeros((len(out[key]['Te_eV']), len(data[key]['data'][-1])))

        for nt in np.arange(len(out[key]['Te_eV'])):
            out[key]['data'][nt,:] = np.asarray(data[key]['data'][nt])

    # Output
    return out