'''

To hack the rocking curve input into XICSRT, this code
rewrites the rocking curve output from ToFu into
the XOP file format

Nov. 29, 2023

'''

# Modules
import tofu.spectro._rockingcurve as rc
from scipy.interpolate import interp1d
import numpy as np
import os

__all__ = [
    'write_xop'
    ]

###########################################
#
#               Main
#
##########################################

def write_xop(
    # File path
    save_path = None,
    save_name = None,
    # Physical parameters
    lamb    = 1.60, # [AA]
    mat     = 'Germanium',
    miller  = np.r_[2, 0, 2],
    d_hkl   = 2.0001752742728575, # [AA]
    ):

    # Material dictonary
    sym = mat[:2]+str(miller[0])+str(miller[1])+str(miller[2])
    dcry = {
        'material': mat,
        'name': mat,
        'symbol': sym,
        'miller': miller,
        'target': {
            'ion': 'Fe25+',
            'lamb': lamb, # e-10
            'units': 'AA',
            },
        'd_hkl': d_hkl,
        }

    # Calculates rocking curve
    dout = rc.compute_rockingcurve(
        crystal = dcry['name'], 
        din=dcry, 
        lamb = dcry['target']['lamb'], 
        plot_power_ratio = False
        )

    # Quantities of interest
    glance = dout['Glancing angles'][0,0,0,:] # dim(angle,)
    sigma = dout['Power ratio'][0,0,0,:] # dim(angle,), perpendicular pol.
    pi = interp1d(
        dout['Glancing angles'][1,0,0,:],
        dout['Power ratio'][1,0,0,:],
        bounds_error = False,
        fill_value = 0.0
        )(glance) # dim(angle,), parallel pol.
    glance -= dout['Bragg angle of reference (rad)']
    glance *= 1e6 # [urad]

    # Opens an ACCII file to write in
    f = open(
        os.path.join(save_path,save_name),
        'w'
        )

    # Header
    f.write("#F \n")
    f.write("#S Material: " + sym + " Wavelength: " + str(lamb) + " AA \n")
    f.write("#C \n")
    f.write("#N 7 \n")

    # Column headers
    f.write(
        "#L  Th-ThB{in} [microrad]"
        + "  Th-ThB{out} [microrad]"
        + "  phase_p[rad]"
        + "  phase_s[rad]"
        + "  Circ Polariz"
        + "  p-polarized"
        + "  s-polarized"
        + "\n"
        )

    # Loop over glancing angles
    for ii in np.arange(len(glance)):
        # Writes data
        f.write(
            "{:1.6E}".format(glance[ii]).rjust(17, ' ')
            + "{:1.6E}".format(glance[ii]).rjust(17, ' ')
            + "{:1.6E}".format(0.0).rjust(17, ' ')
            + "{:1.6E}".format(0.0).rjust(17, ' ')
            + "{:1.6E}".format(0.0).rjust(17, ' ')
            + "{:1.8E}".format(pi[ii]).rjust(17, ' ')
            + "{:1.8E}".format(sigma[ii]).rjust(17, ' ')
            + "\n"
            )

    # Closes file
    f.close()



    



