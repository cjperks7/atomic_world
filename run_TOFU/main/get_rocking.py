'''

Function to get rocking curve from ToFu

cjperks
July 17, 2024

'''

# Modules
import tofu.spectro._rockingcurve as rc
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np

__all__ = [
    '_get_rocking'
    ]

###################################################
#
#           Main
#
###################################################

# Calculates rocking curve
def _get_rocking(
    material = None,
    miller = None,
    lamb0 = None, # [AA]
    mix = 0.5,
    plt_all = False,
    ):

    # Crystal dictionary
    dcry = {
        'material': material,
        'name': 'tmp',
        'symbol': 'tmp',
        'miller': miller,
        'target': {
            'ion': 'tmp',
            'lamb': lamb0, # e-10
            'units': 'AA',
            },
        #'d_hkl': 1.1549,
        }

    # Computes rocking curve
    dout = rc.compute_rockingcurve(
        crystal = dcry['name'], 
        din=dcry, 
        lamb = dcry['target']['lamb'], 
        plot_power_ratio = False
        )

    # Get data arrays
    bragg = dout['Bragg angle of reference (rad)'] # [rad]
    angle = dout['Glancing angles'][0,0,0,:] # [rad]
    sigma = dout['Power ratio'][0,0,0,:]
    pi = interp1d(
        dout['Glancing angles'][1,0,0,:],
        dout['Power ratio'][1,0,0,:],
        bounds_error = False,
        fill_value = 0.0
        )(angle)
    pwr = mix*sigma +(1-mix)*pi

    # Plotting
    if plt_all:
        fig, ax = plt.subplots()

        ax.plot(
            (angle-bragg)*180/np.pi,
            sigma,
            label = 'sigma'
            )
        ax.plot(
            (angle-bragg)*180/np.pi,
            pi,
            label = 'pi'
            )
        ax.plot(
            (angle-bragg)*180/np.pi,
            pwr,
            label = 'mix'
            )
        ax.grid('on')
        ax.legend(labelcolor='linecolor')

        ax.set_ylabel('Reflected power [frac]')
        ax.set_xlabel('Reflected angle - Bragg [deg]')

        ax.set_title(material + ' [%i,%i,%i], %1.5f AA'%(miller[0], miller[1], miller[2], lamb0))

    # Output
    return {
        'angle': angle, # [rad]
        'sigma': sigma,
        'pi': pi,
        'pwr': pwr,
        'bragg': bragg, # [rad]
        }
