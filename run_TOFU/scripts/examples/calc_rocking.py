'''

Calculate rocking curves using ToFu

'''

# Modules
import tofu.spectro._rockingcurve as rc

dcry = {
    'material': 'Germanium',
    'name': 'XRSLR',
    'symbol': 'Ge202',
    'miller': np.r_[2., 0., 2.,],
    'target': {
        'ion': 'Fe25+',
        'lamb': 1.615, # e-10
        'units': 'AA',
        },
    'd_hkl': 2.0003,
    }

dout = rc.compute_rockingcurve(
    crystal = dcry['name'], 
    din=dcry, 
    lamb = dcry['target']['lamb'], 
    plot_power_ratio = True
    )