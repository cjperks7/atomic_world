'''

Script to generate slider plots of LLNL EBIT 
microcalorimeter spectra

cjperks
June 7th, 2024

'''

# Modules
from atomic_world.read_LLNL_EBIT.main import read_histos as rh

# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2

# Make histogram plots
specs, scram = rh.plot_histos(
    x_line = [2179],
    #dE = 5, # [eV]
    )

specs, scram = rh.plot_histos(
    specs=specs,
    scram=scram,
    x_lim = [4300,5000],
    x_line = [2179],
    Te_init = 10.0,
    spec_init = 13,
    plt_scram = True
    )