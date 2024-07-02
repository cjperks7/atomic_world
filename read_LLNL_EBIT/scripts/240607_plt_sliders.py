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
    scram_method = 'cs',
    num_ax = 4,
    )

specs, scram = rh.plot_histos(
    specs=specs,
    scram=scram,
    #x_lim = [4300,5000],
    x_line = [2179, 2459, 2595],
    #Te_init = 10.0,
    cs_init = 14,
    spec_init = 14,
    #plt_scram = True,
    scram_method = 'cs',
    num_ax = 4,
    )

specs, scram = rh.plot_histos(
    specs=specs,
    scram=scram,
    x_units = 'AA',
    x_lim = [2.7, 2.78],
    x_line = [2.72, 2.729, 2.736, 2.747],
    x_highlight = [[2.71, 2.755]],
    #x_lim = [4400, 5000],
    #x_line = [4560],
    #x_highlight = [[4500,4575]],
    scram_method = 'cs',
    cs_init = 15,
    spec_init = 14,
    num_ax = 4,
    )

specs, scram = rh.plot_histos(
    specs=specs,
    scram=scram,
    x_lim=[8000,9250],
    #x_lim = [7600,8200],
    scram_method='cs',
    cs_init = 14,
    spec_init = 14,
    num_ax = 4
    )