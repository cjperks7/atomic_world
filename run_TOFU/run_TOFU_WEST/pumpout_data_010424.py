'''

pumpout_dat_040124.py is a test script to test my understanding of using
the tofu_west module to load XICS data

cjperks
April 1, 2024

'''

# Modules
from atomic_world.run_TOFU import run_TOFU_WEST as rtfw

# Shot #
shot =59984

dout = rtfw.get_west(
    shot = shot,
    plt = False,
    save = True,
    save_path = '/home/CP028959/work/pumpout/shots',
    save_name = str(shot)+'_summary'
    )

'''

# Plots loaded 
dout = np.load(
    '/home/cjperks/2201_Pumpout/WEST/C9_results/shots/'
    +str(shot)+'_summary.npz',
    allow_pickle=True
    )['dout'][()]

dout = rtfw.get_west(
    dout = dout,
    quants = [],
    plt = True,
    t_plt = [5.5,6.5],
    save = False,
    )


'''