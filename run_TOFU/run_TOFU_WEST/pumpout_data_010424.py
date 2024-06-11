'''

pumpout_dat_040124.py is a test script to test my understanding of using
the tofu_west module to load XICS data

cjperks
April 1, 2024

'''

# Modules
from atomic_world.run_TOFU import run_TOFU_WEST as rtfw

# Shot #
shot =59972

dout = rtfw.get_west(
    shot = shot,
    plt = False,
    save = True,
    save_path = '/home/CP028959/work/pumpout/shots',
    save_name = str(shot)+'_summary'
    )

'''
scp CP028959@altair.partenaires.cea.fr:/home/CP028959/work/pumpout/shots/59984_summary.npz .

from atomic_world.run_TOFU import run_TOFU_WEST as rtfw
import numpy as np

shot0 = 59972
shot1 = 60002

for shot in np.arange(shot0, shot1+1):
    dout = rtfw.get_west(
        shot = shot,
        plt = False,
        save = True,
        quants = ['xics'],
        save_path = '/home/CP028959/work/pumpout/XICS',
        save_name = str(shot)+'_raw'
        )


from atomic_world.run_TOFU import run_TOFU_WEST as rtfw
#shot = 59995
shot = 60001

# Plots loaded 
dout = np.load(
    '/home/cjperks/work/2201_Pumpout/WEST/C9_results/shots/'
    +str(shot)+'_summary.npz',
    allow_pickle=True
    )['dout'][()]

dout = rtfw.get_west(
    dout = dout,
    quants = [],
    plt = True,
    t_plt = [4.5,5.5],
    save = False,
    )


'''