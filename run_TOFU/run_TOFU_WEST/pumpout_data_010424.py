'''

pumpout_dat_040124.py is a test script to test my understanding of using
the tofu_west module to load XICS data

cjperks
April 1, 2024

'''

# Modules
import tofu_west as tfw
import imas
import imas_west
import pywed as pw

# Shot number
shot = 59984

# Get equilibrium and wall
equi = imas_west.get(shot, 'equilibrium')
wall = imas_west.get(shot, 'wall')

# Start time
t_ignitron  =  pw.tsmat(shot, 'IGNITRON|1')

# Loads XICS data
ddata_imas = imas_west.get(shot, 'spectrometer_x_ray_crystal')
ddata_tfw = tfw.SpectroX2D.load_data(shot)

# Interferometry
interf = imas_west.get(shot, 'interferometer')

nch_int = len(interf.channel)
t_int = interf.time

ne_m2 = interf.channel[1].n_e_line.data
ne_m3 = interf.channel[1].n_e_line_average.data

# Loads line-integrated interferometry data
ne_m2_tot = np.zeros((len(t_int), nch_int)) # dim(ntime, nch)
for ch in np.arange(nch_int):
    ne_m2_tot[:,ch] = interf.channel[ch].n_e_line.data # [m^-2]

fig, ax = plt.subplots()
ax.plot(t_int, ne_m2_tot)
fig.show()

# ECE
ece = imas_west.get(shot, 'ece')
plt.plot(ece.time,ece.t_e_central.data)

nch_ece = len(ece.channel)
c11 = ece.channel(11)
c11_Te = c11.t_e.data

# Core profiles
corep = imas_west.get(shot, 'core_profiles')

t_core = corep.time

rho_t = corep.profiles_1d[1].grid.rho_tor_norm
ne = corep.profiles_1d[1].electrons.density
te = corep.profiles_1d[1].electrons.temperature

