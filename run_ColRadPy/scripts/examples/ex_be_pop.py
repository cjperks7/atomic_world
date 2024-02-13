'''
In this example Be II (Be+) is used because it is a simple system that
has a parent ion (the next charge state) that has a metastable. This 
allows all of the different functionality to be shown and tested.

'''
 
import sys
sys.path.append('../') #starting in 'examples' so need to go up one
from colradpy import colradpy
import numpy as np

fil = 'input/cpb03_ls#be0.dat' #adf04 file
#fil = 'input/be3_adf04'
temperature_arr = np.linspace(1,100,100) #eV
metastable_levels = np.array([0])   #metastable level, just ground chosen here
density_arr =     np.array([1.e13,4.e14]) # cm-3

#calling the colradpy class with the various inputs
be0 = colradpy(
    fil,
    metastable_levels,
    temperature_arr,
    density_arr,
    use_recombination=True,
    use_recombination_three_body = True,
    use_ionization=True,
    suppliment_with_ecip=True
    )

be0.solve_cr() #solve the CR equations with the quasistatic method

be0.plot_pec_sticks([0],[0],[0]) # dim(Te, ne, metastable)

be0.plot_pec_ratio_temp(0,1,[0],[0]) # dim(pec1, pec2, ne, metastable)

be0.plot_pec_ratio_dens(0,1,[0],[0]) # dim(pec1, pec2, Te, metastable)

import matplotlib.pyplot as plt
plt.ion()
met = 0 #metastable 0, this corresponds to the ground state
te = 0 #first temperature in the grid
ne = 0 #frist density in the grid

'''
Shape of be1.data['processed']['pecs']
    dim 1 -> pec index
    dim 2 -> metastables in ion Z, metastables in ion Z+1
    dim 3 -> Te
    dim 4 -> ne

    Multiply PEC by fractional abundance of relevant ion

'''

#fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
#fig.subplots_adjust(bottom=0.15,top=0.92,left=0.105,right=0.965)
fig, ax1 = plt.subplots()
ax1.vlines(be0.data['processed']['wave_air'],
           np.zeros_like(be0.data['processed']['wave_air']),
           be0.data['processed']['pecs'][:,met,te,ne])
ax1.set_xlim(0,1000)
ax1.set_title('PEC spectrum  T$_e$=' +str(be0.data['user']['temp_grid'][te])+\
              ' eV  ne=' + "%0.*e"%(2,be0.data['user']['dens_grid'][ne]) + ' cm$^{-3}$',size=10)
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('PEC (ph cm$^{-3}$ s$^{-1}$)')

print(np.shape(be0.data['processed']['wave_air']),
      np.shape(be0.data['processed']['pecs']),
      np.shape(be0.data['processed']['pec_levels']))
#(320,) (320, 3, 1, 1) (320, 2)

upper_ind = 7 #ninth excited state
lower_ind = 0  #ground state

pec_ind = np.where( (be0.data['processed']['pec_levels'][:,0] == upper_ind) &\
                    (be0.data['processed']['pec_levels'][:,1] == lower_ind))[0]

#plot the temeprature dependence of the chosen pec at first density in the grid
#fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
#fig.subplots_adjust(bottom=0.15,top=0.93,left=0.105,right=0.965)
fig, ax1 = plt.subplots()
ax1.set_title('Temperature dependence of line ' +\
              str(be0.data['processed']['wave_air'][pec_ind]) +' nm',size=10)
ax1.plot(be0.data['user']['temp_grid'],be0.data['processed']['pecs'][pec_ind[0],met,:,ne])
ax1.set_xlabel('Temperature (eV)')
ax1.set_ylabel('PEC (ph cm$^{-3}$ s$^{-1}$)')

#plot the density dependence of the chosen pec at first density in the grid
#fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
#fig.subplots_adjust(bottom=0.15,top=0.93,left=0.125,right=0.965)
fig, ax1 = plt.subplots()
ax1.set_title('Density dependence of line ' +\
              str(be0.data['processed']['wave_air'][pec_ind]) +' nm',size=10)
ax1.plot(be0.data['user']['dens_grid'],be0.data['processed']['pecs'][pec_ind[0],met,te,:])
ax1.set_xlabel('Density (cm$^{-3}$)')
ax1.set_ylabel('PEC (ph cm$^{-3}$ s$^{-1}$)')

#want to find the index of Be I line at 351.55
pec_ind = np.where( (be0.data['processed']['wave_air'] <352) &\
                    (be0.data['processed']['wave_air'] >351))
print('Wavelength from file ' + str(be0.data['processed']['wave_air'][pec_ind[0]]))
#Wavelength from file [351.55028742]
print('PEC upper and lower levels '+ str(be0.data['processed']['pec_levels'][pec_ind[0]]))
#PEC upper and lower levels [[25  2]]