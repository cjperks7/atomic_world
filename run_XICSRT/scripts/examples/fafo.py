'''

XICSRT test script

Nov. 28, 2023

'''

import numpy as np
import xicsrt
import time
#xicsrt.warn_version('0.8')

# Track time spent
start = time.time()
start0 = start

# 1.
config = dict()

# 2.
config['general'] = {}
config['general']['number_of_iter'] = 5
config['general']['save_images'] = False

# 3.
config['sources'] = {}
config['sources']['source'] = {}
config['sources']['source']['class_name'] = 'XicsrtSourceDirected'
config['sources']['source']['intensity'] = 1e4
config['sources']['source']['wavelength'] = 3.9492
config['sources']['source']['wavelength_dist'] = 'monochrome'
config['sources']['source']['angular_dist'] = 'isotropic_xy'
config['sources']['source']['spread'] = [np.radians(10.0), np.radians(10.0)]
config['sources']['source']['xsize'] = 0.01
config['sources']['source']['ysize'] = 0.01
config['sources']['source']['zsize'] = 0.01


# 4.
config['optics'] = {}
config['optics']['crystal'] = {}
config['optics']['crystal']['class_name'] = 'XicsrtOpticSphericalCrystal'
config['optics']['crystal']['check_size'] = True
config['optics']['crystal']['origin'] = [0.0, 0.0,         0.80374151]
config['optics']['crystal']['zaxis']  = [0.0, 0.59497864, -0.80374151]
config['optics']['crystal']['xsize']  = 0.2
config['optics']['crystal']['ysize']  = 0.2
config['optics']['crystal']['radius'] = 1.0

# Rocking curve FWHM in radians.
# This is taken from x0h for quartz 1,1,-2,0
# Darwin Curve, sigma: 48.070 urad
# Darwin Curve, pi:    14.043 urad
config['optics']['crystal']['crystal_spacing'] = 2.45676
#config['optics']['crystal']['rocking_type'] = 'gaussian'
#config['optics']['crystal']['rocking_fwhm'] = 48.070e-6

config['optics']['crystal']['rocking_type'] = 'file'
config['optics']['crystal']['rocking_filetype'] = 'xop'
config['optics']['crystal']['rocking_file'] = '/home/cjperks/test.dat'

# 5.
config['optics']['detector'] = {}
config['optics']['detector']['class_name'] = 'XicsrtOpticDetector'
config['optics']['detector']['origin'] = [0.0,  0.76871290, 0.56904832]
config['optics']['detector']['zaxis']  = [0.0, -0.95641806, 0.29200084]
config['optics']['detector']['xsize']  = 0.4
config['optics']['detector']['ysize']  = 0.2

# 6.
results = xicsrt.raytrace(config)

# Print the time this took
stop = time.time()
s='Done %10.3E s'%(stop-start)
print(s)

import xicsrt.visual.xicsrt_3d__plotly as xicsrt_3d
fig = xicsrt_3d.figure()
xicsrt_3d.add_rays(results)
xicsrt_3d.add_optics(results['config'])
xicsrt_3d.add_sources(results['config'])
xicsrt_3d.show()

import xicsrt.visual.xicsrt_2d__matplotlib as xicsrt_2d
xicsrt_2d.plot_intersect(results, 'detector')





# Finds origin to localize on
vect = np.array([-0.98966998,  0.14336434,  0.])
x0 = np.array([2.80209255, 0.03308598, 0.018])

xap = np.array([16.17933996, -1.90475217,  0.018])

R0 = 1.85

s = np.linspace(0.5, 1.5, 1000)

tt = s*vect[0] +x0[0]
yy = s*vect[1] +x0[1]
gg = np.sqrt(tt**2+yy**2)

indg = np.argmin(abs(gg-R0))
s_sol = s[indg]
origin = s_sol*vect+x0

# vos phase space parametization
dk = 0.6
nk = 20
kk = np.linspace(-dk, dk, nk)

nperp = 4
dperp = 0.05
perp = np.linspace(-dp,dp, nperp)


dz = 0.06
nz = 4
zz = np.linspace(-dz,dz, nz)

vz = np.array([0,0,1])
vperp = np.cross(vect, vz)

import itertools as itt

for inds in itt.product(range(nk), range(nperp), range(nz)):
    print(inds[0])
    print(inds[1])
    print(inds[2])

    pt = origin + (
        kk[inds[0]]*vect
        + perp[inds[1]]*vperp
        + zz[inds[2]]*vz
        )

    vpt = xap - pt
    vpt /= np.linalg.norm(vpt)


