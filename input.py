# ------------System----------------
# mass
M = 48.0

# number of Particles
Ncube = 4
N = Ncube ** 3

# box side length
L = 14.5369647

# cutoff radius
rc = L/2

# initial temperature
T0 = 1

# system temperature
Ta = 0

# Anderson
anderson = True
eta = 0.3125

# filename
fileName = 'T_0'

# --------------MD-----------------

# total steps
steps = 5000

# gap
h = 0.032

# ------------Metadynamics----------

# Gaussian
meta_w = 0.1 # height
meta_sigma = 0.1 # width

# max number of Gaussian
meta_max = 100

# frequency
meta_tau = 50

# cutoff radius
meta_rc = 1.3469544579712476