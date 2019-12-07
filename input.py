# ------------System----------------
# mass
M = 48.0

# number of Particles
Ncube = 3
N = 4* (Ncube ** 3)

# box side length
#L = 14.5369647
L = 1.56 * Ncube #(density of solid argon) (lattice parameter L/Ncube = 1.56) #sigma argon = 3.4 A

# cutoff radius
#rc = L/2
rc = (3**0.5)*L/(Ncube*1.9)
# initial temperature
T0 = 0.2

# system temperature
Ta = 0.1

# Anderson
anderson = True
eta = 0.3125

# filename
fileName = 'test'

# --------------MD-----------------

# total steps
steps = 2000

# gap
h = 0.032

# ------------Metadynamics----------

# Gaussian
meta_w = 0.1 # height
meta_sigma = 0.1 # width

# max number of Gaussian
meta_max = 100

# frequency
meta_tau = 5
