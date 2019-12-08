# ------------System----------------
# mass
M = 48.0

# number of Particles
Ncube = 4
#N = 4* (Ncube ** 3)
N = Ncube ** 3

# box side length
#L = 1.56 * Ncube #(density of solid argon) (lattice parameter L/Ncube = 1.56) #sigma argon = 3.4 A
L = 4

# cutoff radius
rc = L/2

# initial temperature
T0 = 2

# system temperature
Ta = 0.2

# Anderson
anderson = True
eta = 0.3125

# filename
#fileName = 'T' + str(Ta)
fileName = 'T02'

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
#meta_rc = (3**0.5)*L/(Ncube*1.9)
#meta_rc = L/Ncube/(2**0.5)
meta_rc = 1.2 * (2**(1/6.0))
