# ------------System----------------
# mass
M = 48.0

# number of Particles
Ncube = 4
N = Ncube ** 3

# box side length
L = 4.0

# cutoff radius
rc = 2.0

# initial temperature
T0 = 2

# system temperature
Ta = 0.1

# Anderson
anderson = True
eta = 0.3125


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
meta_tau = 100