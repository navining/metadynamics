# ------------System----------------
# mass
M = 48.0

# number of Particles
Ncube = 4
<<<<<<< HEAD
#N = 4* (Ncube ** 3)
N = Ncube ** 3
=======
N = Ncube**3
#N = 4* (Ncube ** 3)
>>>>>>> 6c1fc448c380bc84f463f16530cda8a0e2b899f7

# box side length
L = 1.56 * Ncube #(density of solid argon) (lattice parameter L/Ncube = 1.56) #sigma argon = 3.4 A


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
fileName = 'T' + str(Ta)

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
meta_tau = 50

# cutoff radius
<<<<<<< HEAD
#meta_rc = (3**0.5)*L/(Ncube*1.9)
#meta_rc = L/Ncube/(2**0.5)
meta_rc = 1.2 * (2**(1/6.0))
=======
meta_rc = 1.34
>>>>>>> 6c1fc448c380bc84f463f16530cda8a0e2b899f7
