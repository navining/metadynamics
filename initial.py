import numpy as np
import random
# Everyone will start their gas in the same initial configuration.
# ----------------------------------------------------------------
def InitPositionCubic(Ncube, L):
    """Places Ncube^3 atoms in a cubic box; returns position vector"""
    N = Ncube ** 3
    position = np.zeros((N, 3))
    rs = L / Ncube
    roffset = L / 2 - rs / 2
    n = 0
    # Note: you can rewrite this using the `itertools.product()` function
    for x in range(0, Ncube):
        for y in range(0, Ncube):
            for z in range(0, Ncube):
                if n < N:
                    position[n, 0] = rs * x - roffset
                    position[n, 1] = rs * y - roffset
                    position[n, 2] = rs * z - roffset
                n += 1
    #print(position)
    return position

def InitPositionFCC(Ncube, L):
    """Places 4 * Ncube^3 atoms in a FCC box; returns position vector"""
    N = 4 * (Ncube ** 3)
    position = np.zeros((N, 3))
    rs = L / Ncube #lattice parameter
    roffset = L / 2 - rs / 2
    n = 0
    # Note: you can rewrite this using the `itertools.product()` function
    for x in range(0, Ncube):
        for y in range(0, Ncube):
            for z in range(0, Ncube):
                if n < N:
                    position[n, 0] = rs * x - roffset
                    position[n, 1] = rs * y - roffset
                    position[n, 2] = rs * z - roffset
                n += 1
                if n < N:
                    position[n, 0] = rs * (x + 0.5) - roffset
                    position[n, 1] = rs * (y + 0.5) - roffset
                    position[n, 2] = rs * z - roffset
                n += 1
                if n < N:
                    position[n, 0] = rs * (x + 0.5) - roffset
                    position[n, 1] = rs * y - roffset
                    position[n, 2] = rs * (z + 0.5) - roffset
                n += 1
                if n < N:
                    position[n, 0] = rs * x - roffset
                    position[n, 1] = rs * (y + 0.5) - roffset
                    position[n, 2] = rs * (z + 0.5) - roffset
                n += 1
    #print("N=", N, len(position))
    #print(position)
    return position

def InitPositionRandom(Ncube, L):
    """Places Ncube^3 atoms in a cubic box; returns position vector"""
    N = Ncube ** 3
    position = np.zeros((N, 3))
    for i in range(N):
        position[i,0] = random.random()*L - 2 / L
        position[i, 1] = random.random() * L - 2 / L
        position[i, 2] = random.random() * L - 2 / L
    return position

def InitVelocity(N, T0, mass=1., seed=1):
    dim = 3
    np.random.seed(seed)
    # generate N x dim array of random numbers, and shift to be [-0.5, 0.5)
    velocity = np.random.random((N, dim)) - 0.5
    sumV = np.sum(velocity, axis=0) / N  # get the average along the first axis
    velocity -= sumV  # subtract off sumV, so there is no net momentum
    KE = np.sum(velocity * velocity)  # calculate the total of V^2
    vscale = np.sqrt(dim * N * T0 / (mass * KE))  # calculate a scaling factor
    velocity *= vscale  # rescale
    return velocity
