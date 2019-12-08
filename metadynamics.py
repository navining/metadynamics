import numpy as np
from properties import *
from input import *
from output import *
import math
import sys
import multiprocessing
import time

from scipy.special import *

def calculate_Qlm(l, m, R, rc, drij, rij):
    Qlm = 0
    bonds = 0
    for i in range(N):
        for j in range(i+1,N):
            r = drij[(i,j)]
            if r <= rc:
                bonds += 1
                rhat = rij[(i,j)] / np.linalg.norm(rij[(i,j)])
                theta = np.arccos(rhat[2])
                phi = np.arctan2(rhat[1], rhat[0])
                Qlm += sph_harm(m, l, phi, theta)
    if bonds!=0:
        Qlm = Qlm / bonds

    return Qlm


def calculate_Ql(l, R, rc ,drij, rij):
    # Input:
    # l (non-negative int): uniquely determines Ql
    # R (numpy.array, size N*3): atom positions
    # rc (float): cutoff distance
    # lbox (float): length of box
    # Output:
    # Q6: the collective variable

    Ql = 0
    for m in range(-l, l, 1):
        Ql += abs(calculate_Qlm(l, m, R, rc, drij, rij)) ** 2

    Ql = math.sqrt(4 * np.pi * Ql / (2 * l + 1))

    return Ql

def calculate_Q6(R,drij, rij):
    # Input:
    # R (numpy.array, size N*3): atom positions
    # Output:
    # Q6: the collective variable

    Q6 = calculate_Ql(6, R, meta_rc, drij, rij)
    return Q6

def update(parameters):
    position, R, s, d = parameters
    nR = R.copy()
    # update atom positions, displacement table and distance table
    nR[position] += d
    nrij = get_displacement_table(N,nR,L)
    ndrij = get_distance_table(N,nrij)
    s_new = calculate_Q6(R, ndrij, nrij)
    dsdr =  (s_new - s) / d
    return position,dsdr

def calculate_ds_dr(R,s,drij,rij):
    # Input:
    # R (numpy.array, size N*3): atom positions
    # Output:
    # ds_dr (numpy.array, size N*3): the derivative of s with respect to atom positions
    d = 0.001
    ds_dr = np.zeros((N,3))

    num_cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(num_cores)

    parameters = []
    for i in range(N):
        for j in range(3):
            parameters.append(((i,j),R,s,d))
    returns = pool.map(update, parameters)

    for x in returns:
        ds_dr[x[0]] = x[1]

    # ds_dr = 0
    return ds_dr

def meta(step, n_gauss, S, force, R, drij, rij):
    # -----------------Parameters----------------
    ## Gaussian
    w = meta_w
    sig = meta_sigma

    ## max number of Gaussian
    max_gauss = meta_max

    ## frequency of Gaussian deposition
    tau = meta_tau

    # ---------------Calculation----------------
    ## calculate CV
    s = calculate_Q6(R,drij, rij)
    ## calculate the derivative of s with respect to atom positions
    ds_dr = calculate_ds_dr(R,s,drij,rij)

    ## every tau step, save the value of s
    if step % tau == 0:
        n_gauss += 1
        if n_gauss < max_gauss:
            S.append(s)
            output(fileName+'_Q6', '%d, %.5f\n' % (step, s))
        else:
            sys.exit("max_gauss exceeded")
    ## calculate the derivative of the history-dependent potential with respect to s
    dV_ds = 0
    for s_tau in S:
        gauss = w * np.exp(-(s - s_tau)**2/2/sig**2)
        dV_ds = dV_ds + gauss * (s - s_tau)**2/2/sig**2

    ## bias the force
    force = - dV_ds * ds_dr

    return force,s

