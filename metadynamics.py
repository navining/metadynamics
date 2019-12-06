import numpy as np
from properties import *
from input import *
import math
import sys
from scipy.special import *

def get_pcle_neighbors(pcle, R, rc, lbox):
    neighbors = []
    for other_pcle in R:
        if np.array_equal(other_pcle, pcle):
            continue
        elif np.linalg.norm(my_disp_in_box(pcle - other_pcle, lbox)) >= rc:
            continue
        else:
            neighbors.append(other_pcle)

    return neighbors


def calculate_Qlm(l, m, R, rc, lbox):
    Qlm = 0
    tot_neighbor_tracker = 0

    for pcle in R:
        neighbors = get_pcle_neighbors(pcle, R, rc, lbox)
        tot_neighbor_tracker += len(neighbors)

        for neighbor in neighbors:
            r = my_disp_in_box(neighbor - pcle, lbox)
            rhat = r / np.linalg.norm(r)
            theta = np.arccos(rhat[2])
            phi = np.arctan2(rhat[1], rhat[0])

            Qlm += sph_harm(m, l, phi, theta)

    Qlm = Qlm / tot_neighbor_tracker

    return Qlm


def calculate_Ql(l, R, rc, lbox):
    # Input:
    # l (non-negative int): uniquely determines Ql
    # R (numpy.array, size N*3): atom positions
    # rc (float): cutoff distance
    # lbox (float): length of box
    # Output:
    # Q6: the collective variable

    Ql = 0
    for m in range(-l, l, 1):
        Ql += abs(calculate_Qlm(l, m, R, rc, lbox)) ** 2

    Ql = math.sqrt(4 * np.pi * Ql / (2 * l + 1))

    return Ql

def calculate_Q6(R):
    # Input:
    # R (numpy.array, size N*3): atom positions
    # Output:
    # Q6: the collective variable

    Q6 = calculate_Ql(6, R, rc, L)
    return Q6

def calculate_ds_dr(R):
    # Input:
    # R (numpy.array, size N*3): atom positions
    # Output:
    # ds_dr (numpy.array, size N*3): the derivative of s with respect to atom positions

    ds_dr = 0
    return ds_dr

def meta(step, n_gauss, S, force, R):
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
    s = calculate_Q6(R)
    print(s)
    ## calculate the derivative of s with respect to atom positions
    ds_dr = calculate_ds_dr(R)

    ## every tau step, save the value of s
    if step % tau == 0:
        n_gauss += 1
        if n_gauss < max_gauss:
            S.append(s)
        else:
            sys.exit("max_gauss exceeded")

    ## calculate the derivative of the history-dependent potential with respect to s
    dV_ds = 0
    for s_tau in S:
        gauss = w * np.exp(-(s - s_tau)**2/2/sig**2)
        dV_ds = dV_ds + gauss * (s - s_tau)**2/2/sig**2

    ## bias the force
    force = force - dV_ds * ds_dr

    return force

