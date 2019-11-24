import numpy as np
from input import *
import sys

def calculate_Q6(R):
    # Input:
    # R (numpy.array, size N*3): atom positions
    # Output:
    # Q6: the collective variable

    Q6 = 0
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

