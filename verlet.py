import numpy as np

# We have written the Verlet time-stepping functions for you below,
# `h` is the time step.
# -----------------------------------------------------------------

def VerletNextR(r_t, v_t, a_t, h):
    """Return new positions after one Verlet step"""
    # Note that these are vector quantities.
    # Numpy loops over the coordinates for us.
    r_t_plus_h = r_t + v_t * h + 0.5 * a_t * h * h
    return r_t_plus_h


def VerletNextV(v_t, a_t, a_t_plus_h, h):
    """Return new velocities after one Verlet step"""
    # Note that these are vector quantities.
    # Numpy loops over the coordinates for us.
    v_t_plus_h = v_t + 0.5 * (a_t + a_t_plus_h) * h
    return v_t_plus_h
