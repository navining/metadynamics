import numpy as np

# The simulation will require most of the functions you have already
# implemented above. If it helps you debug, feel free to copy and
# paste the code here.
def my_distance(drij):
    """
    Compute length of displacement vector drij
    assume drij already accounts for PBC

    Args:
      drij (np.array) : vector(s) of length 3
    Returns:
      float: length (distance) of vector(s)
    """
    return np.linalg.norm(drij, axis=0)


def my_disp_in_box(drij, lbox):
    """
    Impose minimum image condition on displacement vector drij=ri-rj

    Args:
      drij (np.array): length-3 displacement vector ri-rj
      lbox (float): length of cubic cell
    Returns:
      np.array: drij under MIC
    """

    return drij - lbox * np.round(drij / lbox)


def my_pos_in_box(pos, lbox):
    """ wrap positions inside simulation box

    Args:
      pos (np.array): positions, shape (natom, ndim)
      lbox (float): box side length
    Returns:
      np.array: pos in box
    """

    return -lbox / 2 + (pos - lbox / 2) % lbox


def my_kinetic_energy(vel, mass):
    """ Calculate total kinetic energy.

    Args:
      vel (np.array): particle velocities, shape (natom, ndim)
      mass (float): particle mass
    Return:
      float: total kinetic energy
    """
    k = 0.0
    for i in vel:
        k += sum(0.5 * mass * i ** 2)

    return k


def my_potential_energy(rij, rc):
    """ Calculate total potential energy.

    Args:
      rij (np.array): distance table, shape (natom, natom)
    Return:
      float: total potential energy
    """
    vshift = 4 * rc ** (-6) * (rc ** (-6) - 1)
    potential = 0.0
    for i in range(len(rij)):
        for j in range(i + 1, len(rij[0])):
            r = rij[i][j]
            if r <= rc:
                potential += 4 * r ** (-6) * (r ** (-6) - 1) - vshift

    return potential


def my_force_on(i, pos, lbox, rc):
    """
    Compute force on atom i

    Args:
      i (int): particle index
      pos (np.array) : particle positions, shape (natom, ndim)
      lbox (float): side length of cubic box
    Returns:
      np.array: force on atom i, a length-3 vector
    """
    Force = np.zeros(3)
    cur = pos[i]
    for atom in pos:
        if (atom == cur).all():
            continue
        r_ij = my_disp_in_box(cur - atom, lbox)
        r = my_distance(r_ij)
        if r <= rc:
            Force += 24 * r ** (-8) * (2 * r ** (-6) - 1) * r_ij
    return Force


def get_distance_table(N, R, lbox):
    drij = np.zeros((N, N))
    for i in range(N):
        for j in range(i + 1, N):
            disp = my_disp_in_box(R[i] - R[j], lbox)
            distance = my_distance(disp)
            drij[i][j] = distance
            drij[j][i] = distance
    return drij

def my_temperature(k,N):
    """ Calculate system temperature.

           Args:
             Ek (float): kinetic energy
             atoms (integer): number of atoms
           Return:
             float: temperature
           """
    return k/(3*N/2)

def my_pressure(V,N,T,R,F):
    viral = np.mean([np.dot(R[i],F[i]) for i in range(N)])
    return (N*T + viral/3) / V