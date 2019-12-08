import matplotlib.pyplot as plt
import random as rd

import initial as init, verlet as verl
from properties import *
from input import *
from metadynamics import *
from processing import *
from output import *

# Main Loop.
# ------------------------------------------------------------------------

def simulate():
    # -----------------------Initialize--------------------------------
    ## system
    #R = init.InitPositionFCC(Ncube, L)
    #V = init.InitVelocity(N, T0, M)
    R = init.InitPositionFromFile('Liquid_R',N)
    V = init.InitVelocityFromFile('Liquid_V',N)
    E = np.zeros(steps)

    ## metadynamics
    n_gauss = 0
    S = [] # position of the center of the Gaussian

    output(fileName,"steps, temperature, pressure, energy, Q6\n")

    for t in range(0, steps):

        # -----------------------Anderson Thermostat----------------------
        if anderson == True:
            #Tc = T0 + (Ta - T0)/steps * t
            Tc = Ta
            sigma = (Tc / M) ** 0.5
            mean = 0
            for i in V:
                if np.random.random() < eta * h:
                    i[0] = rd.gauss(mean, sigma)
                    i[1] = rd.gauss(mean, sigma)
                    i[2] = rd.gauss(mean, sigma)

        # -----------------------Propagation----------------------------
        ## calculate forces
        F = np.array([my_force_on(i, R, L, rc) for i in range(N)])
        A = F / M

        ## calculate new positions
        nR = verl.VerletNextR(R, V, A, h)
        nR = my_pos_in_box(nR, L)

        ## calculate forces with new positions nR
        nF = np.array([my_force_on(i, nR, L, rc) for i in range(N)])

        ## calculate displacement table
        rij = get_displacement_table(N, nR, L)

        ## calculate distance table
        drij = get_distance_table(N, rij)

        ## bias forces with metadynamics
        if t%50 == 0:
            meta_Q6 = calculate_Q6(nR,drij,rij)
        # metaF, meta_Q6 = meta(t, n_gauss, S, nF, nR, drij, rij)
        # nF = nF + metaF


        ## calculate new velocities
        nA = nF / M
        nV = verl.VerletNextV(V, A, nA, h)

        # update positions:
        R, V = nR, nV

        # -----------------------Measuring Physical Properties----------------------------

        ## calculate kinetic energy contribution
        k = my_kinetic_energy(V, M)

        ## calculate temperature
        T = my_temperature(k, N)

        ## calculate potential energy contribution
        p = my_potential_energy(drij, rc)

        ## calculate total energy
        E[t] = k + p

        ## calculate pressure
        P = my_pressure(L ** 3, N, T, R, nF)


        # ------------------------Output-------------------------------------
        if t % 50 == 0:
            output(fileName, '%d, %.3f, %.3f, %.5f, %.5f\n' % (t, T, P, E[t], meta_Q6))
        if t % 100 == 0:
            write_xyz(fileName + '_' + str(t) + '.xyz', R)
    # write_R(fileName, R)
    # write_V(fileName, V)


    return S


if __name__ == '__main__':
    # Run the simulation
    S = simulate()

    # Post processing
    #postprocessing()
