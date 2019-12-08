import matplotlib.pyplot as plt
import numpy as np
from processing import *

if __name__ == '__main__':
    t = []
    tau = []
    temperature = []
    E = []
    Q6 = []
    Q6_tau = []

    fileName = 'tau100'

    f = open(fileName,'r')
    try:
        lines=f.readlines()
    finally:
        f.close()

    for i in range(1,len(lines)):
        line = lines[i].strip('\n').split(', ')
        t.append(float(line[0]))
        temperature.append(float(line[1]))
        E.append(float(line[3]))
        Q6.append(float(line[4]))

    f = open(fileName + '_Q6', 'r')
    lines = f.readlines()
    f.close()

    for i in range(1, len(lines)):
        line = lines[i].strip('\n').split(', ')
        tau.append(int(line[0]))
        Q6_tau.append(float(line[1]))

    get_potential(Q6, tau)

    # plt.plot(t,Q6)
    # plt.show()
    # plt.plot(t,E)
    # plt.show()
