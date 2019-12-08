from input import *
import numpy as np
import matplotlib.pyplot as plt

def get_potential(S, Tau):
    # Input:
    # S (List): A list of collecive variables at every tau steps
    # tau (number): Frequency of gaussian deposition
    # Output:
    # V (List): The potential trace

    X = np.linspace(0,1,100)
    T = np.arange(Tau[0]-1,Tau[-1]-1,100)
    for t in T:
        V = []
        for x in X:
            v = meta_w*sum([np.exp(-(x-S[tau])**2/(2*meta_sigma**2)) for tau in Tau if tau < t])
            V.append(v)
        plt.plot(X,V)

    plt.show()

    return V

