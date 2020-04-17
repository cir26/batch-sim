import numpy as np


def heat_cap(x, cp):  # J/mol-k
    return np.sum([i*j for i,j in zip(x,cp)])


def cp_v(A, B, C, D, E, T):
    return A+B*(C/T/np.sinh(C/T))**2+D*(E/T/np.cosh(E/T))**2


def cp_l(A, B, C, D, E, T):
    return A+B*T+C*(T**2)+D*(T**3)+E*(T**4)


def heat_vap(A, B, C, D, E, Tc, T):  # j/mol
    Tr = T/Tc
    return A*(1-Tr)**(B+C*Tr+D*(Tr**2)+E*(Tr**3))


