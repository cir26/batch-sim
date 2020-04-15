import numpy as np


def heat_cap(x,cp):  # J/mol-k
    return np.sum([i*j for i,j in zip(x,cp)])

def cp_v(A,B,C,D,E,T):
    return A+B*(C/T/np.sinh(C/T))**2+D*(E/T/np.cosh(E/T))**2

def cp_l(A,B,C,D,E,T):
    return A+B*T+C*(T**2)+D*(T**3)+E*(T**4)

def heat_vap(A, B, C, D, E, Tc, T):  # j/mol
    Tr = T/Tc
    return A*(1-Tr)**(B+C*Tr+D*(Tr**2)+E*(Tr**3))


cpv_const = {'AA_A':4.02e4,
            'AA_B':1.3675e5,
            'AA_C':1.262e3,
            'AA_D':7.003e4,
            'AA_E':5.697e2,
            'CL_A':7.0664e4,
            'CL_B':3.7774e5,
            'CL_C':-1.5631e3,
            'CL_D':2.4215e5,
            'CL_E':7.6957e2,
            'W_A':3.3359e4,
            'W_B':2.6798e4,
            'W_C':2.6093e3,
            'W_D':8.888e3,
            'W_E':1.1676e3,
            'CHA_A':1.9744e5,
            'CHA_B':2.5599e5,
            'CHA_C':1.1911e3,
            'CHA_D':-1.499e7,
            'CHA_E':2.3342e1,
            'N2_A':2.9105e4,
            'N2_B':8.6149e3,
            'N2_C':1.7016e3,
            'N2_D':1.0347e2,
            'N2_E':9.0979e2}

cpl_const = {'AA_A':5.01e4,
            'AA_B':2.456e2,
            'AA_C':0,
            'AA_D':0,
            'AA_E':0,
            'CL_A':6.4384e4,
            'CL_B':5.1457e2,
            'CL_C':0,
            'CL_D':0,
            'CL_E':0,
            'W_A':2.7637e5,
            'W_B':-2.0901e3,
            'W_C':8.125,
            'W_D':-1.4116e-2,
            'W_E':9.3701e-6,
            'CHA_A':2.703e4,
            'CHA_B':5.4915e2,
            'CHA_C':0,
            'CHA_D':0,
            'CHA_E':0,
            'N2_A':-3.34e4,
            'N2_B':3.507e3,
            'N2_C':-4.67e1,
            'N2_D':2.127e-1,
            'N2_E':0}