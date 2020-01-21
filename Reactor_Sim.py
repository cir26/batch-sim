# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 18:05:49 2019

@author: macmc
"""

from scipy.integrate import odeint
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from Kinetics import dif_rxn, kinetic_const, rxn, rxn_BT
from Thermo import Thermo_EquilB

class Polymerization:
    #g/mol
    molar_masses = np.array([18.01528,  # W
                    113.159,  # CL
                    226.318,  # CD
                    60.05256,  # AA
                    131.174,  # P1
                    113.1595,  # BACA
                    114.1674,  # TN
                    130.1668,  # TCO
                    43.04522]) # TA

    def __init__(self, Temperature, Initial_Charge, End_Time, ideal=True, P=None, units=None):
        Start_Time = 0
        N=End_Time*10*3600
        if ideal==True: # assumes all species remain in liquid phase
            state, time = self.reaction_l(Temperature, Initial_Charge, Start_Time, End_Time, N)
            state = np.asarray(state)

            if units == 'kg':
                conversion = np.divide(self.molar_masses, 1000)
                state = state*conversion

            if units == 'lb':
                conversion = np.divide(self.molar_masses, 453.592)
                state = state*conversion


            self.t = np.divide(time, 3600)
            self.W = state[:,0]
            self.CL = state[:,1]
            self.CD = state[:,2]
            self.AA = state[:,3]
            self.P1 = state[:,4]
            self.BACA = state[:,5]
            self.TN = state[:,6]
            self.TCO = state[:,7]
            self.TA = state[:,8]



            num = np.add(np.add(np.add(self.BACA, self.TN), self.TCO), self.P1)
            denom = np.add(self.TCO, self.P1)
            self.DP = np.divide(num, denom)

            self.OP2 = np.multiply(self.TCO, np.divide(self.TN, np.add(self.BACA, self.TN)))
            self.OP3 = np.multiply(self.OP2, np.divide(self.BACA, np.add(self.BACA, self.TN)))

            self.Nylon = np.add(np.add(np.add(self.BACA, self.TN), self.TCO), self.TA)

        if ideal==False: # assumes some species enter vapor phase
            if P is None:
                raise ValueError('Please Specify System Pressure')
            else:
                pass

            state, time = self.reaction_BT(Temperature, Initial_Charge, Start_Time, End_Time, N, P, 100)
            state = np.asarray(state)

            if units == 'kg':
                conversion = np.divide(self.molar_masses, 1000)
                state = state*conversion

            if units == 'lb':
                conversion = np.divide(self.molar_masses, 453.592)
                state = state*conversion

            self.t = np.divide(time, 3600)
            self.W = state[:,0]
            self.CL = state[:,1]
            self.CD = state[:,2]
            self.AA = state[:,3]
            self.P1 = state[:,4]
            self.BACA = state[:,5]
            self.TN = state[:,6]
            self.TCO = state[:,7]
            self.TA = state[:,8]

            num = np.add(np.add(np.add(self.BACA, self.TN), self.TCO), self.P1)
            denom = np.add(self.TCO, self.P1)
            self.DP = np.divide(num, denom)

            self.OP2 = np.multiply(self.TCO, np.divide(self.TN, np.add(self.BACA, self.TN)))
            self.OP3 = np.multiply(self.OP2, np.divide(self.BACA, np.add(self.BACA, self.TN)))

            self.Nylon = np.add(np.add(np.add(self.BACA, self.TN), self.TCO), self.TA)

    def reaction_l(self, pars, initial_charge, start_t, end_t, incr):
        t = np.linspace(start_t, end_t*3600, incr) # initialize time vector in seconds

        '''
        molar_masses = [18.01528,  # Water [dW, dCL, dCD, dAA, dP1, dBACA, dTN, dTC, dTA]
                        113.159,  # CL
                        226.318,  # CD
                        60.05256,  # AA
                        131.174,  # P1
                        113.1595,  # BACA
                        114.1674,  # TN
                        130.1668,  # TCO
                        43.04522]  # TA
        '''
        total_mass = sum(initial_charge) # kg
        initial_conc = [charge * 1000 / self.molar_masses[i] / total_mass for i,charge in enumerate(initial_charge)] # mol of species / kg of mixture

        def dNylon_varT(t,state):
            pass

        def dNylon(state,t):
            T = pars
            #W = state[0]
            #CL = state[1]
            #CD = state[2]
            #AA = state[3]
            #P1 = state[4]
            #BACA = state[5]
            #TN = state[6]
            #TC = state[7]
            #TA = state[8]
            k, K = kinetic_const(T, state[7])
            reactions = rxn(state, k, K)
            current_state = dif_rxn(reactions)
            return current_state


        state_solved = odeint(dNylon, initial_conc, t)
        state_solved = np.multiply(state_solved, total_mass)

        return state_solved, t

    def reaction_BT(self, pars, initial_charge, start_t, end_t, incr, P, max_iter):

        t = np.linspace(start_t, end_t*3600, incr)

        total_mass = sum(initial_charge) # kg
        initial_conc = [charge * 1000 / self.molar_masses[i] / total_mass for i,charge in enumerate(initial_charge)] # mol of species / kg of mixture


        def dNylon(state,t):
            T = pars
            #W = state[0]
            #CL = state[1]
            #CD = state[2]
            #AA = state[3]
            #P1 = state[4]
            #BACA = state[5]
            #TN = state[6]
            #TC = state[7]
            #TA = state[8]
            k, K = kinetic_const(T, state[7])
            reactions = rxn_BT(state, k, K, T, P, max_iter)
            current_state = dif_rxn(reactions)
            return current_state


        state_solved = odeint(dNylon, initial_conc, t)
        state_solved = np.multiply(state_solved, total_mass)

        return state_solved, t



# run simulation
# set initial charges in kg
state_dict={'W':3,
            'CL':3*240,
            'CD':0,
            'AA':0,
            'P1':0,
            'BACA':0,
            'TN':0,
            'TC':0,
            'TA':0,}
state = [i for i in state_dict.values()] # extract initial charges for input
Poly = Polymerization(473, state, 30, ideal=False, P=101325, units='kg')
Poly2 = Polymerization(473, state, 30, ideal=True, units='kg')


fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(Poly.t, Poly.W, 'k', alpha = 0.85, linestyle='dashdot', label='NRTL Binary Model, Pressure = 1 atm')
ax1.plot(Poly2.t, Poly2.W, 'k', label='Ideal')
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.legend()
ax1.set_ylabel('Water (kg)')
ax1.set_xlabel('Time (hours)')
plt.show()