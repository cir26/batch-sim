# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 18:05:49 2019

@author: macmc
"""
#%% Class Import and Run
from scipy.integrate import odeint
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from Kinetics import dif_rxn, kinetic_const, rxn, rxn_BT


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
                    43.04522,  # TA
                    99.17, # CHA
                    98.17]) # TCHA

    def __init__(self, Temperature, Initial_Charge, End_Time, ideal=True, P=None, units=None):
        N=End_Time*10*3600
        if ideal == True:
            state, time, self.T = self.reaction_l(Temperature, Initial_Charge, End_Time, N)
            state = np.asarray(state)

            if units == 'kg':
                conversion = np.divide(self.molar_masses, 1000)
                state = state*conversion
            if units == 'lbs':
                conversion = np.divide(self.molar_masses, 453.592)
                state = state*conversion

            self.t = np.divide(time, 3600)
            self.W = state[:,0] # units of mass
            self.CL = state[:,1]
            self.CD = state[:,2]
            self.AA = state[:,3]
            self.P1 = state[:,4]
            self.BACA = state[:,5]
            self.TN = state[:,6]
            self.TCO = state[:,7]
            self.TA = state[:,8]
            self.CHA = state[:,9]
            self.TCHA = state[:,10]

            num = np.add(np.add(np.add(self.BACA, self.TN), self.TCO), self.P1)
            denom = np.add(self.TCO, self.P1)
            self.DP = np.divide(num, denom)

            self.OP2 = np.multiply(self.TCO, np.divide(self.TN, np.add(self.BACA, self.TN)))
            self.OP3 = np.multiply(self.OP2, np.divide(self.BACA, np.add(self.BACA, self.TN)))

            self.Nylon = np.add(np.add(np.add(np.add(np.add(self.BACA, self.TN), self.TCO), self.TA), self.TCHA), self.P1)

        if ideal==False: # assumes some species enter vapor phase
            if P is None:
                raise ValueError('Please Specify System Pressure')
            else:
                pass

            state, time, self.T, self.P = self.reaction_BT(Temperature, Initial_Charge, End_Time, N, P, 100)
            self.P=np.multiply(np.divide(self.P,101325),100) # convert pressure units to atm e-2
            state = np.asarray(state)

            if units == 'kg':
                conversion = np.divide(self.molar_masses, 1000)
                state = state*conversion

            if units == 'lbs':
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
            self.CHA = state[:,9]
            self.TCHA = state[:,10]

            num = np.add(np.add(np.add(self.BACA, self.TN), self.TCO), self.P1)
            denom = np.add(self.TCO, self.P1)
            self.DP = np.divide(num, denom)

            self.OP2 = np.multiply(self.TCO, np.divide(self.TN, np.add(self.BACA, self.TN)))
            self.OP3 = np.multiply(self.OP2, np.divide(self.BACA, np.add(self.BACA, self.TN)))

            self.Nylon = np.add(np.add(np.add(np.add(np.add(self.BACA, self.TN), self.TCO), self.TA), self.TCHA), self.P1)

    def reaction_l(self, pars, initial_charge, end_t, incr):
        t = np.linspace(0, end_t*3600, incr)
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
        total_mass = sum(initial_charge[:-2])
        initial_conc = [initial_charge[i] * 1000 / self.molar_masses[i] / total_mass for i in range(len(initial_charge)-2)] # mol of species / kg of total mixture
        initial_conc.append(initial_charge[11])

        def dNylon(state, t):
            max_T = pars[0]
            T_rate = pars[1]
            #W = state[0]
            #CL = state[1]
            #CD = state[2]
            #AA = state[3]
            #P1 = state[4]
            #BACA = state[5]
            #TN = state[6]
            #TC = state[7]
            #TA = state[8]
            #CHA = state[9]
            #TCHA = state[10]
            current_state = dif_rxn(state, max_T, T_rate, t=t, term_delay=0)
            return current_state

        state_solved = np.asarray(odeint(dNylon, initial_conc, t))
        materials = state_solved[:,:11]
        T = state_solved[:,11]
        moles = np.multiply(materials, total_mass)

        return moles, t, T

    def reaction_BT(self, pars, initial_charge, end_t, incr, Pres, max_iter):
        t = np.linspace(0, end_t*3600, incr)
        total_mass = sum(initial_charge[:-2])
        initial_conc = [initial_charge[i] * 1000 / self.molar_masses[i] / total_mass for i in range(len(initial_charge)-2)]
        initial_conc.append(initial_charge[11])
        initial_conc.append(initial_charge[12])

        def dNylon(state, t):
            max_T = pars[0]
            T_rate = pars[1]
            #W = state[0]
            #CL = state[1]
            #CD = state[2]
            #AA = state[3]
            #P1 = state[4]
            #BACA = state[5]
            #TN = state[6]
            #TC = state[7]
            #TA = state[8]
            #CHA = state[9]
            #TCHA = state[10]
            current_state = dif_rxn(state, max_T, T_rate, t=t, P_rate=1.3e-4, min_P = Pres, P_delay=3600*1, term_delay=0)
            return current_state

        state_solved = np.asarray(odeint(dNylon, initial_conc, t))
        materials = state_solved[:,:11]
        T = state_solved[:,11]
        P = state_solved[:,12]
        moles = np.multiply(materials, total_mass)

        return moles, t, T,P

# W = state[0]
# CL = state[1]
# CD = state[2]
# AA = state[3]
# P1 = state[4]
# BACA = state[5]
# TN = state[6]
# TC = state[7]
# TA = state[8]
# CHA = state[9]
# TCHA = state[10]
# Temp = state[11]
# Pressure = state[12]

# run simulation
# set initial charges in kg
# state_dict={'W':10, # input mass in kg
#             'CL':750,
#             'CD':0,
#             'AA':720*4*10e-5,
#             'P1':0,
#             'BACA':0,
#             'TN':0,
#             'TC':0,
#             'TA':0,
#             'CHA':0,
#             'TCHA':0,
#             'Temp':273.15+90,
#             'Press':5*101325}
# units='lbs'
# state = [i for i in state_dict.values()] # extract initial conditions for input
# Poly = Polymerization([273.15+255, 1.4e-6], state, 10, ideal=False, P=1*101325, units=units)
# Poly2 = Polymerization([273.15+255, 1.4e-6], state, 10, ideal=True, units=units)


#plot
# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# ax1.plot(Poly.t, Poly.Nylon, 'k', alpha = 0.85, linestyle='dashdot', label='NRTL Binary Nylon-6 ({})'.format(units))
# ax1.plot(Poly2.t, Poly2.Nylon, 'k', label='Ideal Nylon-6 ({})'.format(units))
# ax1.plot(Poly.t, Poly.T,'r',label='Temperature (K)')
# ax1.plot(Poly.t, Poly.P,'b',label='Pressure (atm e-2)')
# ax1.minorticks_on()
# ax1.tick_params(axis='x', which='minor', direction='in')
# ax1.tick_params(axis='y', which='minor', direction='in')
# ax1.yaxis.set_ticks_position('both')
# ax1.xaxis.set_ticks_position('both')
# ax1.tick_params(direction="in")
# ax1.legend()
# ax1.set_xlabel('Time (hours)')
# plt.show()
#
