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
from Kinetics import dif_rxn
from Calculations import attr_calc


class Polymerization:

    # g/mol
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
            state_moles = np.asarray(state)  # moles of species
            if units == 'kg':
                self.molar_masses = np.divide(self.molar_masses, 1000)  # kg/mol
                state = state_moles*self.molar_masses  # kg of species
            if units == 'lbs':
                self.molar_masses = np.divide(self.molar_masses, 453.592)  # lbs/mol
                state = state_moles*self.molar_masses  # lbs of species
            # perform addition and calculation of important data attributes to Polymerization object
            self = attr_calc(self, state_arr=state, time_arr=time)

        if ideal == False:  # assumes some species enter vapor phase
            if P is None:
                raise ValueError('Please Specify System Pressure')
            else:
                pass
            state, time, self.T, self.P = self.reaction_BT(Temperature, Initial_Charge, End_Time, N, P, 100)
            self.P=np.divide(self.P,101325)  # convert pressure units to atm
            self.P=np.multiply(self.P,100)   # convert pressure units to atm e-2
            state_moles = np.asarray(state)  # moles of species
            if units == 'kg':
                self.molar_masses = np.divide(self.molar_masses, 1000)
                state = state_moles * self.molar_masses  # kg of species
            if units == 'lbs':
                self.molar_masses = np.divide(self.molar_masses, 453.592)
                state = state_moles * self.molar_masses  # lbs of species
            # perform addition and calculation of important data attributes to Polymerization object
            self = attr_calc(self, state_arr=state, time_arr=time)

    # ideal
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

    # NRTL
    def reaction_BT(self, pars, initial_charge, end_t, incr, Pres, max_iter):
        t = np.linspace(0, end_t*3600, incr)
        total_mass = sum(initial_charge[:-2])
        initial_conc = [initial_charge[i] * 1000 / self.molar_masses[i] / total_mass for i in range(len(initial_charge)-2)]  # mol species / kg total mixture
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


# test simulation
# set initial charges in kg
state_dict = {'W': 13,
              'CL':750,
              'CD':0,
              'AA':1e-5,
              'P1':0,
              'BACA':0,
              'TN':0,
              'TC':0,
              'TA':0,
              'CHA':0,
              'TCHA':0,
              'Temp':273.15+90,
              'Press':5*101325}
units = 'kg' # convert final units
init_cond = [i for i in state_dict.values()] # extract initial conditions for input
Poly = Polymerization([273.15+255, 1.4e-6], init_cond, 24, ideal=False, P=1*101325, units=units)
Poly2 = Polymerization([273.15+255, 1.4e-6], init_cond, 24, ideal=True, units=units)

#plots
fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.plot(Poly.t, Poly.MWW, 'k', alpha = 0.85, linestyle='dashdot', label='NRTL Binary Nylon-6 MWW (kg/mol)')
ax1.plot(Poly2.t, Poly2.MWW, 'k', label='Ideal Nylon-6 MWW (kg/mol)')
ax1.plot(Poly.t, Poly.MWN, 'g', alpha = 0.85, linestyle='dashdot', label='NRTL Binary Nylon-6 MWN (kg/mol)')
ax1.plot(Poly2.t, Poly2.MWN, 'g', label='Ideal Nylon-6 MWN (kg/mol)')
#ax2.plot(Poly.t, Poly.FAV, 'r', alpha = 0.85, linestyle='dashdot', label='NRTL Binary Nylon-6 FAV')
#ax2.plot(Poly2.t, Poly2.FAV, 'r', label='Ideal Nylon-6 FAV')
#ax1.plot(Poly.t, Poly.T,'r',label='Temperature (K)')
#ax1.plot(Poly.t, Poly.P,'b',label='Pressure (atm)')
# ax1.plot(Poly.t, Poly.TCO, 'k', alpha = 0.85, linestyle='dashdot', label='NRTL Binary Nylon-6 ({})'.format(units))
# ax1.plot(Poly2.t, Poly2.TCO, 'k', label='Ideal Nylon-6 ({})'.format(units))
#ax1.plot(Poly.t, Poly.BACA, 'k', alpha = 0.85, linestyle='dashdot', label='NRTL Binary Nylon-6 ({})'.format(units))
#ax1.plot(Poly2.t, Poly2.BACA, 'k', label='Ideal Nylon-6 ({})'.format(units))
# ax1.plot(Poly.t, Poly.TN, 'r', alpha = 0.85, linestyle='dashdot', label='NRTL Binary Nylon-6 ({})'.format(units))
# ax1.plot(Poly2.t, Poly2.TN, 'r', label='Ideal Nylon-6 ({})'.format(units))
# ax1.plot(Poly.t, Poly.P1, 'b', alpha = 0.85, linestyle='dashdot', label='NRTL Binary Nylon-6 ({})'.format(units))
# ax1.plot(Poly2.t, Poly2.P1, 'b', label='Ideal Nylon-6 ({})'.format(units))
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.legend()
ax1.set_xlabel('Time (hours)')
ax1.set_xlim([1,24])
ax1.set_ylim([0,50])
plt.show()



