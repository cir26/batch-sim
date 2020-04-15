# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 18:05:49 2019

@author: macmc
"""
#%% Class Import and Run
from scipy.integrate import odeint
import pickle
import numpy as np
import matplotlib.pyplot as plt
from Kinetics import dif_rxn
from Calculations import attr_calc
from time import time



class Polymerization:
    count=0
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
        dt = 1  # dt = .1 sec
        N = (End_Time*(1/dt)*3600) + 1  # starting at t=zero
        if ideal is True:
            state, time, self.T = self.reaction_l(Temperature, Initial_Charge, End_Time, N)
            self.state_moles = np.asarray(state)  # moles of species
            if units == 'kg':
                self.molar_masses = np.divide(self.molar_masses, 1000)  # kg/mol
                state = self.state_moles*self.molar_masses  # kg of species
            if units == 'lbs':
                self.molar_masses = np.divide(self.molar_masses, 453.592)  # lbs/mol
                state = self.state_moles*self.molar_masses  # lbs of species
            # perform addition and calculation of important data attributes to Polymerization object
            self = attr_calc(self, state_arr=state, time_arr=time)

        else:  # assumes some species enter vapor phase
            if P is None:
                raise ValueError('Please Specify System Pressure')
            else:
                pass
            '''
            # initiate pickles
            t_check = 0
            with open('t_check', 'wb') as file:
                pickle.dump(t_check, file)
            t12 = []
            with open('t12', 'wb') as file:
                pickle.dump(t12, file)
            t21 = []
            with open('t21', 'wb') as file:
                pickle.dump(t21, file)
            xw = []
            with open('xw', 'wb') as file:
                pickle.dump(xw, file)
            yw = []
            with open('yw', 'wb') as file:
                pickle.dump(yw, file)
            xcap = []
            with open('xcap', 'wb') as file:
                pickle.dump(xcap, file)
            ycap = []
            with open('ycap', 'wb') as file:
                pickle.dump(ycap, file)
            L = []
            with open('L', 'wb') as file:
                pickle.dump(L, file)
            V = []
            with open('V', 'wb') as file:
                pickle.dump(V, file)
            '''
            state, time, self.T, self.P = self.reaction_BT(Temperature, Initial_Charge, End_Time, N, P)
            self.state_moles = np.asarray(state)  # moles of species
            if units == 'kg':
                self.molar_masses = np.divide(self.molar_masses, 1000) # kg/mol
                state = self.state_moles * self.molar_masses  # kg of species
            if units == 'lbs':
                self.molar_masses = np.divide(self.molar_masses, 453.592)
                state = self.state_moles * self.molar_masses  # lbs of species
            # perform addition and calculation of important data attributes to Polymerization object
            self = attr_calc(self, state_arr=state, time_arr=time, ideal=ideal)
            self.P=np.divide(self.P,101325)  # convert pressure units to atm
            self.P=np.multiply(self.P,10)   # convert pressure units to atm e-1

    # ideal
    def reaction_l(self, T, initial_charge, end_t, n):
        t = np.round(np.linspace(0, end_t*3600, int(n)), 3)
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
            max_T = T
            T_rate =  1.4e-6
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
        materials = state_solved[:, :11]
        T = state_solved[:, 11]
        moles = np.multiply(materials, total_mass)

        return moles, t, T

    # NRTL
    def reaction_BT(self, T, initial_charge, end_t, n, Pres):
        t = np.round(np.linspace(0, end_t*3600, int(n)), 3)
        dt = (end_t*3600)/(n-1)
        print('dt: ', dt)
        total_mass = sum(initial_charge[:-2])
        initial_conc = [initial_charge[i] * 1000 / self.molar_masses[i] / total_mass for i in range(len(initial_charge)-2)]  # mol species / kg total mixture
        initial_conc.append(initial_charge[11])
        initial_conc.append(initial_charge[12])

        def dNylon(state, time):
            max_T = T
            T_rate = 1.4e-6
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
            current_state = dif_rxn(state, max_T, T_rate, t=time, P_rate=1.3e-4, min_P=Pres, P_delay=3600*1, term_delay=0)
            return current_state
        start_time = time()
        state_solved = np.asarray(odeint(dNylon, initial_conc, t))
        #state_solved = np.asarray(odeint(dNylon, initial_conc, t, hmin=dt, hmax=dt))
        end_time = time()
        print('solver time elapsed {:.4f} s'.format(end_time - start_time))
        materials = state_solved[:, :11]
        T = state_solved[:, 11]
        P = state_solved[:, 12]
        moles = np.multiply(materials, total_mass)
        return moles, t, T, P


#######################################################################################################################
# test simulation
# set initial charges in kg
cap=400
aa= (2/3)*0.0001*cap
w=0.01*cap
state_dict = {'W':w,
              'CL':cap,
              'CD':0,
              'AA':aa,
              'P1':0,
              'BACA':0,
              'TN':0,
              'TC':0,
              'TA':0,
              'CHA':0,
              'TCHA':0,
              'Temp':273.15+90,
              'Press':5*101325}

units = 'kg'  # convert final units
init_cond = [i for i in state_dict.values()]  # extract initial conditions for input
Poly = Polymerization(273.15+255, init_cond, 12, ideal=False, P=1*101325, units=units)
#Poly2 = Polymerization(273.15+255, init_cond, 12, ideal=True, units=units)

#plots
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(Poly.t, Poly.LV, 'k', alpha=0.45, linestyle='dashdot', label='Liquid volume of binary mixture (m^3)')
ax1.plot(Poly.t, Poly.LV+(Poly.Nylon/1140), 'k', alpha=0.85, linestyle='dashdot', label='Liquid+Nylon volume of binary mixture (m^3)')
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.legend()
ax1.set_xlabel('Time (hours)')
ax1.set_ylabel('Volume (m^3)')
ax1.set_xlim([0, 12])
plt.show()




fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(Poly.T-273.15, Poly.t12, 'k', alpha=0.44, linestyle='dashdot', label='tau (W-CL)')
ax1.plot(Poly.T-273.15, Poly.t21, 'k', alpha=0.88, linestyle='dashdot', label='tau (CL-W)')
ax1.plot([0,300],[0,0], 'k', lw=0.5, linestyle='solid')
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.legend()
ax1.set_ylabel('POLYNRTL Tau values')
ax1.set_xlabel('Temperature (deg C)')
ax1.set_xlim([0, 300])
ax1.set_ylim([-0.1,0.08])
plt.show()




#ax2.plot(Poly.t, Poly.FAV, 'r', alpha = 0.85, linestyle='dashdot', label='NRTL Binary Nylon-6 FAV')
#ax2.plot(Poly2.t, Poly2.FAV, 'r', label='Ideal Nylon-6 FAV')
#ax1.plot(Poly2.t, Poly2.TCO, 'k', label='Ideal Nylon-6 ({})'.format(units))
#ax1.plot(Poly.t, Poly.BACA, 'k', alpha = 0.85, linestyle='dashdot', label='NRTL Binary Nylon-6 ({})'.format(units))
#ax1.plot(Poly2.t, Poly2.BACA, 'k', label='Ideal Nylon-6 ({})'.format(units))
# ax1.plot(Poly.t, Poly.TN, 'r', alpha = 0.85, linestyle='dashdot', label='NRTL Binary Nylon-6 ({})'.format(units))
# ax1.plot(Poly2.t, Poly2.TN, 'r', label='Ideal Nylon-6 ({})'.format(units))
# ax1.plot(Poly.t, Poly.P1, 'b', alpha = 0.85, linestyle='dashdot', label='NRTL Binary Nylon-6 ({})'.format(units))
# ax1.plot(Poly2.t, Poly2.P1, 'b', label='Ideal Nylon-6 ({})'.format(units))
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(Poly.t, Poly.VV, 'k', alpha=0.85, linestyle='dashdot', label='Vapor volume of binary mixture (m^3)')
ax1.plot(Poly.t, Poly.P, 'b',label='Pressure (atm e-1)')
#ax1.plot(Poly.t, Poly.T,'r',label='Temperature (K)')
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.legend()
ax1.set_xlabel('Time (hours)')
ax1.set_xlim([0, 12])
plt.show()


fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(Poly.t, Poly.MWW, 'k', alpha = 0.85, linestyle='dashdot', label='NRTL Binary Nylon-6 MW-W (kg/mol)')
ax1.plot(Poly.t, Poly.MWW2, 'r', alpha = 0.85, linestyle='dashdot', label='NRTL Binary Nylon-6 MW-W-2 (kg/mol)')
ax1.plot(Poly.t, Poly.MWZ, 'k', alpha = 0.5, linestyle='dashdot', label='NRTL Binary Nylon-6 MW-Z (kg/mol)')
ax1.plot(Poly.t, Poly.MWZ2, 'r', alpha = 0.5, linestyle='dashdot', label='NRTL Binary Nylon-6 MW-Z-2 (kg/mol)')
ax1.plot(Poly.t, Poly.MWN, 'k', alpha = 0.25, linestyle='dashdot', label='NRTL Binary Nylon-6 MW-N (kg/mol)')
ax1.plot(Poly.t, Poly.MWN2, 'r', alpha = 0.25, linestyle='dashdot', label='NRTL Binary Nylon-6 MW-N-2 (kg/mol)')
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.legend()
ax1.set_xlabel('Time (hours)')
ax1.set_ylabel('Molecular weight (kg/mol)')
ax1.set_xlim([1,12])
ax1.set_ylim([0,50])
plt.show()


fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(Poly.t, Poly.Nylon, 'k', alpha=0.85, linestyle='dashdot', label='NRTL Binary Nylon-6 (kg)')
ax1.plot(Poly.t, Poly.CL, 'g', alpha=0.85, linestyle='dashdot', label='NRTL Binary Caprolactam(kg)')
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
lines = ax1.get_lines()
ax1.set_xlabel('Time (hours)')
ax1.set_ylabel('Nylon & Caprolactam mass (kg)', color='k')
ax1.set_xlim([0,12])

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(Poly.t, Poly.W, 'b', alpha=0.5, linestyle='dashdot', label='NRTL Binary Water (kg)')
ax2.plot(Poly.t, Poly.CD+Poly.P1, 'r', alpha=0.85, linestyle='dashdot', label='NRTL Binary Cyclic Dimer (kg)')
lines += ax2.get_lines()
ax2.add_artist(plt.legend([lines[i] for i in [0,1,2,3]],['Nylon-6 (kg)','Caprolactam (kg)','Water (kg)','Cyclic Dimer (kg)'],loc=1))
ax2.tick_params(axis='y', labelcolor='k')
ax2.set_ylim([0,20])
ax2.set_ylabel('Water & Cyclic Dimer mass (kg)', color='k')
plt.title('NRTL Binary Model')
fig.tight_layout()

plt.show()
