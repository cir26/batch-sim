# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 18:05:49 2019
Edited through: Fri May 14:47:38 2020
@author: Robert McMillin (mcmillinre3) and Cristian Romero (cir26)
"""
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from Kinetics import dif_rxn
from Calculations import attr_calc
from time import time



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

    def __init__(self, Temperature, Initial_Cond, End_Time, ideal=True, P=None, units=None):
        dt = 1  # dt = .1 sec
        N = (End_Time*(1/dt)*3600) + 1  # starting at t=zero
        if ideal is True:
            state, time, self.T = self.reaction_l(Temperature, Initial_Cond, End_Time, N)
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

            state, time, self.T, self.P = self.reaction_BT(Temperature, Initial_Cond, End_Time, N, P)
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

    # ideal, DEPRECATED
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

    # NRTL model
    def reaction_BT(self, T, initial_charge, end_t, n, Pres):
        t = np.round(np.linspace(0, end_t*3600, int(n)), 3)
        dt = (end_t*3600)/(n-1)
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
# set initial charge in kg for caprolactam (cap)
# ideal ratios of other reagents already included
cap = 1500
cha = 0.0001*cap
w = 0.01*cap
state_dict = {'W':w,                # Water
              'CL':cap,             # Caprolactam
              'CD':0,               # Cyclic dimer
              'AA':0,               # Acetic acid
              'P1':0,               # Polymer, terminated 1 chain length
              'BACA':0,             # Bounded n6 polymer
              'TN':0,               # Terminating amine group
              'TC':0,               # Terminating carboxylic acid
              'TA':0,               # Terminating acetic acid
              'CHA':cha,            # Cyclohexylamine
              'TCHA':0,             # Terminating cyclohexylamine
              'Temp':273.15+90,     # Temperature in Kelvin
              'Press':5*101325}     # Pressure in Pascal

units = 'kg'  # convert final units
init_cond = [i for i in state_dict.values()]  # extract initial conditions for input
# initialize polymerization object to run simulation with given conditions
Poly = Polymerization(Temperature=273.15+255,
                      Initial_Cond=init_cond,
                      End_Time=10,
                      ideal=False,
                      P=1*101325,
                      units=units)


#######################################################################################################################
# suggested plots
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(Poly.t[115:-1], Poly.Enthalpy[115:]/1000, 'k', alpha=0.85, linestyle='dashdot', label='$\mathregular{\Delta}$ Enthalpy')
#ax1.plot(Poly.t, Poly.H_r/1000, 'k', alpha=0.4, linestyle='dashdot', label='$\mathregular{\Delta}$ Heat of Reaction')
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.ticklabel_format(axis='y',style='sci',scilimits=(5,5))
ax1.legend()
ax1.set_xlabel('Time (hours)')
ax1.set_ylabel('Energy (kJ)')
ax1.set_xlim([0, 10])
plt.show()


fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(Poly.t[115:3444], Poly.Enthalpy[115:3444]/1000, 'k', alpha=0.85, linestyle='dashdot', label='$\mathregular{\Delta}$ Enthalpy')
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
lines = ax1.get_lines()
ax1.set_xlabel('Time (hours)')
ax1.set_ylabel('Energy (kJ)')
ax1.set_xlim([0, 1])



ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(Poly.t, Poly.T-273.15, 'darkred', alpha=0.85, linestyle='solid', label='Temperature (deg C)')
lines += ax2.get_lines()
ax2.add_artist(plt.legend([lines[i] for i in [0, 1]], ['$\mathregular{\Delta}$ Enthalpy', 'Temperature (deg C)'], loc=1))
ax2.tick_params(axis='y', labelcolor='darkred')
ax2.minorticks_on()
ax2.set_ylabel('Temperature (deg C)', color='darkred', rotation=-90, labelpad=15)
ax2.set_ylim([0, 350])
fig.tight_layout()
plt.show()

heat_density = 687  # heat density of Therminol hot oil (kg/m^3)
oil = ((Poly.heat_mass/heat_density)*264.172)  # gal
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(Poly.t[115:3444], oil[115:], 'k', alpha=0.85, linestyle='dashdot', label="Therminol XP heating oil @ 304 deg C")
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.legend()
ax1.set_xlabel('Time (hours)')
ax1.set_ylabel('Volumetric flow rate (gal/sec)')
ax1.set_xlim([0, 1])
plt.show()



fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(Poly.t[3444:-1], Poly.coolingW_mass*0.264172, 'k', alpha=0.85, linestyle='dashdot', label="Cooling water @ 4 deg C")
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.ticklabel_format(axis='y',style='sci',scilimits=(2,2))
ax1.legend()
ax1.set_xlabel('Time (hours)')
ax1.set_ylabel('Volumetric flow rate (gal/sec)')
ax1.set_xlim([1, 10])
plt.show()

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(Poly.t, Poly.LV, 'k', alpha=0.45, linestyle='dashdot', label='Liquid volume of binary mixture')
ax1.plot(Poly.t, Poly.LV+(Poly.Nylon/1140), 'k', alpha=0.85, linestyle='dashdot', label='Liquid+Nylon volume of binary mixture')
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.legend()
ax1.set_xlabel('Time (hours)')
ax1.set_ylabel('Volume ($\mathregular{m^3}$)')
ax1.set_xlim([0, 10])
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




fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(Poly.t, Poly.VV, 'k', alpha=0.75, linestyle='dashdot', label='Vapor volume of binary mixture $\mathregular{m^3}$')
ax1.plot(Poly.t, Poly.P, 'b',label='Pressure (atm e-1)')
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
lines = ax1.get_lines()
ax1.set_xlabel('Time (hours)')
ax1.set_xlim([0, 10])
ax1.set_ylim([0, 100])

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(Poly.t, Poly.T-273.15, 'darkred', alpha=0.85, linestyle='solid', label='Temperature (deg C)')
lines += ax2.get_lines()
ax2.add_artist(plt.legend([lines[i] for i in [0, 1, 2]], ['Vapor volume of binary mixture ($\mathregular{m^3}$)', 'Pressure (atm $\cdot\mathregular{10^{-1}}$)', 'Temperature (deg C)'], loc=2))
ax2.tick_params(axis='y', labelcolor='darkred')
ax2.minorticks_on()
ax2.set_ylabel('Temperature (deg C)', color='darkred', rotation=-90, labelpad=15)
ax2.set_ylim([0, 350])
fig.tight_layout()

plt.show()


fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(Poly.t, Poly.MWW, 'r', alpha = 0.85, linestyle='dashdot', label='Nylon-6 MW-Weight Average (kg/mol)')
ax1.plot(Poly.t, Poly.MWZ, 'k', alpha = 0.65, linestyle='dashdot', label='Nylon-6 MW-Z Average (kg/mol)')
ax1.plot(Poly.t, Poly.MWN, 'k', alpha = 0.3, linestyle='dashdot', label='Nylon-6 MW-Number Average (kg/mol)')
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.legend()
ax1.set_xlabel('Time (hours)')
ax1.set_ylabel('Molecular weight (kg/mol)')
ax1.set_xlim([1,10])
ax1.set_ylim([0,30])
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
ax1.set_xlim([0,10])

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(Poly.t, Poly.W, 'b', alpha=0.5, linestyle='dashdot', label='NRTL Binary Water (kg)')
ax2.plot(Poly.t, Poly.CD, 'r', alpha=0.85, linestyle='dashdot', label='NRTL Binary Cyclic Dimer (kg)')
lines += ax2.get_lines()
ax2.add_artist(plt.legend([lines[i] for i in [0,1,2,3]],['Nylon-6','Caprolactam','Water','Cyclic Dimer'],loc=1))
ax2.tick_params(axis='y', labelcolor='k')
ax2.minorticks_on()
ax2.set_ylim([0,20])
ax2.set_ylabel('Water & Cyclic Dimer mass (kg)', color='k', rotation=-90, labelpad=15)
plt.title('Nonrandom Two-Liquid (NRTL) Binary Model')

plt.show()
