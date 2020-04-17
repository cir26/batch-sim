import numpy as np
from Kinetics import load_data, rxn, kinetic_const
from Thermo import Thermo_EquilB
from Enthalpy import heat_cap, cp_v, cp_l, heat_vap

# function to perform addition and calculation of important data attributes to Polymerization object
def attr_calc(self, state_arr, time_arr, ideal=True):
    # set species mass data as object attributes
    self.t = np.divide(time_arr, 3600)  # convert time to units of hours
    self.W = state_arr[:, 0]
    self.CL = state_arr[:, 1]
    self.CD = state_arr[:, 2]
    self.AA = state_arr[:, 3]
    self.P1 = state_arr[:, 4]
    self.BACA = state_arr[:, 5]
    self.TN = state_arr[:, 6]
    self.TCO = state_arr[:, 7]
    self.TA = state_arr[:, 8]
    self.CHA = state_arr[:, 9]
    self.TCHA = state_arr[:, 10]
    # species moles
    self.W_m = self.state_moles[:, 0]
    self.CL_m = self.state_moles[:, 1]
    self.CD_m = self.state_moles[:, 2]
    self.AA_m = self.state_moles[:, 3]
    self.P1_m = self.state_moles[:, 4]
    self.BACA_m = self.state_moles[:, 5]
    self.TN_m = self.state_moles[:, 6]
    self.TCO_m = self.state_moles[:, 7]
    self.TA_m = self.state_moles[:, 8]
    self.CHA_m = self.state_moles[:, 9]
    self.TCHA_m = self.state_moles[:, 10]

    # calculate mass of oligomers
    self.OP2 = np.multiply(self.TCO, np.divide(self.TN, np.add(self.BACA, self.TN)))
    self.OP3 = np.multiply(self.OP2, np.divide(self.BACA, np.add(self.BACA, self.TN)))
    # calculate polymer (N6) mass
    self.Nylon = np.add(np.add(np.add(np.add(np.add(self.BACA, self.TN), self.TCO), self.TA), self.TCHA), self.P1)
    # calculate polymer (N6) moles
    self.Nylon_m = np.add(np.add(np.add(np.add(np.add(self.BACA_m, self.TN_m), self.TCO_m), self.TA_m), self.TCHA_m), self.P1_m)

    # calculate number-average degree of polymerization
    num = np.add(np.add(np.add(self.BACA, self.TN), self.TCO), self.P1)
    denom = np.add(self.TCO, self.P1)
    self.DP = np.divide(num, denom)
    # calculate weight-average molecular weight of polymer (N6)
    num = np.add(np.add(np.add(np.add(np.add(self.molar_masses[6]**2 * self.TN,
                                             self.molar_masses[7]**2 * self.TCO),
                                             self.molar_masses[5]**2 * self.BACA),
                                             self.molar_masses[4] ** 2 * self.P1),
                                             self.molar_masses[8]**2 * self.TA),
                                             self.molar_masses[10]**2 * self.TCHA)
    denom = np.multiply(np.add(np.add(np.add(self.molar_masses[6]*self.TN,
                                             self.molar_masses[7]*self.TCO),
                                             self.molar_masses[8]*self.TA),
                                             self.molar_masses[10]*self.TCHA),
                                             0.5)
    self.MWW = np.divide(num,denom)
    # calculate weight-average molecular weight of polymer (N6) using moles
    num = np.add(np.add(np.add(np.add(np.add(self.molar_masses[6] ** 2 * self.TN_m,
                                             self.molar_masses[7] ** 2 * self.TCO_m),
                                             self.molar_masses[5] ** 2 * self.BACA_m),
                                             self.molar_masses[4] ** 2 * self.P1_m),
                                             self.molar_masses[8] ** 2 * self.TA_m),
                                             self.molar_masses[10] ** 2 * self.TCHA_m)
    denom = np.multiply(np.add(np.add(np.add(self.molar_masses[6] * self.TN_m,
                                             self.molar_masses[7] * self.TCO_m),
                                             self.molar_masses[8] * self.TA_m),
                                             self.molar_masses[10] * self.TCHA_m),
                                             0.5)
    self.MWW2 = np.divide(num, denom)
    # calculate number-average molecular weight of polymer (N6)
    num = np.add(np.add(np.add(np.add(np.add(self.molar_masses[6]*self.TN,
                                             self.molar_masses[7]*self.TCO),
                                             self.molar_masses[5]*self.BACA),
                                             self.molar_masses[4]*self.P1),
                                             self.molar_masses[8]*self.TA),
                                             self.molar_masses[10]*self.TCHA)
    denom = np.multiply(np.add(np.add(np.add(self.TN,self.TCO),self.TA),self.TCHA),0.5)
    self.MWN = np.divide(num,denom)
    # calculate number-average molecular weight of polymer (N6) using moles
    num = np.add(np.add(np.add(np.add(np.add(self.molar_masses[6] * self.TN_m,
                                             self.molar_masses[7] * self.TCO_m),
                                             self.molar_masses[5] * self.BACA_m),
                                             self.molar_masses[4] * self.P1_m),
                                             self.molar_masses[8] * self.TA_m),
                                             self.molar_masses[10] * self.TCHA_m)
    denom = np.multiply(np.add(np.add(np.add(self.TN_m, self.TCO_m), self.TA_m), self.TCHA_m), 0.5)
    self.MWN2 = np.divide(num, denom)
    # calculate Z-average molecular weight of polymer (N6)
    num = np.add(np.add(np.add(np.add(np.add(self.molar_masses[6] ** 3 * self.TN,
                                             self.molar_masses[7] ** 3 * self.TCO),
                                             self.molar_masses[5] ** 3 * self.BACA),
                                             self.molar_masses[4] ** 3 * self.P1),
                                             self.molar_masses[8] ** 3 * self.TA),
                                             self.molar_masses[10] ** 3 * self.TCHA)
    denom = np.multiply(np.add(np.add(np.add(self.molar_masses[6] ** 2 * self.TN,
                                             self.molar_masses[7] ** 2 * self.TCO),
                                             self.molar_masses[8] ** 2 * self.TA),
                                             self.molar_masses[10] ** 2 * self.TCHA),
                                             0.5)
    self.MWZ = np.divide(num, denom)
    # calculate Z-average molecular weight of polymer (N6) using moles
    num = np.add(np.add(np.add(np.add(np.add(self.molar_masses[6] ** 3 * self.TN_m,
                                             self.molar_masses[7] ** 3 * self.TCO_m),
                                             self.molar_masses[5] ** 3 * self.BACA_m),
                                             self.molar_masses[4] ** 3 * self.P1_m),
                                             self.molar_masses[8] ** 3 * self.TA_m),
                                             self.molar_masses[10] ** 3 * self.TCHA_m)
    denom = np.multiply(np.add(np.add(np.add(self.molar_masses[6] ** 2 * self.TN_m,
                                             self.molar_masses[7] ** 2 * self.TCO_m),
                                             self.molar_masses[8] ** 2 * self.TA_m),
                                             self.molar_masses[10] ** 2 * self.TCHA_m),
                                             0.5)
    self.MWZ2 = np.divide(num, denom)
    # calculate polydispersity index (PDI)
    self.PDI = np.divide(self.MWW,self.MWN)
    # calculating FAV
    self.FAV = np.multiply(9.38e-9,(2*self.MWN)**(2.15))

    if ideal is False:  # perform phase equilibria and enthalpy calculations
        # phase equilibrium constants
        const = {'A_w':5.459,
                 'B_w':3.0542e-1,
                 'C_w':6.4713e2,
                 'D_w':8.1e-2,
                 'A_cap':7.1180e-1,
                 'B_cap':2.54e-1,
                 'C_cap':8.06e2,
                 'D_cap':2.857e-1}
        V, L, xw, xcap, yw, ycap, t12, t21 =[], [], [], [], [], [], [], []
        for w,cl,temp,pres in zip(self.state_moles[:, 0], self.state_moles[:, 1], self.T, self.P):
            zw = w / (w + cl)
            zcap = cl / (w + cl)
            F = w + cl  # FEED
            # Adjust Water and Caprolactam Concentrations with NRTL
            V_t, L_t, F, xw_t, xcap_t, yw_t, ycap_t, _, _, t12_t, t21_t = Thermo_EquilB(temp, pres, F, zw, zcap, 100)
            V.append(V_t[0])
            L.append(L_t[0])
            xw.append(xw_t[0])
            xcap.append(xcap_t[0])
            yw.append(yw_t[0])
            ycap.append(ycap_t[0])
            t12.append(t12_t)
            t21.append(t21_t)
        self.V = np.asarray(V)
        self.L = np.asarray(L)
        self.xw = np.asarray(xw)  # liquid mole fraction
        self.xcap = np.asarray(xcap)
        self.yw = np.asarray(yw)  # vapor mole fraction
        self.ycap = np.asarray(ycap)
        self.t12 = np.asarray(t12)
        self.t21 = np.asarray(t21)
        # phase equilibrium data
        '''
        self.t12 = np.asarray(load_data('t12'))
        self.t21 = np.asarray(load_data('t21'))
        self.xw = np.asarray(load_data('xw'))  # liquid mole fraction
        self.xcap = np.asarray(load_data('xcap'))
        self.yw = np.asarray(load_data('yw'))  # vapor mole fraction
        self.ycap = np.asarray(load_data('ycap'))
        self.L = np.asarray(load_data('L'))
        self.V = np.asarray(load_data('V'))
        '''
        #self.xw_mole = np.asarray((self.xw/self.molar_masses[0])/((self.xw/self.molar_masses[0])+(self.xcap/self.molar_masses[1])))
        #self.xcap_mole = np.asarray((self.xcap/self.molar_masses[1])/((self.xw/self.molar_masses[0])+(self.xcap/self.molar_masses[1])))
        # Liquid Molar Volume of mixture (m^3/mol)
        self.LMV = np.asarray(self.xw * (1/(const['A_w']/(const['B_w']**(1+(1-(self.T/const['C_w']))**const['D_w']))))/1000 + \
                   self.xcap * (1/(const['A_cap']/(const['B_cap']**(1+(1-(self.T/const['C_cap']))**const['D_cap']))))/1000)
        self.LV = np.asarray(self.LMV * self.L)  # m^3
        # Liquid density (kg/m^3)
        self.density_L = ((self.xw * self.W_m * self.molar_masses[0]) +
                          (self.xcap * self.CL_m * self.molar_masses[1])) / (self.LMV * ((self.xw * self.W_m) + (self.xcap * self.CL_m)))
        # Vapor molar volume (m^3/mol)
        self.VMV = 8.314*self.T/self.P
        self.VV = self.VMV * self.V  # m^3

        # enthalpy calculations
        cpv_const = {'AA_A': 4.02e4,
                     'AA_B': 1.3675e5,
                     'AA_C': 1.262e3,
                     'AA_D': 7.003e4,
                     'AA_E': 5.697e2,
                     'CL_A': 7.0664e4,
                     'CL_B': 3.7774e5,
                     'CL_C': -1.5631e3,
                     'CL_D': 2.4215e5,
                     'CL_E': 7.6957e2,
                     'W_A': 3.3359e4,
                     'W_B': 2.6798e4,
                     'W_C': 2.6093e3,
                     'W_D': 8.888e3,
                     'W_E': 1.1676e3,
                     'CHA_A': 1.9744e5,
                     'CHA_B': 2.5599e5,
                     'CHA_C': 1.1911e3,
                     'CHA_D': -1.499e7,
                     'CHA_E': 2.3342e1,
                     'N2_A': 2.9105e4,
                     'N2_B': 8.6149e3,
                     'N2_C': 1.7016e3,
                     'N2_D': 1.0347e2,
                     'N2_E': 9.0979e2}

        # calculate specific heat for vapor phase of AA, CL, W, CHA, ?N2?
        self.cpv_AA = (1/1000)*np.asarray([
            cp_v(cpv_const['AA_A'], cpv_const['AA_B'], cpv_const['AA_C'], cpv_const['AA_D'], cpv_const['AA_E'], temp)
            for temp in self.T])
        self.cpv_CL = (1/1000)*np.asarray([
            cp_v(cpv_const['CL_A'], cpv_const['CL_B'], cpv_const['CL_C'], cpv_const['CL_D'], cpv_const['CL_E'], temp)
            for temp in self.T])
        self.cpv_W = (1/1000)*np.asarray([
            cp_v(cpv_const['W_A'], cpv_const['W_B'], cpv_const['W_C'], cpv_const['W_D'], cpv_const['W_E'], temp) for
            temp in self.T])
        self.cpv_CHA = (1/1000)*np.asarray([
            cp_v(cpv_const['CHA_A'], cpv_const['CHA_B'], cpv_const['CHA_C'], cpv_const['CHA_D'], cpv_const['CHA_E'],
                 temp) for temp in self.T])
        self.cpv_N2 = (1/1000)*np.asarray([
            cp_v(cpv_const['N2_A'], cpv_const['N2_B'], cpv_const['N2_C'], cpv_const['N2_D'], cpv_const['N2_E'], temp)
            for temp in self.T])

        cpl_const = {'AA_A': 5.01e4,
                     'AA_B': 2.456e2,
                     'AA_C': 0,
                     'AA_D': 0,
                     'AA_E': 0,
                     'CL_A': 6.4384e4,
                     'CL_B': 5.1457e2,
                     'CL_C': 0,
                     'CL_D': 0,
                     'CL_E': 0,
                     'W_A': 2.7637e5,
                     'W_B': -2.0901e3,
                     'W_C': 8.125,
                     'W_D': -1.4116e-2,
                     'W_E': 9.3701e-6,
                     'CHA_A': 2.703e4,
                     'CHA_B': 5.4915e2,
                     'CHA_C': 0,
                     'CHA_D': 0,
                     'CHA_E': 0,
                     'N2_A': -3.34e4,
                     'N2_B': 3.507e3,
                     'N2_C': -4.67e1,
                     'N2_D': 2.127e-1,
                     'N2_E': 0}

        # calculate specific heat for liquid phase of AA, CL, W, CHA, ?N2?(J/mol-K)
        self.cpl_AA = (1/1000)*np.asarray([
            cp_l(cpl_const['AA_A'], cpl_const['AA_B'], cpl_const['AA_C'], cpl_const['AA_D'], cpl_const['AA_E'], temp)
            for temp in self.T])
        self.cpl_CL = (1/1000)*np.asarray([
            cp_l(cpl_const['CL_A'], cpl_const['CL_B'], cpl_const['CL_C'], cpl_const['CL_D'], cpl_const['CL_E'], temp)
            for temp in self.T])
        self.cpl_W = (1/1000)*np.asarray([
            cp_l(cpl_const['W_A'], cpl_const['W_B'], cpl_const['W_C'], cpl_const['W_D'], cpl_const['W_E'], temp) for
            temp in self.T])
        self.cpl_CHA = (1/1000)*np.asarray([
            cp_l(cpl_const['CHA_A'], cpl_const['CHA_B'], cpl_const['CHA_C'], cpl_const['CHA_D'], cpl_const['CHA_E'],
                 temp) for temp in self.T])
        self.cpl_N2 = (1/1000)*np.asarray([
            cp_v(cpl_const['N2_A'], cpl_const['N2_B'], cpl_const['N2_C'], cpl_const['N2_D'], cpl_const['N2_E'], temp)
            for temp in self.T])

        # calculate specific heat for Nylon (J/mol-K)
        self.cp_nylon = [mwn * (0.1526 * temp + 223.95) for temp,mwn in zip(self.T, self.MWN)]

        # calculate vapor heat capacity (J/mol-K)
        #self.cpv = [cl * x_cl + w * x_w for cl, x_cl, w, x_w in zip(self.cpv_CL, self.ycap, self.cpv_W, self.yw)]
        # calculate liquid heat capacity (J/mol-K)
        #self.cpl = [cl * x_cl + w * x_w for cl, x_cl, w, x_w in zip(self.cpl_CL, self.xcap, self.cpl_W, self.xw)]

        # calculate heat of vaporization for W (J/mol)
        self.heat_vap_W = (1/1000)*np.asarray([heat_vap(5.2053e7,3.199e-1,-2.12e-1,2.5795e-1,0,647,temp) for temp in self.T])

        # calculate heat of reaction
        rxn_const = {
            1: [1.66e2, 8.32e4, 1.2e4, 7.87e4, 8.03e3, -33.01],
            2: [5.26e6, 9.74e4, 3.37e6, 8.65e4, -2.49e4, 3.951],
            3: [7.93e5, 9.56e4, 4.55e6, 8.42e4, -1.69e4, -29.08],
            4: [2.38e8, 1.76e5, 6.47e8, 1.57e5, -4.02e4, -60.79],
            5: [7.14e4, 8.92e4, 8.36e5, 8.54e4, -1.33e4, 2.439]
        }
        '''
        W = state[0]
        CL = state[1]
        CD = state[2]
        AA = state[3]
        P1 = state[4]
        BACA = state[5]
        TN = state[6]
        TC = state[7]
        TA = state[8]
        CHA = state[9]
        TCHA = state[10]
        '''


        self.H_r = []
        H1 = 8.0287e3
        H2 = -2.4889e4
        H3 = -1.6927e4
        H4 = -4.0186e4
        H5 = -1.3266e4
        for temp, pres, cond in zip(self.T, self.P, self.state_moles):
            cond = np.append(cond, [temp, pres])
            k, K = kinetic_const(temp, cond[7])
            Reactions = rxn(cond, k, K, temp, pres, ideal=False, max_iter=100)
            R1 = Reactions[0]
            R2 = Reactions[1]
            R3 = Reactions[2]
            R4 = Reactions[3]
            R5 = Reactions[4]
            R6 = Reactions[5]
            R7 = Reactions[6]
            R8 = Reactions[7]
            R9 = Reactions[8]
            R10 = Reactions[9]
            R11 = Reactions[10]
            R12 = Reactions[11]
            R13 = Reactions[12]
            R14 = Reactions[13]
            R15 = Reactions[14]
            self.H_r.append(R1*H1 + H2*(R2 + R3 + R4 + R5 + R11 + R12 + R14 + R15) + H3*(R6 + R7 + R13) + H4*R8 + H5*(R9 + R10))
        self.H_r = np.asarray(self.H_r)  # J/(kg-s)
        total_mass = np.sum(state_arr[0])
        # calculate net system enthalpy (J)
        '''
        self.term1=[]
        self.term2 = []
        self.term3 = []
        self.term4 = []
        self.term5 = []
        self.term6 = []
        for i in range(0, len(self.T) - 1):
            self.term1.append(self.H_r[i] * total_mass/3600)
            self.term2.append((self.T[i+1]-self.T[i]) * (self.Nylon_m[i] * self.cp_nylon[i]))
            self.term3.append((self.T[i+1]-self.T[i]) * (self.CL_m[i] * self.cpv_CL[i] * self.ycap[i]))
            self.term4.append((self.T[i+1]-self.T[i]) * (self.W_m[i] * self.cpv_W[i] * self.yw[i]))
            self.term5.append((self.T[i+1]-self.T[i]) * (self.CL_m[i] * self.cpl_CL[i] * self.xcap[i]))
            self.term6.append((self.T[i+1]-self.T[i]) * (self.W_m[i] * self.cpl_W[i] * self.xw[i]))
        self.term1 = np.asarray(self.term1)
        self.term2 = np.asarray(self.term2)
        self.term2[np.isnan(self.term2)] = 0
        self.term3 = np.asarray(self.term3)
        self.term4 = np.asarray(self.term4)
        self.term5 = np.asarray(self.term5)
        self.term6 = np.asarray(self.term6)
        '''
        self.Enthalpy = np.asarray(
            [(self.H_r[i] * self.density_L * self.LV/3600) + (self.T[i+1]-self.T[i]) *
                                               ((self.Nylon_m[i] * self.cp_nylon[i]) +
                                                (self.CL_m[i] * self.cpv_CL[i] * self.ycap[i]) +
                                                (self.W_m[i] * self.cpv_W[i] * self.yw[i]) +
                                                (self.CL_m[i] * self.cpl_CL[i] * self.xcap[i]) +
                                                (self.W_m[i] * self.cpl_W[i] * self.xw[i]))
             for i in range(0, len(self.T)-1)]
        )
        check=0
        for i in range(0, len(self.t)):
            if round(self.T[i],0) == round(100+273.15,0) and check == 0:
                self.Enthalpy[i] += self.W_m[i] * self.heat_vap_W[i]
                #print('time: ',self.t[i])
                check += 1
                print('Heat of Vaporization added')
        cp_w = -(1/1000)*cp_l(cpl_const['W_A'], cpl_const['W_B'], cpl_const['W_C'], cpl_const['W_D'], cpl_const['W_E'], 25+273)  # j/mol
        cp_steam = 2782  # kj/kg
        coolingW_mass = []
        steam_mass = []
        for i in range(0, 3452):
            steam_mass.append(self.Enthalpy[i]/1000/cp_steam)
        self.steam_mass = np.asarray(steam_mass)  # kg
        for i in range(3452, len(self.t)-1):
            coolingW_mass.append((self.Enthalpy[i]/cp_w)*18.01528/1000)  # kg
        self.coolingW_mass = np.asarray(coolingW_mass)

    # return object
    return self
