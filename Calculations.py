import numpy as np
from Kinetics import load_data
from Thermo import Thermo_EquilB

# function to perform addition and calculation of important data attributes to Polymerization object
def attr_calc(self, state_arr, time_arr, ideal=True):
    # set species mass data as object attributes
    self.t = np.divide(time_arr, 3600)
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

    if ideal is False:
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
        self.LMV = self.xw * (1/(const['A_w']/(const['B_w']**(1+(1-(self.T/const['C_w']))**const['D_w']))))/1000 + \
                   self.xcap * (1/(const['A_cap']/(const['B_cap']**(1+(1-(self.T/const['C_cap']))**const['D_cap']))))/1000
        self.LV = self.LMV * self.L  # m^3
        # Vapor molar volume (m^3/mol)
        self.VMV = 8.314*self.T/self.P
        self.VV = self.VMV * self.V  # m^3
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
    # calculate mass of oligomers
    self.OP2 = np.multiply(self.TCO, np.divide(self.TN, np.add(self.BACA, self.TN)))
    self.OP3 = np.multiply(self.OP2, np.divide(self.BACA, np.add(self.BACA, self.TN)))
    # calculate polymer (N6) mass
    self.Nylon = np.add(np.add(np.add(np.add(np.add(self.BACA, self.TN), self.TCO), self.TA), self.TCHA), self.P1)

    # return object
    return self