import numpy as np

# function to perform addition and calculation of important data attributes to Polymerization object
def attr_calc(self, state_arr, time_arr):
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
    # calculate number-average degree of polymerization
    num = np.add(np.add(np.add(self.BACA, self.TN), self.TCO), self.P1)
    denom = np.add(self.TCO, self.P1)
    #denom = self.TCO
    self.DP = np.divide(num, denom)
            #W = [0]
            #CL = [1]
            #CD = [2]
            #AA = [3]
            #P1 = [4]
            #BACA =[5]
            #TN = [6]
            #TCO = [7]
            #TA = [8]
            #CHA = [9]
            #TCHA = [10]
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
    # calculate number-average molecular weight of polymer (N6)
    num = np.add(np.add(np.add(np.add(np.add(self.molar_masses[6]*self.TN,
                                             self.molar_masses[7]*self.TCO),
                                             self.molar_masses[5]*self.BACA),
                                             self.molar_masses[4]*self.P1),
                                             self.molar_masses[8]*self.TA),
                                             self.molar_masses[10]*self.TCHA)
    denom = np.multiply(np.add(np.add(np.add(self.TN,self.TCO),self.TA),self.TCHA),0.5)
    self.MWN = np.divide(num,denom)
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