import numpy as np
from Thermo import Thermo_EquilB
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


def dif_rxn(Reactions):
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

	dW = R2 + R3 + R4 + R5 + R11 + R12 - (R1 + R8)
	dCL = -(R1 + R6 + R7)
	dCD = -(R8 + R9 + R10)
	dAA = -(R11 + R12)
	dP1 = R1 - (2*R2 + R3 + R4 + R6 + R9 + R11)
	dBACA = R3 + R4 + 2*R5 + R7 + R9 + 2*R10 + R12
	dTN = R2 + R6 + R8 + R9 - (R5 + R12)
	dTC = R2 + R6 + R8 + R9 + R11 - R5
	dTA = R11 + R12
	
	return [dW, dCL, dCD, dAA, dP1, dBACA, dTN, dTC, dTA]	

rxn_const = {
				1: [1.66e2, 8.32e4, 1.2e4, 7.87e4, 8.03e3, -33.01],
				2: [5.26e6, 9.74e4, 3.37e6, 8.65e4, -2.49e4, 3.951],
				3: [7.93e5, 9.56e4, 4.55e6, 8.42e4, -1.69e4, -29.08],
				4: [2.38e8, 1.76e5, 6.47e8, 1.57e5, -4.02e4, -60.79],
				5: [7.14e4, 8.92e4, 8.36e5, 8.54e4, -1.33e4, 2.439]
			}


def kinetic_const(T, TC):
	rxn_const = {
					1: [1.66e2, 8.32e4, 1.2e4, 7.87e4, 8.03e3, -33.01],
					2: [5.26e6, 9.74e4, 3.37e6, 8.65e4, -2.49e4, 3.951],
					3: [7.93e5, 9.56e4, 4.55e6, 8.42e4, -1.69e4, -29.08],
					4: [2.38e8, 1.76e5, 6.47e8, 1.57e5, -4.02e4, -60.79],
					5: [7.14e4, 8.92e4, 8.36e5, 8.54e4, -1.33e4, 2.439]
				}

	def k_func(T, R, TC, rxn_const, flag):
		a0 = rxn_const[flag][0]
		e0 = rxn_const[flag][1]
		ac = rxn_const[flag][2]
		ec = rxn_const[flag][3]
		
		return a0*np.exp(-e0/R/T)+ac*np.exp(-ec/R/T)*TC
		
	def K_func(T, R, rxn_const, flag):
		
		H = rxn_const[flag][4]
		S = rxn_const[flag][5]
		
		return np.exp((S-H/T)/R)

	R = 8.314
	
	lil_k = []
	K = []
	for i in range(1,6):
		lil_k.append(k_func(T, R, TC, rxn_const, i))
		K.append(K_func(T, R, rxn_const, i))
	return lil_k, K
			
def rxn(conc, k, K):
	
	k1 = k[0]
	k2 = k[1]
	k3 = k[2]
	k4 = k[3]
	k5 = k[4]
	
	K1 = K[0]
	K2 = K[1]
	K3 = K[2]
	K4 = K[3]
	K5 = K[4]
	
	W = conc[0]
	CL = conc[1]
	CD = conc[2]
	AA = conc[3]
	P1 = conc[4]
	BACA = conc[5]
	TN = conc[6]
	TC = conc[7]
	TA = conc[8]



	if BACA != 0 or TN !=0:

		P2 = TC*(BACA/(BACA+TN))
		P3 = TC*BACA*TN/(BACA+TN)/(BACA+TN)

	else:

		P2 = 0
		P3 = 0
																																																																																													
	def R_1(CL, W, P1, k1, K1):
		
		return k1*CL*W-k1/K1*P1
		
	def R_2(P1, P2, W, k2, K2):   
		
		return k2*P1*P1 - k2/K2*P2*W
		
	def R_3(P1, TC, W, BACA, TN, k2, K2):

		if BACA != 0 or TN != 0:
			return k2*P1*TC-k2/K2*W*TC*(BACA/(BACA+TN))
		else:
			return k2*P1*TC
		
	def R_4(TN, P1, W, BACA, TC, k2, K2):

		if BACA != 0 or TC != 0:
			return k2*TN*P1-k2/K2*W*TN*(BACA/(BACA+TC))
		else:
			return k2*TN*P1
		
	def R_5(TN, TC, W, BACA, k2, K2):

		if BACA != 0 or TN != 0:
			return k2*TN*TC-k2/K2*W*BACA*(BACA/(BACA+TN))
		else:
			return k2*TN*TC
		
	def R_6(P1, CL, P2, k3, K3):
		
		return k3*P1*CL-k3/K3*P2
		
	def R_7(TN, CL, BACA, TC, k3, K3):
		if BACA != 0 or TN != 0:
			return k3*TN*CL-k3/K3*TN*(BACA/(BACA+TN))
		else:
			return k3*TN*CL
		
	def R_8(CD, W, P2, k4, K4):
		
		return k4*CD*W-k4/K4*P2
		
	def R_9(P1, CD, P3, k5, K5):
		
		return k5*P1*CD-k5/K5*P3
		
	def R_10(TN, CD, BACA, TC, k5, K5):

		if BACA != 0 or TC != 0:
			return k5*TN*CD-k5/K5*TN*(BACA/(BACA+TC))
		else:
			return k5*TN*CD
		
	def R_11(AA, P1, W, TA, TC, BACA, k2, K2):

		if BACA !=0 or TC != 0:
			return k2*AA*P1-k2/K2*W*TA*(TC/(TC+BACA))
		else:
			return k2*AA*P1
		
	def R_12(AA, TN, W, TA, BACA, TC, k2, K2):

		if BACA != 0 or TC != 0:
			return k2*AA*TN-k2/K2*W*TA*(BACA/(TC+BACA))
		else:
			return k2*AA*TN
	
	R1 = R_1(CL, W, P1, k1, K1)
	R2 = R_2(P1, P2, W, k2, K2)
	R3 = R_3(P1, TC, W, BACA, TN, k2, K2)
	R4 = R_4(TN, P1, W, BACA, TC, k2, K2)
	R5 = R_5(TN, TC, W, BACA, k2, K2)
	R6 = R_6(P1, CL, P2, k3, K3)
	R7 = R_7(TN, CL, BACA, TC, k3, K3)
	R8 = R_8(CD, W, P2, k4, K4)
	R9 = R_9(P1, CD, P3, k5, K5)
	R10 = R_10(TN, CD, BACA, TC, k5, K5)
	R11 = R_11(AA, P1, W, TA, TC, BACA, k2, K2)
	R12 = R_12(AA, TN, W, TA, BACA, TC, k2, K2)
	
	return [R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12]


def rxn_BT(conc, k, K, T, P, max_iter):
	k1 = k[0]
	k2 = k[1]
	k3 = k[2]
	k4 = k[3]
	k5 = k[4]

	K1 = K[0]
	K2 = K[1]
	K3 = K[2]
	K4 = K[3]
	K5 = K[4]

	W = conc[0]
	CL = conc[1]
	CD = conc[2]
	AA = conc[3]
	P1 = conc[4]
	BACA = conc[5]
	TN = conc[6]
	TC = conc[7]
	TA = conc[8]

	zw = W/(W+CL)
	zcap = CL/(W+CL)
	F = W + CL

	#Adjust Water and Caprolactam Concentrations with NRTL
	V, L, F, xw, xcap, yw, ycap, g1, g2 = Thermo_EquilB(T, P, F, zw, zcap, 100)
	W = L[0]*xw[0]
	CL = L[0]*xcap[0]



	if BACA != 0 or TN != 0:

		P2 = TC * (BACA / (BACA + TN))
		P3 = TC * BACA * TN / (BACA + TN) / (BACA + TN)

	else:

		P2 = 0
		P3 = 0

	def R_1(CL, W, P1, k1, K1):

		return k1 * CL * W - k1 / K1 * P1

	def R_2(P1, P2, W, k2, K2):

		return k2 * P1 * P1 - k2 / K2 * P2 * W

	def R_3(P1, TC, W, BACA, TN, k2, K2):

		if BACA != 0 or TN != 0:
			return k2 * P1 * TC - k2 / K2 * W * TC * (BACA / (BACA + TN))
		else:
			return k2 * P1 * TC

	def R_4(TN, P1, W, BACA, TC, k2, K2):

		if BACA != 0 or TC != 0:
			return k2 * TN * P1 - k2 / K2 * W * TN * (BACA / (BACA + TC))
		else:
			return k2 * TN * P1

	def R_5(TN, TC, W, BACA, k2, K2):

		if BACA != 0 or TN != 0:
			return k2 * TN * TC - k2 / K2 * W * BACA * (BACA / (BACA + TN))
		else:
			return k2 * TN * TC

	def R_6(P1, CL, P2, k3, K3):

		return k3 * P1 * CL - k3 / K3 * P2

	def R_7(TN, CL, BACA, TC, k3, K3):
		if BACA != 0 or TN != 0:
			return k3 * TN * CL - k3 / K3 * TN * (BACA / (BACA + TN))
		else:
			return k3 * TN * CL

	def R_8(CD, W, P2, k4, K4):

		return k4 * CD * W - k4 / K4 * P2

	def R_9(P1, CD, P3, k5, K5):

		return k5 * P1 * CD - k5 / K5 * P3

	def R_10(TN, CD, BACA, TC, k5, K5):

		if BACA != 0 or TC != 0:
			return k5 * TN * CD - k5 / K5 * TN * (BACA / (BACA + TC))
		else:
			return k5 * TN * CD

	def R_11(AA, P1, W, TA, TC, BACA, k2, K2):

		if BACA != 0 or TC != 0:
			return k2 * AA * P1 - k2 / K2 * W * TA * (TC / (TC + BACA))
		else:
			return k2 * AA * P1

	def R_12(AA, TN, W, TA, BACA, TC, k2, K2):

		if BACA != 0 or TC != 0:
			return k2 * AA * TN - k2 / K2 * W * TA * (BACA / (TC + BACA))
		else:
			return k2 * AA * TN

	R1 = R_1(CL, W, P1, k1, K1)
	R2 = R_2(P1, P2, W, k2, K2)
	R3 = R_3(P1, TC, W, BACA, TN, k2, K2)
	R4 = R_4(TN, P1, W, BACA, TC, k2, K2)
	R5 = R_5(TN, TC, W, BACA, k2, K2)
	R6 = R_6(P1, CL, P2, k3, K3)
	R7 = R_7(TN, CL, BACA, TC, k3, K3)
	R8 = R_8(CD, W, P2, k4, K4)
	R9 = R_9(P1, CD, P3, k5, K5)
	R10 = R_10(TN, CD, BACA, TC, k5, K5)
	R11 = R_11(AA, P1, W, TA, TC, BACA, k2, K2)
	R12 = R_12(AA, TN, W, TA, BACA, TC, k2, K2)

	return [R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12]

def main():
	pass









