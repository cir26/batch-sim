import numpy as np
from Thermo import Thermo_EquilB



def dif_rxn(condition, max_T, T_rate, t, P_rate=0, min_P=None, P_delay=0, term_delay=0): # returns list of reaction interaction differential equation for each species
	T = condition[11]
	if min_P is None:
		pass
	else:
		P = condition[12]
	k, K = kinetic_const(T, condition[7])

	if min_P is None:
		Reactions = rxn(condition, k, K)
	else:
		Reactions = rxn(condition, k, K, T, P, ideal=False, max_iter=100)

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

	def forcing(t, delay):
		if t >= delay:
			return 1
		else:
			return 0
	dW = R2 + R3 + R4 + R5 + R11 + R12 - (R1 + R8)
	dCL = -(R1 + R6 + R7)
	dCD = -(R8 + R9 + R10)
	dAA = -(R11 + R12)*forcing(t,term_delay)
	dP1 = R1 - (2*R2 + R3 + R4 + R6 + R9 + R11)
	dBACA = R3 + R4 + 2*R5 + R7 + R9 + 2*R10 + R12
	dTN = R2 + R6 + R8 + R9 - (R5 + R12)
	dTC = R2 + R6 + R8 + R9 + R11 - R5
	dTA = R11 + R12*forcing(t,term_delay)
	dCHA = -(R13 + R14 + R15)
	dTCHA = R13 + R14 + R15
	dT = T_rate*T*(max_T-T)
	if min_P==None:
		return [dW, dCL, dCD, dAA, dP1, dBACA, dTN, dTC, dTA, dCHA, dTCHA, dT]
	else:
		dP = P_rate * P * (1 - (P / min_P))*forcing(t,P_delay)
		return [dW, dCL, dCD, dAA, dP1, dBACA, dTN, dTC, dTA, dCHA, dTCHA, dT, dP]


def kinetic_const(T, TC): # returns list of forward rxn constant, and list of reaction constant equilibrium (K=k/k_prime)
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


def rxn(conc, k, K, T=None, P=None, ideal=True, max_iter=None):
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

	for i in range(len(conc)):
		if conc[i] < 0:
			conc[i] = 0
		else:
			pass

	# mol spec / kg tot mix
	W = conc[0]
	CL = conc[1]
	CD = conc[2]
	AA = conc[3]
	P1 = conc[4]
	BACA = conc[5]
	TN = conc[6]
	TC = conc[7]
	TA = conc[8]
	CHA = conc[9]
	TCHA = conc[10]

	if ideal is True:
		yw = [0]
		ycap = [0]
		V = [0]
		xw = [W / (W + CL)]
		xcap = [CL / (W + CL)]
		L = [W + CL]
	else:
		zw = W / (W + CL)
		zcap = CL / (W + CL)
		F = W + CL  # FEED
		# Adjust Water and Caprolactam Concentrations with NRTL
		V, L, F, xw, xcap, yw, ycap, g1, g2, _, _ = Thermo_EquilB(T, P, F, zw, zcap, max_iter)
		W = L[0] * xw[0]
		CL = L[0] * xcap[0]



	def R_1(CL, W, P1, k1, K1):

		return k1 * CL * W - k1 / K1 * P1

	def R_2(P1, W, TN, TC, BACA, TCHA, k2, K2):

		if any([TC != 0, BACA != 0, TCHA != 0]):
			return k2 * P1 * P1 - k2 / K2 * W * TN * TC / (TC + BACA + TCHA)
		else:
			return k2 * P1 * P1

	def R_3(P1, TC, W, BACA, TN, TA, k2, K2):

		if any([BACA != 0, TN != 0, TA != 0]):
			return k2 * P1 * TC - k2 / K2 * W * TC * BACA / (BACA + TN + TA)
		else:
			return k2 * P1 * TC

	def R_4(TN, P1, W, BACA, TC, TCHA, k2, K2):

		if any([BACA != 0, TC != 0, TCHA != 0]):
			return k2 * TN * P1 - k2 / K2 * W * TN * BACA / (BACA + TC + TCHA)
		else:
			return k2 * TN * P1

	def R_5(TN, TC, W, BACA, TA, k2, K2):

		if any([BACA != 0, TN != 0, TA != 0]):
			return k2 * TN * TC - k2 / K2 * W * BACA * BACA / (BACA + TN + TA)
		else:
			return k2 * TN * TC

	def R_6(P1, CL, TN, TC, BACA, TCHA, k3, K3):
		if any([TC != 0, BACA != 0, TCHA != 0]):
			return k3 * P1 * CL - k3 / K3 * TN * TC / (TC + BACA + TCHA)
		else:
			return k3 * P1 * CL

	def R_7(TN, CL, BACA, TC, TCHA, k3, K3):
		if any([BACA != 0, TC != 0, TCHA != 0]):
			return k3 * TN * CL - k3 / K3 * TN * BACA / (BACA + TC + TCHA)
		else:
			return k3 * TN * CL

	def R_8(CD, W, TN, TC, BACA, TCHA, k4, K4):
		if any([TC != 0, BACA != 0, TCHA != 0]):
			return k4 * CD * W - k4 / K4 * TN * TC / (TC + BACA + TCHA)
		else:
			return k4 * CD * W

	def R_9(P1, CD, TN, BACA, TC, TCHA, k5, K5):
		if any([BACA != 0, TC != 0, TCHA != 0]):
			return k5 * P1 * CD - k5 / K5 * TN * BACA * TC / (BACA + TC + TCHA) / (BACA + TC + TCHA)
		else:
			return k5 * P1 * CD

	def R_10(TN, CD, BACA, TC, TCHA, k5, K5):
		if any([BACA != 0, TC != 0, TCHA != 0]):
			return k5 * TN * CD - k5 / K5 * TN * (BACA / (BACA + TC + TCHA)) ** 2
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

	def R_13(CHA, CL, TCHA, TN, BACA, k3, K3):

		if TN != 0 or BACA != 0:
			return k3 * CHA * CL - k3 / K3 * TCHA * (TN / (TN + BACA))
		else:
			return k3 * CHA * CL

	def R_14(CHA, P1, W, TCHA, TN, BACA, k2, K2):

		if TN != 0 or BACA != 0:
			return k2 * CHA * P1 - k2 / K2 * W * TCHA * (TN / (TN + BACA))
		else:
			return k2 * CHA * P1

	def R_15(CHA, TC, W, TCHA, BACA, TN, k2, K2):

		if TN != 0 or BACA != 0:
			return k2 * CHA * TC - k2 / K2 * W * TCHA * (BACA / (TN + BACA))

		else:
			return k2 * CHA * TC

	R1 = R_1(CL, W, P1, k1, K1)
	R2 = R_2(P1, W, TN, TC, BACA, TCHA, k2, K2)
	R3 = R_3(P1, TC, W, BACA, TN, TA, k2, K2)
	R4 = R_4(TN, P1, W, BACA, TC, TCHA, k2, K2)
	R5 = R_5(TN, TC, W, BACA, TA, k2, K2)
	R6 = R_6(P1, CL, TN, TC, BACA, TCHA, k3, K3)
	R7 = R_7(TN, CL, BACA, TC, TCHA, k3, K3)
	R8 = R_8(CD, W, TN, TC, BACA, TCHA, k4, K4)
	R9 = R_9(P1, CD, TN, BACA, TC, TCHA, k5, K5)
	R10 = R_10(TN, CD, BACA, TC, TCHA, k5, K5)
	R11 = R_11(AA, P1, W, TA, TC, BACA, k2, K2)
	R12 = R_12(AA, TN, W, TA, BACA, TC, k2, K2)
	R13 = R_13(CHA, CL, TCHA, TN, BACA, k3, K3)
	R14 = R_14(CHA, P1, W, TCHA, TN, BACA, k2, K2)
	R15 = R_15(CHA, TC, W, TCHA, BACA, TN, k2, K2)

	return [R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15]

def main():
	pass









