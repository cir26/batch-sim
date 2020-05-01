import numpy as np
from scipy.optimize import fsolve

# functions used for phase equilbira calculations
def psat(A, B, C, D, E, T):
    return np.exp(A + B / T + C * np.log(T) + D * T ** E)


def gamma1(x1, x2, g12, g21, t12, t21):

    gamma = np.exp((x2) * (x2) * (t21 * (g21 / (x1 + (x2) * g21)) ** 2 + t12 * g12 / (((x2) + x1 * g12) ** 2)))

    return gamma


def gamma2(x1, x2, g12, g21, t12, t21):
    gamma = np.exp(x1 * x1 * (t12 * (g12 / ((x2) + x1 * g12)) ** 2 + t21 * g21 / ((x1 + (x2) * g21) ** 2)))

    return gamma


def tau(a, b, c, T):
    return a + b / T + c * np.log(T)


def VF_Ratio(V, z1, z2, K1, K2, F):
    a = V / F
    k1 = 1 / (K1 - 1)
    k2 = 1 / (K2 - 1)
    val = z1 / (k1 + a) + z2 / (k2 + a)
    return val


def xy(V, F, z1, z2, K1, K2):
    a = V / F
    x1 = z1 / (1 + a * (K1 - 1))
    x2 = z2 / (1 + a * (K2 - 1))
    y1 = K1 * x1
    y2 = K2 * x2
    return x1, x2, y1, y2


def Thermo_EquilB(T, P, F, zw, zcap, max_iter):
    a12 = -0.313
    a21 = 0.628
    b12 = -15.4
    b21 = -13.7
    c12 = 0.0495
    c21 = -0.0898

    t12 = tau(a12, b12, c12, T)
    t21 = tau(a21, b21, c21, T)

    G12 = np.exp(-.3 * t12)
    G21 = np.exp(-.3 * t21)

    A1 = 7.3649e1
    A2 = 7.4172e1

    B1 = -7.2582e3
    B2 = -1.0469e4

    C1 = -7.3037
    C2 = -6.8944

    D1 = 4.1653e-6
    D2 = 1.2113e-18

    E1 = 2.
    E2 = 6.

    ps1 = psat(A1, B1, C1, D1, E1, T)
    ps2 = psat(A2, B2, C2, D2, E2, T)

    g1_BP = gamma1(zw, zcap, G12, G21, t12, t21)
    g2_BP = gamma2(zw, zcap, G12, G21, t12, t21)

    ideal_DP = 1 / (zw / ps1 + zcap / ps2)
    ideal_BP = zw * ps1 * g1_BP + zcap * ps2 * g2_BP

    if P > ideal_BP:

        V = 0
        L = F
        return [V], [L], [F], [zw], [zcap], [0], [0], [g1_BP], [g2_BP], t12, t21

    elif P < ideal_DP:
        V = F
        L = 0
        return [V], [L], [F], [0], [0], [zw], [zcap], [1], [1], t12, t21

    else:
        g1_guess = 1
        g2_guess = 1

    guess_K1 = ps1 * g1_guess / P
    guess_K2 = ps2 * g2_guess / P

    x_old = 1
    for i in range(max_iter):

        V, _, ier, _ = fsolve(VF_Ratio, 0.25, args=(zw, zcap, guess_K1, guess_K2, F), full_output=1)

        if ier != 1:

            V, _, ier, _ = fsolve(VF_Ratio, 0.5, args=(zw, zcap, guess_K1, guess_K2, F), full_output=1)

            if ier != 1:

                V, _, ier, _ = fsolve(VF_Ratio, 0.75, args=(zw, zcap, guess_K1, guess_K2, F), full_output=1)

            else:
                V, _, ier, _ = fsolve(VF_Ratio, 0.75, args=(zw, zcap, guess_K1, guess_K2, F), full_output=1)
                print('Maximum Function Calls Reached')


        xw, xcap, yw, ycap = xy(V, F, zw, zcap, guess_K1, guess_K2)
        g1_guess = gamma1(xw, xcap, G12, G21, t12, t21)
        g2_guess = gamma2(xw, xcap, G12, G21, t12, t21)

        guess_K1 = ps1 * g1_guess / P
        guess_K2 = ps2 * g2_guess / P
        err = np.fabs(x_old - xw)

        if err < 1e-8:
            L = F - V
            return V, L, F, xw, xcap, yw, ycap, g1_guess, g2_guess, t12, t21
        x_old = xw

    L = F - V
    return V, L, F, xw, xcap, yw, ycap, g1_guess, g2_guess, t12, t21
