import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

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

def Thermo_EquilB(T, P, zcap):
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

    zw=1-zcap

    ps1 = psat(A1, B1, C1, D1, E1, T)
    ps2 = psat(A2, B2, C2, D2, E2, T)
    g1_BP = gamma1(zw, zcap, G12, G21, t12, t21)
    g2_BP = gamma2(zw, zcap, G12, G21, t12, t21)

    ycap = zcap*g2_BP*ps2/P

    def f(xw,P,T):
        xcap = 1-xw
        g1_DP = gamma1(xw, 1-xw, G12, G21, t12, t21)
        g2_DP = gamma2(xw, 1-xw, G12, G21, t12, t21)
        return (xw*g1_DP*ps1)+(xcap*g2_DP*ps2)-P

    for i in range(100):

        xw, _, ier, _ = fsolve(f, i/100, args=(P, T), full_output=1)

        if xw > 1e-7 and xw < 1:
            print(xw,' ',i)
        else:
            pass
    return ycap, 1 - ycap, xw, 1 - xw
zcap=np.linspace(0,1,50)
T=np.linspace(60+273.15,120+273.15,50)
P=24000 # Pa
y_cap=[]
y_w=[]
x_cap=[]
x_w=[]
for i,j in zip(T,zcap):
    ycap,yw,xw,xcap= Thermo_EquilB(i,P,j)
    print(" ")
    y_cap.append(ycap)
    y_w.append(yw)
    x_cap.append(xcap)
    x_w.append(xw)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x_w, T, 'k', alpha = 0.85, linestyle='dashdot', label='x - Water')
ax1.plot(y_w, T, 'k', label='y - Water')
ax1.set_ylabel("Temperature (K)")
plt.legend()
plt.show()
