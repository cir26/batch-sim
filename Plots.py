import numpy as np
import matplotlib.pyplot as plt
from Reactor_Sim import Polymerization
from mpl_toolkits import mplot3d
from matplotlib import animation
from timeit import default_timer as timer
############################################################
# Nylon-6 production playground
############################################################

state_dict={'W':12.9375, # input mass in kg
            'CL':500,
            'CD':0,
            'AA':1e-6,
            'P1':0,
            'BACA':0,
            'TN':0,
            'TC':0,
            'TA':0,
            'CHA':0,
            'TCHA':0,
            'Temp':273.15+90,
            'Press':5*101325}
# units='lbs'
# state = [i for i in state_dict.values()] # extract initial conditions for input
# Poly = Polymerization([273.15+255, 1.4e-6], state, 10, ideal=False, P=1*101325, units=units)
# Poly2 = Polymerization([273.15+255, 1.4e-6], state, 10, ideal=True, units=units)
#
# # plot N6 production, temp, and press vs time
# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# ax1.plot(Poly.t, Poly.Nylon, 'k', alpha = 0.85, linestyle='dashdot', label='NRTL Binary Nylon-6 ({})'.format(units))
# ax1.plot(Poly2.t, Poly2.Nylon, 'k', label='Ideal Nylon-6 ({})'.format(units))
# ax1.plot(Poly.t, Poly.T,'r',label='Temperature (K)')
# ax1.plot(Poly.t, Poly.P,'b',label='Pressure (atm e-2)')
# ax1.minorticks_on()
# ax1.tick_params(axis='x', which='minor', direction='in')
# ax1.tick_params(axis='y', which='minor', direction='in')
# ax1.yaxis.set_ticks_position('both')
# ax1.xaxis.set_ticks_position('both')
# ax1.tick_params(direction="in")
# ax1.legend()
# ax1.set_xlabel('Time (hours)')
# plt.show()

# 3D plots
def surface_plot(feed1,arr1,feed2,arr2,n):
    Z=[0]*n
    Z=[Z]*n
    Z=np.array(Z)
    for i,(x1,y1) in enumerate(zip(arr1,arr2)):
        for j, (x2,y2) in enumerate(zip(x1,y1)):
            state_dict[feed1]=x2
            state_dict[feed2] = y2
            units='lbs'
            state = [i for i in state_dict.values()] # extract initial conditions for input
            Poly = Polymerization([273.15+255, 1.4e-6], state, 10, ideal=False, P=1*101325, units=units)
            Z[i,j]=Poly.Nylon[-1]
    return Z

def find_optimal(metric='max',feed1,arr1,feed2,arr2,opt,arr3):
    if metric == 'max':
        max = 0
        optimal = {feed1: 0, feed2: 0, opt:0}
        for (a1, b1, c1) in zip(arr1, arr2, arr3):
            for (a2, b2, c2) in zip(a1,b1,c1):
                if c2 > max:
                    max = c2
                    optimal[feed1] = a2
                    optimal[feed2] = b2
                    optimal[opt] = c2
    else:
        print("Metric not available")
    print('Optimal settings \n {}: {} \n {}: {} \n {}: {} \n'.format(feed1, optimal[feed1],
                                                                     feed2, optimal[feed2],
                                                                     opt, optimal[opt]))

n=10
X_label = 'W'
Y_label = 'CL'
Z_label = 'Nylon-6'
x = np.linspace(5,20, n)
y = np.linspace(10,600, n)

X, Y = np.meshgrid(x, y)
start_time=timer()
Z=surface_plot(X_label,X,Y_label,Y,n)
end_time=timer()
print("Minutes elapsed: ", (end_time - start_time)/60)
print("Run Time Rate: ", (end_time - start_time)/(n**2), " sec/iter")
#X=np.multiply(X,2.20462)   # convert to lbs
#Y=np.multiply(Y,2.20462)
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_title('Nylon-6 production with initial feeds')
ax.set_xlabel('{} (kg)'.format(X_label))
ax.set_ylabel('{} (kg)'.format(Y_label))
ax.set_zlabel('{} (kg)'.format(Z_label))
ax.view_init(25, 30)
plt.show()

find_optimal(X_label,X,Y_label,Y,Z_label,Z)
