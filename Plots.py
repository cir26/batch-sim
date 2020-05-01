import numpy as np
import matplotlib.pyplot as plt
from Reactor_Sim import Polymerization
from mpl_toolkits import mplot3d
from matplotlib import animation
from timeit import default_timer as timer
############################################################
# Nylon-6 production playground
############################################################

cap=1500
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
# 3D plots
def surface_plot(feed1,arr1,feed2,arr2,n):
    Z=[0]*n
    Z=[Z]*n
    Z=np.array(Z)
    for i,(x1,y1) in enumerate(zip(arr1,arr2)):
        for j, (x2,y2) in enumerate(zip(x1,y1)):
            state_dict[feed1]=x2
            state_dict[feed2] = y2
            units = 'kg'  # convert final units
            init_cond = [i for i in state_dict.values()]  # extract initial conditions for input
            Poly = Polymerization(273.15 + 255, init_cond, 10, ideal=False, P=1 * 101325, units=units)
            Z[i,j]=Poly.Nylon[-1]
    return np.asarray(Z)



n=10
X_label = 'W'
Y_label = 'CL'
Z_label = 'Nylon-6'
x = np.linspace(2,100, n)
y = np.linspace(100,2000, n)
X, Y = np.meshgrid(x, y)
print('X: ',X)
print('Y: ',Y)
start_time=timer()
Z=surface_plot(X_label,X,Y_label,Y,n)
print('Z: ',Z)
end_time=timer()
print("Minutes elapsed: ", (end_time - start_time)/60)
print("Run Time Rate: ", (end_time - start_time)/(n**2), " sec/iter")
#X=np.multiply(X,2.20462)   # convert to lbs
#Y=np.multiply(Y,2.20462)
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, Z/(X+Y), rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_title('Nylon-6 conversion with initial feeds')
ax.set_xlabel('{} (kg)'.format(X_label))
ax.set_ylabel('{} (kg)'.format(Y_label))
ax.set_zlabel('{} (% conversion)'.format(Z_label))
ax.view_init(30, 30)
plt.show()

conv = Z/(X+Y)

fig = plt.figure()
ax1 = fig.add_subplot(111)
colors = ['black', 'red', 'blue', 'green', 'orange', 'black', 'red', 'blue', 'green', 'orange']
alp = [0.5, 0.5, 0.5, 0.5, 0.5, 1 ,1 ,1 ,1 , 1]
for i in range(1,n,2):
    ax1.plot(X[i], conv[i], colors[i], alpha=alp[i], linestyle='dashdot', label='Initial CL: {}'.format(round(Y[i][0],0)))
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.set_ylabel('Nylon (% conversion)', color='k')
ax1.set_xlabel('Water (kg)')
ax1.legend()
plt.show()


fig = plt.figure()
ax1 = fig.add_subplot(111)
colors = ['black', 'red', 'blue', 'green', 'orange', 'black', 'red', 'blue', 'green', 'orange']
alp = [0.35, 0.35, 0.35, 0.35, 0.35, 1 ,1 ,1 ,1 , 1]
for i in range(1,n,2):
    ax1.plot(np.asarray([Y[j][i] for j in range(0,n)]), [conv[j][i] for j in range(0,n)], colors[i], alpha=alp[i], linestyle='dashdot', label='Initial W: {}'.format(round(X[0][i],0)))
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.set_ylabel('Nylon (% conversion)', color='k')
ax1.set_xlabel('Caprolactam (kg)')
ax1.legend()
plt.show()
#find_optimal(X_label,X,Y_label,Y,Z_label,Z)



# CHA & AA varying
cap=1500
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
n=10
X_label = 'AA'
Y_label = 'CHA'
Z_label = 'Nylon-6'
x = np.array([.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 15])
y = np.array([.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 15])
X, Y = np.meshgrid(x, y)
print('X: ',X)
print('Y: ',Y)
start_time=timer()
Z=surface_plot(X_label,X,Y_label,Y,n)
print('Z: ',Z)
end_time=timer()
print("Minutes elapsed: ", (end_time - start_time)/60)
print("Run Time Rate: ", (end_time - start_time)/(n**2), " sec/iter")
conv = Z/(cap+w)

#X=np.multiply(X,2.20462)   # convert to lbs
#Y=np.multiply(Y,2.20462)
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, conv, rstride=1, cstride=1,
                cmap='plasma', edgecolor='none')
ax.set_title('Nylon-6 conversion with varying terminators')
ax.set_xlabel('{} (kg)'.format(X_label))
ax.set_ylabel('{} (kg)'.format(Y_label))
ax.set_zlabel('{} (% conversion)'.format(Z_label))
ax.view_init(30, 30)
plt.show()


fig = plt.figure()
ax1 = fig.add_subplot(111)
colors = ['black', 'red', 'blue', 'green', 'orange', 'black', 'red', 'blue', 'green', 'orange']
alp = [0.5, 0.5, 0.5, 0.5, 0.5, 1 ,1 ,1 ,1 , 1]
for i in range(0,n,2):
    ax1.plot(X[i], conv[i], colors[i], alpha=alp[i], linestyle='dashdot', label='Initial CHA: {}'.format(round(Y[i][0],4)))
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.set_ylabel('Nylon (% conversion)', color='k')
ax1.set_xlabel('Acetic Acid (kg)')
ax1.legend()
plt.show()


fig = plt.figure()
ax1 = fig.add_subplot(111)
colors = ['black', 'red', 'blue', 'green', 'orange', 'black', 'red', 'blue', 'green', 'orange']
alp = [0.35, 0.35, 0.35, 0.35, 0.35, 1 ,1 ,1 ,1 , 1]
for i in range(0,n,2):
    ax1.plot(np.asarray([Y[j][i] for j in range(0,n)]), [conv[j][i] for j in range(0,n)], colors[i], alpha=alp[i], linestyle='dashdot', label='Initial AA: {}'.format(round(X[0][i],4)))
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.set_ylabel('Nylon (% conversion)', color='k')
ax1.set_xlabel('Cyclohexylamine (kg)')
ax1.legend()
plt.show()




# optimizing MW
cap=1500
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
def surface_plot(feed1,arr1,feed2,arr2,n):
    Z=[0]*n
    Z=[Z]*n
    Z=np.array(Z)
    for i,(x1,y1) in enumerate(zip(arr1,arr2)):
        for j, (x2,y2) in enumerate(zip(x1,y1)):
            state_dict[feed1]=x2
            state_dict[feed2] = y2
            units = 'kg'  # convert final units
            init_cond = [i for i in state_dict.values()]  # extract initial conditions for input
            Poly = Polymerization(273.15 + 255, init_cond, 10, ideal=False, P=1 * 101325, units=units)
            Z[i,j]=Poly.MWW[-100]
    return np.asarray(Z)


n=10
X_label = 'W'
Y_label = 'CL'
Z_label = 'Nylon-6'
x = np.linspace(2,100, n)
y = np.linspace(100,2000, n)
X, Y = np.meshgrid(x, y)
print('X: ',X)
print('Y: ',Y)
start_time=timer()
Z=surface_plot(X_label,X,Y_label,Y,n)
print('Z: ',Z)
end_time=timer()
print("Minutes elapsed: ", (end_time - start_time)/60)
print("Run Time Rate: ", (end_time - start_time)/(n**2), " sec/iter")
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_title('Nylon-6 MWW with initial feeds')
ax.set_xlabel('{} (kg)'.format(X_label))
ax.set_ylabel('{} (kg)'.format(Y_label))
ax.set_zlabel('{} MW (kg/mol)'.format(Z_label))
ax.view_init(30, 30)
plt.show()

conv = Z

fig = plt.figure()
ax1 = fig.add_subplot(111)
colors = ['black', 'red', 'blue', 'green', 'orange', 'black', 'red', 'blue', 'green', 'orange']
alp = [0.5, 0.5, 0.5, 0.5, 0.5, 1 ,1 ,1 ,1 , 1]
for i in range(1,n,2):
    ax1.plot(X[i], conv[i], colors[i], alpha=alp[i], linestyle='dashdot', label='Initial CL: {}'.format(round(Y[i][0],0)))
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.set_ylabel('Nylon (kg/mol)', color='k')
ax1.set_xlabel('Water (kg)')
ax1.legend()
plt.show()


fig = plt.figure()
ax1 = fig.add_subplot(111)
colors = ['black', 'red', 'blue', 'green', 'orange', 'black', 'red', 'blue', 'green', 'orange']
alp = [0.35, 0.35, 0.35, 0.35, 0.35, 1 ,1 ,1 ,1 , 1]
for i in range(1,n,2):
    ax1.plot(np.asarray([Y[j][i] for j in range(0,n)]), [conv[j][i] for j in range(0,n)], colors[i], alpha=alp[i], linestyle='dashdot', label='Initial W: {}'.format(round(X[0][i],0)))
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.set_ylabel('Nylon (kg/mol)', color='k')
ax1.set_xlabel('Caprolactam (kg)')
ax1.legend()
plt.show()
#find_optimal(X_label,X,Y_label,Y,Z_label,Z)



# CHA & AA varying
cap=1500
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
n=10
X_label = 'AA'
Y_label = 'CHA'
Z_label = 'Nylon-6'
x = np.array([.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 15])
y = np.array([.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 15])
X, Y = np.meshgrid(x, y)
print('X: ',X)
print('Y: ',Y)
start_time=timer()
Z=surface_plot(X_label,X,Y_label,Y,n)
print('Z: ',Z)
end_time=timer()
print("Minutes elapsed: ", (end_time - start_time)/60)
print("Run Time Rate: ", (end_time - start_time)/(n**2), " sec/iter")
#X=np.multiply(X,2.20462)   # convert to lbs
#Y=np.multiply(Y,2.20462)
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                cmap='plasma', edgecolor='none')
ax.set_title('Nylon-6 molecular weight with varying terminators')
ax.set_xlabel('{} (kg)'.format(X_label))
ax.set_ylabel('{} (kg)'.format(Y_label))
ax.set_zlabel('{} MW (kg/mol)'.format(Z_label))
ax.view_init(30, 30)
plt.show()

conv = Z

fig = plt.figure()
ax1 = fig.add_subplot(111)
colors = ['black', 'red', 'blue', 'green', 'orange', 'black', 'red', 'blue', 'green', 'orange']
alp = [0.5, 0.5, 0.5, 0.5, 0.5, 1 ,1 ,1 ,1 , 1]
for i in range(0,n,2):
    ax1.plot(X[i], conv[i], colors[i], alpha=alp[i], linestyle='dashdot', label='Initial CHA: {}'.format(round(Y[i][0],4)))
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.set_ylabel('Nylon (kg/mol)', color='k')
ax1.set_xlabel('Acetic Acid (kg)')
ax1.legend()
plt.show()


fig = plt.figure()
ax1 = fig.add_subplot(111)
colors = ['black', 'red', 'blue', 'green', 'orange', 'black', 'red', 'blue', 'green', 'orange']
alp = [0.35, 0.35, 0.35, 0.35, 0.35, 1 ,1 ,1 ,1 , 1]
for i in range(0,n,2):
    ax1.plot(np.asarray([Y[j][i] for j in range(0,n)]), [conv[j][i] for j in range(0,n)], colors[i], alpha=alp[i], linestyle='dashdot', label='Initial AA: {}'.format(round(X[0][i],4)))
ax1.minorticks_on()
ax1.tick_params(axis='x', which='minor', direction='in')
ax1.tick_params(axis='y', which='minor', direction='in')
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(direction="in")
ax1.set_ylabel('Nylon (kg/mol)', color='k')
ax1.set_xlabel('Cyclohexylamine (kg)')
ax1.legend()
plt.show()



