Authors: Robert McMillin (mcmillinre3) and Cristian Romero (cir26)

A class named Polymerization is developed to assist in simulating Nylon-6 production.

This project is split into 6 modules including:
Enthalpy - functions defined for enthalpy calculations.
Thermo - functions defined for thermodynamic phase equilibrium calculations.
Kinetics - functions defined for reaction kinetics calculations.
Calculations - Post-simulation calculations can be found here (molecular weight averages, polymer mass, system net enthalpy, etc).
Reactor_Sim - MAIN module calling all other modules used for setting reaction conditions and running simulation with suggested plots.
Plots - other plots explored can be found here.


To run a simulation execute Reactor_Sim.py.
To change reactor conditions edit Reactor_Sim.py starting at line 145.
Suggested plots are included in Reactor_Sim starting at line 177.

Example simulation run:
Simulation_Object = Polymerization(Temperature=273.15+255,
                                   Initial_Charge=init_cond,
                                   End_Time=10,
                                   ideal=False,
                                   P=1*101325,
                                   units=units)

Class arguments explained:
Temperature - Set system final temperature in Kelvin. Temperature ramp is currently artificially imposed on system to occur within first hour.
Initial_Charge - Set initial conditions of system as list. This includes initial reactant charges in kg, temperature in Kelvin, and Pressure in Pascal.
                 Order is important here. Use the included dictionary named state_dict (line 152) as a guide.
                 state_dict = {'W':w,                # Water
                               'CL':cap,             # Caprolactam
                               'CD':0,               # Cyclic dimer
                               'AA':0,               # Acetic acid
                               'P1':0,               # Polymer, terminated 1 chain length (P2 and P3 are also calculated)
                               'BACA':0,             # Bounded n6 polymer
                               'TN':0,               # Terminating amine group
                               'TC':0,               # Terminating carboxylic acid
                               'TA':0,               # Terminating acetic acid
                               'CHA':cha,            # Cyclohexylamine
                               'TCHA':0,             # Terminating cyclohexylamine
                               'Temp':273.15+90,     # Temperature in Kelvin
                               'Press':5*101325}     # Pressure in Pascal
End_Time - Set duration of reaction in hours
ideal - Set to False. Originally used to test reaction with ideal liquid system but this function is now deprecated.
P - Final system pressure in Pascal. Pressure ramp is currently artificially imposed on system to occur between first and second hour.
units - Set to 'kg'. The other option is 'lbs' but suggested plots are not guaranteed to be produced with specified units.