import numpy as np
import matplotlib as plt
import scipy as sp
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

"""Differential Equation"""
def dSdt(t,mf):

    ### Mass
    m_base = 1800                  # Base mass of the car (kg)
    m_cargo = 150                  # Mass of cargo in the car (kg)
    mc_tot = m_base + m_cargo + mf # Total car mass
    
    ### Forces
    # Drag
    vel = 70              # Velocity (mph)
    conv = 2.23693629     # mph/m/s
    v = vel/conv          # Velocity (m/s)
    Cd = 0.33             # Drag coefficient
    roe = 1.2             # density of air at room temp (~70-75⁰F) (kg/m^3)
    A = 2                 # Cross sectional area of vehicle (m^2)
    Fd = 0.5*Cd*roe*A*v^2 # Drag force
    
    # Gravity
    g = 9.81                          # Gravitational acceleration
    theta = 2                         # Driving slope angle (deg)
    N = mc_tot*g*cos((pi/180)*theta)  # normal force
    Fg = mc_tot*g*sin((pi/180)*theta) # Horrizontal gravitaional force
    
    # Friction
    ur = 0.01 # coeff. rolling friction
    f = ur*N  # friction due to tires
    
    F_neg = Fd + f + Fg # Total Negative Force
    
    ### Fuel Mass ODE
    e_e = 0.3           # engine efficiency
    E = 46*10^6         # energy density (J/kg) for gasoline 
    Pe = F_neg*v        # Power exerted by the car
    dmfdt = -Pe/(e_e*E) # Derivative of fuel mass with respect to time

    return [dmfdt]


"""Car Properties"""
### Fuel
V0 = 10 # Initial volume of fuel (gal)
roeg = 800 # density of gasoline (kg/m^3)
convV = 0.00378541 # m^3/gal
mf0 = V0*convV*roeg # Initial mass of fuel (kg)

