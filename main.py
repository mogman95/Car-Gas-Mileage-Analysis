from cmath import pi
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

"""Differential Equation"""
def dmfdt(t,mf):

    ### Mass
    m_base = 1800                  # Base mass of the car (kg)
    m_cargo = 110                  # Mass of cargo in the car (kg)
    mc_tot = m_base + m_cargo + mf # Total car mass
    
    ### Forces
    # Drag
    vel = 70              # Velocity (mph)
    conv = 2.23693629     # mph/m/s
    v = vel/conv          # Velocity (m/s)
    Cd = 0.33             # Drag coefficient
    roe = 1.2             # density of air at room temp (~70-75⁰F) (kg/m^3)
    A = 2                 # Cross sectional area of vehicle (m^2)
    Fd = 0.5*Cd*roe*A*v**2 # Drag force
    
    # Gravity
    g = 9.81                          # Gravitational acceleration
    theta = 0                         # Driving slope angle (deg)
    N = mc_tot*g*np.cos((pi/180)*theta)  # Normal force
    Fg = mc_tot*g*np.sin((pi/180)*theta) # Horrizontal gravitaional force
    
    # Friction
    ur = 0.01 # coeff. rolling friction
    f = ur*N  # friction due to tires
    
    F_neg = Fd + f + Fg # Total Negative Force
    
    ### Fuel Mass ODE
    e_e = 0.3           # engine efficiency
    E = 46*10^6         # energy density (J/kg) for gasoline 
    Pe = F_neg*v        # Power exerted by the car (W)

    return -Pe/(e_e*E) # Derivative of fuel mass with respect to time

"""Initail Values"""
# region

### Fuel
V0 = 10 # Initial volume of fuel (gal)
roeg = 800 # density of gasoline (kg/m^3)
convV = 0.00378541 # m^3/gal
mf0 = V0*convV*roeg # Initial mass of fuel (kg)


# endregion

"""Solving"""
# region

### Time
t_total = 1 # Total time run (hrs)
t = np.linspace(0, int(t_total*3600), int(t_total*3600)) # Creates evenly spaced time values up to t_total
t_min = t#/60
### Solver
sol = odeint(dmfdt, y0=mf0, t=t, tfirst=True)
#print(sol.T[0])

# endregion

"""Plotting"""
plt.figure(figsize=(10, 10))
plt.title('Fuel Mass vs Time')
plt.xlabel('Time (min)')
plt.ylabel('Fuel Mass (kg)')
plt.scatter(t_min,sol.T[0], color='blue')
#plt.scatter(x_d,y_d, color='red')
#plt.legend(['No Drag','With Drag','Terminal Velocity'])
plt.show()