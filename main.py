from xml.sax.handler import feature_validation
import numpy as np
from cmath import pi
import matplotlib.pyplot as plt
from scipy.integrate import odeint

"""Initail Values"""
#region

V0 = 24 # Initial volume of fuel (gal)
roeg = 800 # density of gasoline (kg/m^3)
convV = 0.00378541 # m^3/gal
mf0 = V0*convV*roeg # Initial mass of fuel (kg)
m_base = 1800                  # Base mass of the car (kg)
m_cargo = 110                  # Mass of cargo in the car (kg)
vel = 70              # Velocity (mph)
Cd = 0.33             # Drag coefficient
roe = 1.2             # density of air at room temp (~70-75⁰F) (kg/m^3)
A = 2                 # Cross sectional area of vehicle (m^2)
g = 9.81                          # Gravitational acceleration
theta = 0                         # Driving slope angle (deg)
ur = 0.015 # coeff. rolling friction
e_e = 0.3           # engine efficiency
E = 46*(10**6)      # energy density (J/kg) for gasoline

# endregion

"""Differential Equation"""
def dmfdt(t,mf):
    mc_tot = m_base + m_cargo + mf # Total car mass
    conv = 2.23693629                    # mph/m/s
    v = vel/conv                         # Velocity (m/s)
    Fd = 0.5*Cd*roe*A*v**2               # Drag force
    N = mc_tot*g*np.cos((pi/180)*theta)  # Normal force
    Fg = mc_tot*g*np.sin((pi/180)*theta) # Horrizontal gravitaional force
    f = ur*N                             # friction due to tires
    F_neg = Fd + f + Fg                  # Total Negative Force
    Pe = F_neg*v                         # Power exerted by the car (W)

    return -Pe/(e_e*E) # Derivative of fuel mass with respect to time

"""Solving & Calculations"""
# region

### Time
t_total = 20 # Total time run (hrs)
t = np.linspace(0, int(t_total*3600), int(t_total*3600)) # Creates evenly spaced time values up to t_total
t_min = t/60 # Converts the time array into minutes

### Solution
sol = odeint(dmfdt, # Solves the given function
    y0=mf0,         # Initial y value
    t=t,            # Time array
    tfirst=True)    # Tells the solver the the x-variable comes first in the function
Vf = ((sol.T[0]/roeg)/convV) # Converts the solution array into volume (gal)
Vf = Vf[Vf >= 0]             # Removes negative values

### Mileage
mlg = [] # Initialize list

for i in range(1, len(Vf)):
    mlg_i = (vel*(1/3600)) / (Vf[i-1] - Vf[i]) # mileage = distance travelled / volume of fuel used
    mlg.append(mlg_i) # Add the value calculated above to the list

# endregion

"""Plotting"""
# region
plt.figure(num=1, figsize=(10, 10))
plt.title('Fuel Volume vs Time')
plt.xlabel('Time (min)')
plt.ylabel('Fuel Volume (gal)')
plt.scatter(t_min[0:len(Vf)], Vf, color='black')

plt.figure(num=2, figsize=(10, 10))
plt.title('Mileage vs Fuel Volume')
plt.xlabel('Fuel Volume (gal)')
plt.ylabel('Mileage (mpg)')
plt.scatter(Vf[0:-1], mlg, color='blue')
plt.show()
# endregion