import numpy as np
import pandas as pd
from cmath import pi
import matplotlib.pyplot as plt
from scipy.integrate import odeint

"""Initail Values"""
#region

V0 = 10                    # Initial volume of fuel (gal)
roeg = 800                 # density of gasoline (kg/m^3)
convV = 0.00378541         # m^3/gal
mf0 = V0*convV*roeg        # Initial mass of fuel (kg)
m_base = 1800              # Base mass of the car (kg)
m_cargo = 110              # Mass of cargo in the car (kg)
vel = np.arange(65,76,0.5) # Velocity (mph)
Cd = 0.33                  # Drag coefficient
roe = 1.2                  # density of air at room temp (~70-75⁰F) (kg/m^3)
A = 2                      # Cross sectional area of vehicle (m^2)
g = 9.81                   # Gravitational acceleration
theta = 0                  # Driving slope angle (deg)
ur = 0.015                 # coeff. rolling friction
e_e = 0.3                  # engine efficiency
E = 46*(10**6)             # energy density (J/kg) for gasoline

# endregion

"""Differential Equation"""
def dmfdt(t,mf):
    mc_tot = m_base + m_cargo + mf # Total car mass
    conv = 2.23693629                    # mph/m/s
    v = vel_i/conv                         # Velocity (m/s)
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
t_total = 10 # Total time run (hrs)
t = np.linspace(0, int(t_total*3600), int(t_total*3600 + 1)) # Creates evenly spaced time values up to t_total
t_min = t/60 # Converts the time array into minutes

Vf_v = pd.DataFrame(index=t, columns=vel)
mlg_v = pd.DataFrame(index=vel, columns=['Average Mileage'])
# print(mlg_v)

for vel_i in vel:
    ### Solution
    sol_i = odeint(dmfdt,          # Solves the given function
        y0=mf0,         # Initial y value
        t=t,            # Time array
        tfirst=True)    # Tells the solver the the x-variable comes first in the function
    Vf = ((sol_i.T[0]/roeg)/convV) # Converts the solution array into volume (gal)
    Vf = Vf[Vf >= 0]               # Removes negative values
    # for i,value in enumerate(Vf):
    #     Vf_v.loc[i,vel_i] = value  # Add the Vf values to the dataframe
    # print(vel_i,1)
    t_fin = t[len(Vf)-1] # Finds time when car runs out of gas

    avg_mlg = (vel_i * (t_fin/3600)) / V0 # Finds average mileage
    mlg_v.loc[vel_i] = avg_mlg
    ### Mileage
    # mlg = [] # Initialize list
    # for i in range(0, len(Vf)-1):
    #     mlg_i = (vel_i*(1/3600)) / (Vf[i] - Vf[i+1]) # mileage = distance travelled / volume of fuel used
    #     mlg_v.loc[i,vel_i] = mlg_i                   # Add the value calculated above to the dataframe
    # print(vel_i,2)
print(mlg_v)

# endregion

"""Plotting"""
# region

plot = mlg_v.plot(title = 'Average Mileage vs Velocity',
    xlabel = 'Velocity (mph)',
    ylabel = 'Average Mileage (mpg)',
    grid=True)
plt.show()

# endregion