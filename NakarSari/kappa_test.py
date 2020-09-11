import numpy as np
import matplotlib.pyplot as plt
from function_declaration import *
import csv


mass_solar = 2e30
radius_solar = 6.96e8
luminosity_solar = 3.8e26
c = 3e8


n = 3
kappa = np.array([0.001,0.005,0.01,0.015,0.02])
mu = 0.19
kappat = np.array([0.01,0.05,0.1,0.15,0.2])
acc = np.array([1,0.1, 0.01, 0.001,0.0001,0.00001])


M_15 = 1
R_5 = 1
E_51 = 1
M = M_15*15*mass_solar
R =  R_5*5*radius_solar


for i in range(5):

    r_0,rho_0,m_0,t_0,v_0,t_s,rho_star,E_0,eta_0,tau_0 = breakout_values(R,M,n,kappa[i],mu,E_51,M_15)

    print('Finished calculating Breakout Values')

    print(r_0,rho_0,m_0,t_0,v_0,t_s,rho_star,E_0,eta_0,tau_0 )

    runtime = int(2*24*3600)
    if t_0<1:
        time = 1
    else:
        time = int(t_0)
    print(time)
    luminosity_obs = np.zeros(runtime-int(t_0))
    temp = np.zeros(runtime-int(t_0))
    print_time = np.zeros(runtime-int(t_0))


    print('Assigned Memory for Arrays')

    for time in range(int(t_0),runtime):

        print_time[time-int(t_0)] = time
        if time>int(t_0):
            temp[time-int(t_0)] = temp_observable(time,t_0,t_s,v_0,n, m_0, R,mu, rho_0,E_0,eta_0,tau_0,r_0,temp[time-int(t_0)-1])
        else:
            temp[time-int(t_0)] = temp_observable(time,t_0,t_s,v_0,n, m_0, R,mu, rho_0,E_0,eta_0,tau_0,r_0,1)
        luminosity_obs[time-int(t_0)] = luminosity(time,t_0,t_s,E_0,n)*10**7


    plt.loglog(print_time,temp, label = '$\kappa = $' + str(kappat[i]) +'$cm^2g^{-1}$')

plt.ylabel('Temperature [K]')
plt.xlabel('Time [s]')
plt.yscale('log')
plt.legend()
plt.grid()
plt.show()
