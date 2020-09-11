import numpy as np
import matplotlib.pyplot as plt
from function_declaration import *
import csv


mass_solar = 2e30
radius_solar = 6.96e8
luminosity_solar = 3.8e26
c = 3e8



datatime = np.zeros(60)
datapoints = np.zeros(60)

datatim = np.zeros(60)
datalum = np.zeros(60)


with open('Data/lum_7.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        datatim[i] = row[0]
        datalum[i] = row[1]
        i += 1


with open('Data/temp_7.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        datatime[i] = row[0]
        datapoints[i] = row[1]
        i += 1


n = 1.5
kappa = 0.034
mu = 0.19

R = 1.6*500*radius_solar
M = 1.2*15*mass_solar
M_15 = M/(15*mass_solar)
E_51 = 1.2

r_0,rho_0,m_0,t_0,v_0,t_s,rho_star,E_0,eta_0,tau_0 = breakout_values(R,M,n,kappa,mu,E_51,M_15)


print('Finished calculating Breakout Values')
print(tau_0,t_0,E_0)


runtime = int(1e6)
time = int(t_0)

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

plt.grid()
plt.loglog(print_time,luminosity_obs, label = 'Described Implementation')
plt.loglog(print_time,luminosity_obs/1.2, label = 'Adjusted Implementation', color = 'green')
plt.scatter(datatim,datalum*1.4,marker = 'x',color = 'black', label = 'Nakar and Sari (2010); Figure 6')
plt.ylabel('Luminosity [erg]')
plt.xlabel('Time [s]')
plt.legend()
#plt.xlim(left = 20,right = 30000)
plt.show()

plt.loglog(print_time, temp, label = 'Described Implementation', color = 'red')
plt.loglog(datatim,datapoints, label = 'Nakar and Sari (2010); Figure 7',color = 'green')
plt.legend()
plt.grid()
plt.ylabel('Temperature [K]')
plt.xlabel('Time [s]')
plt.xlim(left = 2e3,right = 1e6)
plt.show()
