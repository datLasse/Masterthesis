import numpy as np
import matplotlib.pyplot as plt
from function_declaration import *
import csv


datatime = np.zeros(59)
datapoints = np.zeros(59)
datatim = np.zeros(59)
datalum = np.zeros(59)

with open('temp_6.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        datatime[i] = row[0]
        datapoints[i] = row[1]
        i += 1

with open('lum_6.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        datatim[i] = row[0]
        datalum[i] = row[1]
        i += 1

mass_solar = 2e30
radius_solar = 6.96e8
luminosity_solar = 3.8e26
c = 3e8

R = 45*radius_solar
M = 16*mass_solar
M_15 = M/(15*mass_solar)
E_51 = 2.3
n = 3
kappa = 0.034
mu = 0.19

r_0,rho_0,m_0,t_0,v_0,t_s,rho_star,E_0,eta_0 = breakout_values(R,M,n,kappa,mu,E_51,M_15)

r_0 = 3.11*10**(10)
tau_0 =9.34
t_0 = 5.032
m_0 = 3.37e24
v_0 =31044.83*10**3
t_s = 1064.96
rho_0 = 1.534e-6
E_0 = 4.25e39
eta_0 = 12.3

time = int(t_0)
luminosity_obs = np.zeros(3600*24-int(t_0))
temp = np.zeros(3600*24-int(t_0))
print_time = np.zeros(3600*24-int(t_0))
for time in range(int(t_0),3600*24):
    print_time[time-int(t_0)] = time
    temp[time-int(t_0)] = temp_observable(time,t_0,t_s,v_0,n, m_0, R,mu, rho_0,E_0,eta_0,tau_0,r_0)
    luminosity_obs[time-int(t_0)] = luminosity(time,t_0,t_s,E_0,n)*10**(7)



print(r_0,
rho_0,
m_0,
t_0
,v_0
,t_s
,rho_star
,E_0)

plt.loglog(print_time,temp, label = "me")
plt.loglog(datatime,datapoints, label = "N+S")
plt.legend()
plt.xlabel("Time [s]")
plt.ylabel("Temperature [K]")
#plt.ylim(10**4,6e6)
#plt.xlim(right = 3e4)

plt.savefig("Plots\luminosity.png")
plt.clf()

plt.loglog(print_time,luminosity_obs, label = "me")
plt.loglog(datatim,datalum, label = "nakar")
plt.legend()
#plt.xlim(left= 2e1,right = 3e4)
plt.savefig("Plots\Badn.png")
