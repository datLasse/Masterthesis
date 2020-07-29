import numpy as np
import matplotlib.pyplot as plt
from function_declaration import *
import csv


datatime = np.zeros(60)
datapoints = np.zeros(60)
datatim = np.zeros(60)
datalum = np.zeros(60)
datatimeCTIO = np.zeros(11)
datapointsCTIO = np.zeros(11)
datatimeSAAO = np.zeros(8)
datapointsSAAO = np.zeros(8)


with open('Data/temp_7.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        datatime[i] = row[0]
        datapoints[i] = row[1]
        i += 1

with open('Data/lum_7.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        datatim[i] = row[0]
        datalum[i] = row[1]
        i += 1

with open('Data/Sn1987_CTIO.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        datatimeCTIO[i] = row[0]
        datapointsCTIO[i] = row[1]
        i += 1

with open('Data/Sn1987_SAAO.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        datatimeSAAO[i] = row[0]
        datapointsSAAO[i] = row[1]
        i += 1


mass_solar = 2e30
radius_solar = 6.96e8
luminosity_solar = 3.8e26
c = 3e8

R = 1.6*500*radius_solar
M = 1.2*15*mass_solar
M_15 = M/(15*mass_solar)
E_51 = 1.2
n = 1.5
kappa = 0.034
mu = 0.19

r_0,rho_0,m_0,t_0,v_0,t_s,rho_star,E_0,eta_0,tau_0 = breakout_values(R,M,n,kappa,mu,E_51,M_15)

print('Finished calculating Breakout Values')

#print(E_0)

#r_0 = 3.11*10**(10)
#tau_0 =9.34
#t_0 = 5.032
#m_0 = 3.37e24
#v_0 =31044.83*10**3
#t_s = 1064.96
#rho_0 = 1.534e-6
#E_0 = 2.09*10**41
#eta_0 = 12.3

time_1 = (t_s
        *np.power(eta_0*np.power(t_0/t_s,1/6),3*(1.19*n+1)/(9.88*n+5)))
time2 = (t_s
        *np.power(eta_0*np.power(t_0/t_s,1/6),6*(1.19*n+1)/(12.48*n+1)))

runtime = int(10*3600*24)

#print(runtime)
time = int(t_0)
print(t_0)

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
    #print(temp[time-int(t_0)-1])
    luminosity_obs[time-int(t_0)] = luminosity(time,t_0,t_s,E_0,n)*10**7

print('Calculated Temperature and Luminosity')



plt.plot(print_time,temp, label = "My code")
plt.yscale('log')
plt.plot(datatime,datapoints, label = "Nakar and Saris Fig.7")
plt.legend()
plt.xscale('log')
#plt.axvspan(t_0, t_s, color='grey', alpha=0.25)
#plt.axvspan(t_s, runtime, color='green', alpha=0.25)
#plt.axvspan(time_1, time2, color='red', alpha=0.25)
#plt.axvspan(time2, runtime, color='yellow', alpha=0.25)
plt.xlabel("Time [s]")
plt.ylabel("Temperature [K]")
#plt.xlim(left = 100)
#plt.ylim(10**4,1e6)
#plt.xlim(right = 10**6)
plt.grid()
plt.xscale('log')

plt.show()
plt.clf()

print('Plotted Temperature')

plt.plot(print_time,luminosity_obs, label = "My code")
plt.plot(datatim,datalum, label = "Nakar and Saris Fig.7")
#plt.scatter(datatimeCTIO,10**datapointsCTIO, label = "Observational Data from CTIO",color='black',marker = '.')
#plt.scatter(datatimeSAAO,10**datapointsSAAO, label = "Observational Data from SAAO",color = 'black',marker = 'x')
#plt.xlim(right = 140)
plt.yscale('log')
plt.xscale('log')
plt.xlabel("Time [s]")
plt.ylabel("Luminosity [erg]")
plt.legend()
#plt.xlim(left= 0,right = 140)
plt.grid()
plt.show()

print('Plotted Luminosity')
