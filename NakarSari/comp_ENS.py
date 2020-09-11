import numpy as np
import matplotlib.pyplot as plt
from function_declaration import *
import csv


mass_solar = 2e30
radius_solar = 6.96e8
luminosity_solar = 3.8e26
c = 3e8



datatimeens1 = np.zeros(23)
datapointsens1 = np.zeros(23)

datatimeens2 = np.zeros(28)
datapointsens2 = np.zeros(28)

datatime5001 = np.zeros(13)
datapoints5001 = np.zeros(13)


datatime5002 = np.zeros(14)
datapoints5002 = np.zeros(14)

with open('Data/Ensman_5001.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        datatimeens1[i] = row[0]
        datapointsens1[i] = row[1]
        i += 1

with open('Data/Ensman_5002.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        datatimeens2[i] = row[0]
        datapointsens2[i] = row[1]
        i += 1


with open('Data/Ensman_5001_temp2.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        datatime5001[i] = row[0]
        datapoints5001[i] = row[1]
        i += 1

with open('Data/Ensman_5002_temp2.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        datatime5002[i] = row[0]
        datapoints5002[i] = row[1]
        i += 1




n = 3
kappa = 0.034
mu = 0.19

R = 45*radius_solar
M = 16*mass_solar
M_15 = M/(15*mass_solar)
E_51 = 2.3

acc = np.array([1000,10,1,0.1, 0.01, 0.001,0.0001,0.0001])


for i in range(8):

    r_0,rho_0,m_0,t_0,v_0,t_s,rho_star,E_0,eta_0,tau_0 = breakout_values(R,M,n,kappa,mu,E_51,M_15)


    print('Finished calculating Breakout Values')

    print(eta_0 )

    runtime = int(t_s+10)
    time = int(t_0)

    luminosity_obs = np.zeros(runtime-int(t_0))
    temp = np.zeros(runtime-int(t_0))
    print_time = np.zeros(runtime-int(t_0))


    print('Assigned Memory for Arrays')

    for time in range(int(t_0),runtime):

        print_time[time-int(t_0)] = time
        if time>int(t_0):
            temp[time-int(t_0)] = temp_observable(time,t_0,t_s,v_0,n, m_0, R,mu, rho_0,E_0,eta_0,tau_0,r_0,temp[time-int(t_0)-1],acc[i])
        else:
            temp[time-int(t_0)] = temp_observable(time,t_0,t_s,v_0,n, m_0, R,mu, rho_0,E_0,eta_0,tau_0,r_0,1,acc[i])
        luminosity_obs[time-int(t_0)] = luminosity(time,t_0,t_s,E_0,n)*10**7

    plt.loglog(print_time,temp, label = '$Accuracy = $' + str(acc[i]))

#plt.grid()

#plt.plot(datatime5001-6953,10**datapoints5001,color = 'red',linestyle = '--', label = 'Hydrodynamic Model for $E_{51} = 1$ by Ensman and Burrows (1992)')
#plt.plot(datatime5002-4635,10**datapoints5002,color = 'green',linestyle = '--', label = 'Hydrodynamic Model for $E_{51} = 2.3$ by Ensman and Burrows (1992)')


#plt.scatter(datatimeens1-6979,10**datapointsens1,color = 'black',marker = 'x', label = 'Hydrodynamic Model for $E_{51} = 1$ by Ensman and Burrows (1992)')
#plt.scatter(datatimeens2-4707,10**datapointsens2,color = 'black',marker = '.', label = 'Hydrodynamic Model for $E_{51} = 2.03$ by Ensman and Burrows (1992)')
plt.ylabel('Temperature [K]')
plt.xlabel('Time [s]')
plt.yscale('log')
#plt.xscale('log')
plt.legend()
#plt.xlim(left = 20, right = 2e4)
plt.show()
