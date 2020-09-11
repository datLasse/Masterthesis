import numpy as np
import matplotlib.pyplot as plt
from function_declaration import *
import csv


mass_solar = 2e30
radius_solar = 6.96e8
luminosity_solar = 3.8e26
c = 3e8



datatimeens1 = np.zeros(38)
datapointsens1 = np.zeros(38)

datatimeens2 = np.zeros(38)
datapointsens2 = np.zeros(38)

datatimetemp = np.zeros(64)
datapointstemp = np.zeros(64)

datatime14E13 = np.zeros(34)
datalum14E13 = np.zeros(34)

datatim14E07 = np.zeros(31)
datalum14E07 = np.zeros(31)

datatimeCTIO = np.zeros(11)
datapointsCTIO = np.zeros(11)

datatimeSAAO = np.zeros(8)
datapointsSAAO = np.zeros(8)


with open('Data/Luminosity_Blin_141.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        datatimeens1[i] = row[0]
        datapointsens1[i] = row[1]
        i += 1

with open('Data/Temp_Blin_141.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        datatimetemp[i] = row[0]
        datapointstemp[i] = row[1]

        i += 1

with open('Data/lum_blin_2.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        datatimeens2[i] = row[0]
        datapointsens2[i] = row[1]
        i += 1

with open('Data/Blin_14E1.3.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        datatime14E13[i] = row[0]
        datalum14E13[i] = row[1]
        i += 1

with open('Data/Blin_14E07.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        datatim14E07[i] = row[0]
        datalum14E07[i] = row[1]
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



n = 3
kappa = 0.034
mu = 0.19

R = np.array([40*radius_solar,45*radius_solar,50*radius_solar])
M = np.array([28*mass_solar,25*mass_solar,16.3*mass_solar])
M_15 = M/(15*mass_solar)
E_51 = np.array([0.7,1,1.3])


for i in range(3):

    r_0,rho_0,m_0,t_0,v_0,t_s,rho_star,E_0,eta_0,tau_0 = breakout_values(R[i],M[i],n,kappa,mu,E_51[i],M_15[i])

    print('Finished calculating Breakout Values')


    print(eta_0)

    runtime = int(10*3600*24)
    time = int(t_0)

    luminosity_obs = np.zeros(runtime-int(t_0))
    temp = np.zeros(runtime-int(t_0))
    print_time = np.zeros(runtime-int(t_0))


    print('Assigned Memory for Arrays')

    for time in range(int(t_0),runtime):

        print_time[time-int(t_0)] = time

        if time>int(t_0):
            temp[time-int(t_0)] = temp_observable(time,t_0,t_s,v_0,n, m_0, R[i],mu, rho_0,E_0,eta_0,tau_0,r_0,temp[time-int(t_0)-1])
        else:
         temp[time-int(t_0)] = temp_observable(time,t_0,t_s,v_0,n, m_0, R[i],mu, rho_0,E_0,eta_0,tau_0,r_0,1)

        luminosity_obs[time-int(t_0)] = luminosity(time,t_0,t_s,E_0,n)*10**7

    plt.loglog(print_time ,temp, label = 'Parameter Set E'+str(E_51[i]) )





plt.grid()


plt.scatter(datatimetemp - 6627,datapointstemp, label = "Model 14E1 by Blinnikov et al. (2000)",color='black',marker = '.')

#plt.plot(datatimeens1-6583,10**datapointsens1,color = 'black', label = 'Hydrodynamic Model 14E1 by Blinnikov et al. (2000)')
#plt.plot(datatimeens2,10**datapointsens2,color = 'black', label = 'Hydrodynamic Model E14E1 by Blinnikov et al. (2000)')
#plt.plot(datatime14E13,10**datalum14E13, color = 'black',linestyle = '-.',label = 'Hydrodynamic Model E14E1.3 by Blinnikov et al. (2000)')
#plt.plot(datatim14E07,10**datalum14E07,color = 'black',linestyle = '--', label = 'Hydrodynamic Model E14E0.7 by Blinnikov et al. (2000)')


#plt.scatter(datatimeCTIO,10**datapointsCTIO, label = "Observational Data from CTIO",color='black',marker = '.')
#plt.scatter(datatimeSAAO,10**datapointsSAAO, label = "Observational Data from SAAO",color = 'black',marker = 'x')
#plt.axvspan(t_0, t_s, color='grey', alpha=0.25)
#plt.axvspan(t_s, runtime, color='green', alpha=0.25)
plt.ylabel('Temperature [K]')
plt.xlabel('Time [s]')
plt.yscale('log')

plt.legend()
#plt.xlim(left = 20,right = 30000)
plt.show()
