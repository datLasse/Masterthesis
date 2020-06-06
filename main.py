#this file is meant to implement the relevent functions
import numpy as np
from function_declaration import *
import matplotlib.pyplot as plt

#importing .txt file with values such as stellar radius and stellar mass
input_ar = [[500,100,600,800,10,80,50,20],[15,7.5,25,40,10,200,150,70],[1,0.5,3,5,7,10,8,2],[1.5,1.5,1.5,1.5,3,3,3,3],[0.034,0.034,0.034,0.034,0.02,0.02,0.034,0.034]
,[0.19,0.19,0.19,0.19,0.19,0.19,0.19,0.19]]

solar_M = 2e30
solar_R = 6.96e8
solar_L = 3.8e26
c = 3e8
n_d = 20

for i in range(8):
    R = input_ar[0][i]*solar_R
    M = input_ar[1][i]*solar_M
    M_15 = M/(15*solar_M)
    E_51 = input_ar[2][i]
    n = input_ar[3][i]
    kappa = input_ar[4][i]
    mu = input_ar[5][i]

    r_0,rho_0,m_0,t_0,v_0,t_s,rho_star = breakout_values(R,M,n,kappa,mu,E_51,M_15)

    print("This is run" + str(i))

    if t_0 <1:
        t_0 = 1

    time = int(t_0)
    luminosity_obs = np.zeros(n_d*24*3600-int(t_0))
    temp = np.zeros(n_d*24*3600-int(t_0))
    print_time = np.zeros(n_d*24*3600-int(t_0))
    for time in range(int(t_0),n_d*24*3600):
        luminosity_obs[time-int(t_0)] = luminosity(time,t_0,t_s,m_0,v_0,n)
        print_time[time-int(t_0)] = time
        temp[time-int(t_0)] = temp_observable(time,t_0,t_s,v_0,n, m_0, R,mu, rho_0)




    fig, (ax1,ax2) = plt.subplots(1,2, constrained_layout=True)
    ax1.loglog(print_time,temp)
    ax1.set_title("Tempertaure")
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Temperature [K]")
    fig.suptitle('Plot for Parameter Set'+str(i), fontsize=16)
    ax2.loglog(print_time,luminosity_obs)
    ax2.set_title("Luminosity")
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Luminosity [W]")

    plt.savefig("Plots\plot"+str(i)+".png")

    plt.clf()




#first the breakout point and values need to be obtained


#iterating over time steps through both luminosity and the observable temperature


#converting luminsoisty to apparent magnitude and exporting a plot
