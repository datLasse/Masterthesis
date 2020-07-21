#this file is meant to implement the relevent functions
import numpy as np
from function_declaration import *
import matplotlib.pyplot as plt


#importing .txt file with values such as stellar radius and stellar mass
input_ar = [[500,100,600,800,10,80,50,20],[15,7.5,25,40,10,200,150,70],[1,0.5,3,5,7,10,8,2],[1.5,1.5,1.5,1.5,3,3,3,3],[0.034,0.034,0.034,0.034,0.02,0.02,0.034,0.034]
,[0.19,0.19,0.19,0.19,0.19,0.19,0.19,0.19]]


mass_solar = 2e30
radius_solar = 6.96e8
luminosity_solar = 3.8e26
c = 3e8


for i in range(1):
    R = 45*radius_solar
    M = 16*mass_solar
    M_15 = M/(15*mass_solar)
    E_51 = 2.3
    n = 3
    kappa = 0.02
    mu = 0.19

    r_0,rho_0,m_0,t_0,v_0,t_s,rho_star = breakout_values(R,M,n,kappa,mu,E_51,M_15)

    print("This is run" + str(i))

    if t_0 <1:
        t_0 = 1

    time = int(t_0)
    luminosity_obs = np.zeros(10**4-int(t_0))
    temp = np.zeros(10**4-int(t_0))
#    Mag = np.zeros(n_d*24*3600-int(t_0))
    print_time = np.zeros(10**4-int(t_0))
    for time in range(int(t_0),10**4):
        print_time[time-int(t_0)] = time
        temp[time-int(t_0)] = temp_observable(time,t_0,t_s,v_0,n, m_0, R,mu, rho_0)
        luminosity_obs[time-int(t_0)] = luminosity(time,t_0,t_s,m_0,v_0,n)




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
#    ax2.plot(print_time,Mag)
#    plt.xscale("log")
#    ax2.set_title("Absolute Magnitude")
#    ax2.set_xlabel("Time [s]")
#    ax2.set_ylabel("Magnitude")

    plt.savefig("Plots\plot_"+str(i)+".png")

    plt.clf()
#
#    plt.plot(print_time,Mag)
#    plt.xscale("log")
#    plt.title("Magnitude Parameter Set"+str(i))
#    plt.xlabel("Time [s]")
#    plt.ylabel("Magnitude")

#    plt.savefig("Plots\Magnitude_plot_"+str(i)+".png")




#first the breakout point and values need to be obtained


#iterating over time steps through both luminosity and the observable temperature


#converting luminsoisty to apparent magnitude and exporting a plot
