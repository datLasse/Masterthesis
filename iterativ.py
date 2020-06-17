import numpy as np
import function_declaration as fd


radius_solar = 6.96e8
R = 45*radius_solar
c = 3e8

r_0 = 3.11*10**(10)
tau_0 =9.34
t_0 = 5.032
m_0 = 3.37e24
v_0 =31044.83*10**3
t_s = 1064.96
rho_0 = 1.534e-6
E_0 = 4.25e39
eta_0 = 12.3
n = 3

def thompson_coefficent(temp_obs,density_breakout):
    boltzmann_constant = 8.617*10**(-5)
    density_unit = density_breakout*10**(-3)
    temp_unit = temp_obs*boltzmann_constant
    yMAX = 3*np.power(density_unit/10e-9,-1/2)*np.power(temp_unit/100,9/4)
    thompson = (1/2)*np.log(yMAX)*(1.6+np.log(yMAX))

    if thompson < 1:
        thompson = 1
        return thompson

    else:
        thompson = 1
        return thompson

def temp_observable_breakout(coupling_breakout,time_breakout,time_transition,energy_breakout,poly_index,density_breakout,depth_breakout,radius_breakout):
    stefan_boltzmann_constant = 5.67e-8
    temp_BB = (np.power((fd.luminosity(time_breakout,time_breakout,time_transition,energy_breakout,poly_index)*depth_breakout)
            /(radius_breakout**2*stefan_boltzmann_constant*c),
            1/4))

    temp_obs = coupling_breakout**2*temp_BB
    itcount = 0
    terminator = 1

    while terminator == 1:

        if thompson_coefficent(temp_obs,density_breakout) == 1:

            terminator = 0
            print("thompson coefficent = 1 | " + str(itcount)+ " | "+str(temp_obs))
            return temp_obs

        else:

            temp_iter = temp_BB*coupling_breakout**2*thompson_coefficent(temp_obs,density_breakout)**(-2)

            if np.abs(temp_iter-temp_obs) < 0.1:

                temp_obs = temp_iter
                terminator = 0
                print("solved by iteration | " + str(itcount)+ " | " +str(temp_obs))
                return temp_obs

            else:

                temp_obs = temp_iter
                itcount += 1


temp = temp_observable_breakout(eta_0,t_0,t_s,E_0,n,rho_0,tau_0,r_0)

print(temp)
