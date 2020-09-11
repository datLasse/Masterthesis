import numpy as np
import matplotlib.pyplot as plt

#-----------------------------GENERAL CONSTANTS-----------------------------------------

boltzmann_constant = 8.62*10**(-5)
stefan_boltzmann_constant = 5.67*10**(-8)
c = 3e8
a = 7.56*10**(-16)


mass_solar = 2e30
radius_solar = 6.96e8
#------------------------------------PROGENITOR VALUES-------------------------------------------

# RSG --- THIS PROGENITOR HAS BEEN CHOSSEN AS IT IS KNOWN TO BE IN THERMAL EQUILIBIRUM
n = 3
kappa = 0.034
velo_index = 0.19

M_15 = 1.2
R_500 = 1.6
E_51 = 1.2
M = M_15*15*mass_solar
R =  R_500*500*radius_solar


#-------------------------ONLY A TEST FILE - - FUNCTIONS ARE DECLARED DIRECTLY IN FILE----------------------------------------------
def breakout_values(stellar_radius,stellar_mass,poly_index,opacity,velo_index,E_51,M_15):

    density_avg = stellar_mass/(stellar_radius**3)

    density_breakout = density_avg * np.power((c*(poly_index+1)/(opacity*stellar_radius*density_avg*1800e3))*np.power(M_15/E_51,0.5),poly_index/(poly_index*(1-velo_index)+1))

    radius_breakout = stellar_radius * (1 - np.power(density_breakout/density_avg,1/poly_index))

    mass_breakout = ((4*np.pi*stellar_radius**3*density_avg)/(poly_index+1))*np.power(density_breakout/density_avg,(poly_index+1)/poly_index)

    velocity_breakout = 1800e3*(E_51/M_15)**(1/2)*np.power(density_breakout/density_avg, -velo_index)

    time_breakout = (stellar_radius - radius_breakout) / velocity_breakout

    time_transition = stellar_radius / velocity_breakout

    energy_breakout = mass_breakout*np.power(velocity_breakout,2)

    coupling_breakout = 0.2 * np.power(velocity_breakout/(10**7), 15 / 4) * np.power(density_breakout/(10**(-6)), -1 / 8)

    opticaldepth_breakout = opacity*stellar_radius*density_avg/(poly_index+1)*np.power(density_breakout/density_avg,(poly_index+1)/poly_index)


    return radius_breakout, density_breakout, mass_breakout, time_breakout, velocity_breakout, time_transition, density_avg, energy_breakout,coupling_breakout, opticaldepth_breakout   #BREAKOUT VALUE CALCULATION IS UNCHANGED

#-------------ETA IS CALCULATED WITH THE DESNITY MODEL FOLLOWING RABINAK AND WAXMAN-------------------------------------------------------

def coupling_coefficentRW(mass_shell, mass_breakout, time, time_breakout, time_transition, poly_index, velocity_breakout, density_breakout, coupling_breakout):
    coupling_shell = 0

    if time == time_breakout:
        coupling_shell = coupling_breakout
    if time<time_transition:
        coupling_shell = (coupling_breakout * np.power(mass_shell / mass_breakout, -(17 * poly_index +16*velo_index*poly_index+ 25)/(8 * poly_index + 8)) *
         np.power(time / time_breakout, -1/6))
    return coupling_shell

#-------------ETA IS CALCULATED WITH NAKAR AND SARIS ASSUMED DENSITY---------------------------------------------------------------------------------------------

def coupling_coefficentNS(mass_shell, mass_breakout, time, time_breakout, time_transition, poly_index, velocity_breakout, density_breakout, coupling_breakout):
    coupling_shell = 0

    if time == time_breakout:
        coupling_shell = coupling_breakout
    if time<time_transition:
        coupling_shell = (coupling_breakout * np.power(mass_shell / mass_breakout, -(17 * poly_index + 9)/(8 * poly_index + 8)) *
         np.power(time / time_breakout, -1/6))
    return coupling_shell

#---------------------------THE CHANGES APPLIED TO THE DENSITY MODEL ONLY HOLD IN THE PLANAR PHASE AND THE CASE THAT ETA_0 < 1 -----------------------------
def temp_observableNS(time,time_breakout,time_transition,velocity_breakout,poly_index, mass_breakout, stellar_radius, velo_index,density_breakout,energy_breakout,coupling_breakout, depth_breakout,radius_breakout,temp_BB_breakout):
    mass_s = mass_breakout

    coupling_shell = coupling_coefficentNS(mass_s, mass_breakout, time, time_breakout, time_transition, poly_index, velocity_breakout, density_breakout,coupling_breakout)

    if time < time_transition:

        if coupling_shell < 1:#-------------------------------------------------------------------------------------------------------------------------------------------------------

            temp_obs = temp_BB_breakout*np.power(coupling_breakout,(2*(poly_index+1)/(17*poly_index+9)))*np.power(time/time_breakout,-(2*(9*poly_index+5))/(3*(17*poly_index+9)))
            return temp_obs
        else:
            print('$\eta>1$')

    else:
        print('$t>t_{s}$')

def temp_observableRW(time,time_breakout,time_transition,velocity_breakout,poly_index, mass_breakout, stellar_radius, velo_index,density_breakout,energy_breakout,coupling_breakout, depth_breakout,radius_breakout,temp_BB_breakout):
    mass_s = mass_breakout

    coupling_shell = coupling_coefficentRW(mass_s, mass_breakout, time, time_breakout, time_transition, poly_index, velocity_breakout, density_breakout,coupling_breakout)

    if time < time_transition:

        if coupling_shell < 1:#-------------------------------------------------------------------------------------------------------------------------------------------------------

            temp_obs = temp_BB_breakout*np.power(coupling_breakout,(2*(poly_index+1)/(17 * poly_index +16*velo_index*poly_index+ 25)))*np.power(time/time_breakout,-(2*(9*poly_index+13+8*velo_index*poly_index))/(3*(17 * poly_index +16*velo_index*poly_index+ 25)))
            return temp_obs
        else:
            print('$\eta>1$')

    else:
        print('$t>t_{s}$')

#--------------------------------CALCULATING VALUES ANALOG TO ANY OTHER FILE---------------------------------------------------------------------------------------------


r_0,rho_0,m_0,t_0,v_0,t_s,rho_star,E_0,eta_0,tau_0 = breakout_values(R,M,n,kappa,velo_index,E_51,M_15)

T_BB0 = ((E_0*tau_0)/(a*c*r_0**2*t_0))**(1/4)

print('Finished calculating Breakout Values')

runtime = int(t_s)
time = int(t_0)

tempNS = np.zeros(runtime-int(t_0))
tempRW = np.zeros(runtime-int(t_0))
print_time = np.zeros(runtime-int(t_0))

print('Assigned Memory for Arrays')

for time in range(int(t_0),runtime):

    print_time[time-int(t_0)] = time
    tempNS[time-int(t_0)] = temp_observableNS(time,t_0,t_s,v_0,n, m_0, R,velo_index, rho_0,E_0,eta_0,tau_0,r_0,T_BB0)
    tempRW[time-int(t_0)] = temp_observableRW(time,t_0,t_s,v_0,n, m_0, R,velo_index, rho_0,E_0,eta_0,tau_0,r_0,T_BB0)

print('Calculated Temperature')


plt.plot(print_time, tempNS,label = 'Approximated Density', color = 'red')
plt.plot(print_time, tempRW,label = 'Derived Density Density',color = 'green')
plt.legend()
plt.xlabel('Time [s]')
plt.ylabel('Temperature [K]')
plt.grid()
plt.ticklabel_format(axis='both', style='sci', scilimits = (0,0))


plt.show()
