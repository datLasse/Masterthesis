#This file is meant to store the declarations of all used functions


import numpy as np


#--------------------------------------------------------------------------------------------------

#the breakout point and the relevant values are calculated based on the Appendix A

def breakout_values(stellar_radius,stellar_mass,poly_index,opacity,velo_index,E_51,M_15):

    density_avg = stellar_mass/(stellar_radius**3)

    density_breakout = density_avg * np.power((c*(poly_index+1)/(opacity*stellar_radius*density_avg*1800e3))*np.power(M_15/E_51,0.5),poly_index/(poly_index*(1-velo_index)+1))

    radius_breakout = stellar_radius * (1 - np.power(density_breakout/density_avg,1/poly_index))

    mass_breakout = ((4*np.pi*stellar_radius**3*density_avg)/(poly_index+1))*np.power(density_breakout/density_avg,(poly_index+1)/poly_index)

    velocity_breakout = 1800e3*(E_51/M_15)**(1/2)*np.power(density_breakout/density_avg, -velo_index)

    time_breakout = (stellar_radius - radius_breakout) / velocity_breakout

    time_transition = stellar_radius / velocity_breakout

    thickness = (stellar_radius - radius_breakout)/stellar_radius

    return radius_breakout, density_breakout, mass_breakout, time_breakout, velocity_breakout, time_transition, density_avg, thickness

print(breakout_values(R,M,n,kappa,mu,E_51,M_15))

#---------------------------------------------------------------------------------------------------

#the mass of the luminosity shell is calculated for each specific given moment in time compare eq.6

def mass_shell(mass_breakout,time,time_transition,poly_index,velo_index):
    if time>time_transition:

        mass_shell = np.power(time/time_transition
        ,(2*(poly_index+1))
        /((1+velo_index)*poly_index+1))*mass_breakout #hhh

        return mass_shell

    else:

        mass_shell = mass_breakout

        return mass_shell

#----------------------------------------------------------------------------------------------------

#the luminosity of the luminosity shell is calculated based on eq.4

def luminosity(time,time_breakout,time_transition,mass_breakout,velocity_breakout,poly_index):
    if time>time_transition:

        luminosity_obs = (mass_breakout*velocity_breakout**2)
        /(time_breakout)
        *np.power(time_transition/time_breakout
        ,-4/3)
        *np.power(time/time_transition
        ,-(2.28*poly_index-2)
        /(3*(1.19*poly_index+1)))

        return luminosity_obs

    else:

        luminosity_obs = (mass_breakout*velocity_breakout**2)
        /(time_breakout)
        *np.power(time/time_breakout
        ,-4/3)

        return luminosity_obs

#------------------------------------------------------------------------------------------------------

#the thermal coupling coefficent referred to as eta in the paper is calculated based on eq.10,16,18

def coupling_coefficent(mass_shell, mass_breakout, time, time_breakout, time_transition, poly_index, velocity_breakout, density_breakout):
    coupling_breakout = 0, 2 * np.power(velocity_breakout / 10 ** 4, 15 / 4) * np.power(density_breakout / 10 ** (-9), -1 / 8)
    coupling_shell = 0

    if time<time_transition:
        coupling_shell = eta_breakout * np.power(mass_shell / mass_breakout, -(17 * poly_index + 9) / (8 * poly_index + 8)) * np.power(time / time_breakout, -1 / 6)
    if time>time_transition:
        coupling_shell = eta_breakout * np.power(mass_shell / mass_breakout, -(22.32 * poly_index + 17) / (8 * poly_index + 8)) * np.power(time_transition / time_breakout, -1 / 6) * np.power(time / time_transition, (42 * poly_index + 49) / (12 * (1, 19 * poly_index + 1)))
    return coupling_shell

#----------------------------------------------------------------------------------------------------------

#the coefficent for inverse thomson scattering (referred to as xi) is calculated by first finding the relevant blackbody temperature based
# on the stefan-boltzmann law and the luminosity then calculating the temperature if inverse compton scattering wouldn't make a contribution
# and checking if this is indeed the case if so the coefficent and the temperature is returned. If inverse thomspn scattering makes a non-neglible contribution
# the numercial solver scipy fsolve is used to calculate the temperature which then is used to calculate the thompson coefficent. Based on eq. 13.
#This function also returns the calculated blackbody temperature.

def thompson_coefficent(luminosity_obs, stellar_radius,coupling_shell,density_shell):
    stefan_boltzmann_constant = 5.67*10**(-8)
    boltzmann_constant = 8.62*10**(-5)
    temp_BB = np.power(luminosity_obs/(4*np.pi*stellar_radius**2*stefan_boltzmann_constant),1/4)
    temp_breakout = coupling_shell**2*temp_BB
    yMAX = 3*np.power(density_shell/(10**(-6)),-1/2)*np.power(temp_breakout*boltzmann_constant/100)

    if yMAX<1:
        thompson_coefficent = 1
        return thompson_coefficent, temp_breakout, temp_BB
    if 1>0.5*np.log(yMAX)*(1.6+np.log(yMAX)):
        thompson_coefficent = 1
        return thompson_coefficent, temp_breakout, temp_BB
    else:
        condition = 1
        while condition == 1:
            temp_new = coupling_shell**2*temp_BB*np.power(0.5*np.log(3*np.power(density_shell/(10**(-6)),-1/2)*np.power(temp_breakout*boltzmann_constant/100))*(1.6+np.log(3*np.power(density_shell/(10**(-6)),-1/2)*np.power(temp_breakout*boltzmann_constant/100))),-2)
            if ((temp_new-temp_breakout)**2)**(0.5)<0.01:
                condition = 0
                temp_breakout = temp_new
            else:
                temp_breakout = temp_new
        thompson_coefficent = np.power(0.5*np.log(3*np.power(density_shell/(10**(-6)),-1/2)*np.power(temp_breakout*boltzmann_constant/100))*(1.6+np.log(3*np.power(density_shell/(10**(-6)),-1/2)*np.power(temp_breakout*boltzmann_constant/100))),-2)
        return thompson_coefficent, temp_breakout, temp_BB

#-----------------------------------------------------------------------------------------------------------------------------

def temp_observable(temp_BB,coupling_shell,thompson_shell,time,time_breakout,time_transition):
    return temp_obs


def density_shell(time, time_breakout, time_transition, velocity_breakout, radius_breakout, ):
    if time<time_transition:
        

    return densisty
