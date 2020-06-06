#This file is meant to store the declarations of all used functions


import numpy as np


c = 3e8
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


    return radius_breakout, density_breakout, mass_breakout, time_breakout, velocity_breakout, time_transition, density_avg



#---------------------------------------------------------------------------------------------------

#the mass of the luminosity shell is calculated for each specific given moment in time compare eq.6

def mass_shell(mass_breakout,time,time_transition,poly_index,velo_index):
    if time>time_transition:

        mass_shell = np.power(time/time_transition
        ,(2*(poly_index+1))
        /((1+velo_index)*poly_index+1))*mass_breakout

        return mass_shell

    else:

        mass_shell = mass_breakout

        return mass_shell

#----------------------------------------------------------------------------------------------------

#the luminosity of the luminosity shell is calculated based on eq.4

def luminosity(time,time_breakout,time_transition,mass_breakout,velocity_breakout,poly_index):
    if time>time_transition:

        luminosity_obs = ((mass_breakout*velocity_breakout**2)
        /(time_breakout)
        *np.power(time_transition/time_breakout
        ,-4/3)
        *np.power(time/time_transition
        ,-(2.28*poly_index-2)
        /(3*(1.19*poly_index+1))))

        return luminosity_obs

    else:

        luminosity_obs = ((mass_breakout*velocity_breakout**2)
        /(time_breakout)
        *np.power(time/time_breakout
        ,-4/3))

        return luminosity_obs

#------------------------------------------------------------------------------------------------------

#the thermal coupling coefficent referred to as eta in the paper is calculated based on eq.10,16,18

def coupling_coefficent(mass_shell, mass_breakout, time, time_breakout, time_transition, poly_index, velocity_breakout, density_breakout):
    coupling_breakout = 0.2 * np.power(velocity_breakout*10**-7, 15 / 4) * np.power(density_breakout*10**(6), -1/8)
    coupling_shell = 0

    if time == time_breakout:
        coupling_shell = coupling_breakout
    if time<time_transition:
        coupling_shell = (coupling_breakout * np.power(mass_shell / mass_breakout, -(17 * poly_index + 9)/(8 * poly_index + 8)) *
         np.power(time / time_breakout, -1/6))
    if time>time_transition:
        coupling_shell = coupling_breakout * np.power(mass_shell / mass_breakout, -(22.32 * poly_index + 17)/(8 * poly_index + 8)) * np.power(time_transition / time_breakout, -1 / 6) * np.power(time / time_transition, (42 * poly_index + 49) / (12 * (1.19 * poly_index + 1)))
    return coupling_shell

#-----------------------------------------------------------------------------------------------------------------------

def temp_blackbody(stellar_radius,time,time_breakout,time_transition,mass_breakout,velocity_breakout,poly_index):
    stefan_boltzmann_constant = 5.67*10**(-8)
    temp_BB = np.power(luminosity(time,time_breakout,time_transition,mass_breakout,velocity_breakout,poly_index)/(4*np.pi*stellar_radius**2*stefan_boltzmann_constant),1/4)
    return temp_BB


#----------------------------------------------------------------------------------------------------------

#the coefficent for inverse thomson scattering (referred to as xi) is calculated by first finding the relevant blackbody temperature based
# on the stefan-boltzmann law and the luminosity then calculating the temperature if inverse compton scattering wouldn't make a contribution
# and checking if this is indeed the case if so the coefficent and the temperature is returned. If inverse thomspn scattering makes a non-neglible contribution
# the numercial solver scipy fsolve is used to calculate the temperature which then is used to calculate the thompson coefficent. Based on eq. 13.
#This function also returns the calculated blackbody temperature.

def thompson_coefficent(luminosity_obs, stellar_radius,coupling_shell,density_breakout,time_breakout, time):
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
            temp_new = coupling_shell**2*temp_BB*np.power(0.5*np.log(3*np.power((density_breakout*(time_breakout/time))/(10**(-6)),-1/2)*np.power(temp_breakout*boltzmann_constant/100))*(1.6+np.log(3*np.power((density_breakout*(time_breakout/time))/(10**(-6)),-1/2)*np.power(temp_breakout*boltzmann_constant/100))),-2)
            if ((temp_obs_new-temp_obs)**2)**(0.5)<0.001:
                condition = 0
                temp_breakout = temp_new
                return thompson_coefficent, temp_breakout, temp_BB
            else:
                temp_breakout = temp_new
                thompson_coefficent = np.power(0.5*np.log(3*np.power((density_breakout*(time_breakout/time))/(10**(-6)),-1/2)*np.power(temp_breakout*boltzmann_constant/100))*(1.6+np.log(3*np.power((density_breakout*(time_breakout/time))/(10**(-6)),-1/2)*np.power(temp_breakout*boltzmann_constant/100))),-2)
        return thompson_coefficent, temp_breakout, temp_BB

#-----------------------------------------------------------------------------------------------------------------------------

def temp_observable(time,time_breakout,time_transition,velocity_breakout,poly_index, mass_breakout, stellar_radius, velo_index,density_breakout):
    mass_s = mass_shell(mass_breakout,time,time_transition,poly_index,velo_index)
    boltzmann_constant = 8.62*10**(-5)
    coupling_breakout = 0.2 * np.power(velocity_breakout*10**-7, 15 / 4) * np.power(density_breakout*10**(6), -1 / 8)
    stefan_boltzmann_constant = 5.67*10**(-8)
    temp_BB = np.power((mass_breakout*velocity_breakout**2)/(time_breakout*4*np.pi*stellar_radius**2*stefan_boltzmann_constant),1/4)
    coupling_shell = coupling_coefficent(mass_s, mass_breakout, time, time_breakout, time_transition, poly_index, velocity_breakout, density_breakout)
    if time <= time_transition:
        if coupling_shell < 1:
            temp_obs = temp_BB*np.power(coupling_breakout,(2*(poly_index+1)/(17*poly_index+9)))*np.power(time/time_breakout,-(2*(9*poly_index+5))/(3*(17*poly_index+9)))
            return temp_obs

        if coupling_shell > 1:
            temp_obs_breakout = coupling_breakout**2*temp_BB
            print(temp_obs_breakout)
            condition = 1
            while condition == 1:
                A = 3*np.power((density_breakout*(time_breakout/time_breakout))/(10**(-6)),-1/2)*np.power(temp_obs_breakout*boltzmann_constant/100,9/4)
                if A< 1 or 0.5*np.log(A)*(1.6+np.log(A))<1:
                    thompson_coefficent_breakout = 1
                    temp_obs_breakout = coupling_breakout**2*temp_BB
                    condition = 0
                    print("first loop xi is smaller than 1"+ str(temp_obs_breakout))

                else:

                    temp_new = coupling_shell**2*temp_BB*np.power(0.5*np.log(3*np.power((density_breakout*(time_breakout/time))/(10**(-6)),-1/2)*np.power(temp_obs_breakout*boltzmann_constant/100,9/4))*(1.6+np.log(3*np.power((density_breakout*(time_breakout/time))/(10**(-6)),-1/2)*np.power(temp_obs_breakout*boltzmann_constant/100,9/4))),-2)
                #    print("first loop" + str(temp_new))
                    if ((temp_new-temp_obs_breakout)**2)**(0.5)<0.001:
                        condition = 0
                        temp_obs_breakout = temp_new
                        thompson_coefficent_breakout = 0.5*np.log(A)*(1.6+np.log(A))


                    else:
                        temp_obs_breakout = temp_new

            temp_obs = temp_obs_breakout * np.power(time/time_breakout,-2/3)
            condition1 = 1
            while condition1 == 1:
                A = 3*np.power((density_breakout*(time_breakout/time_breakout))/(10**(-6)),-1/2)*np.power(temp_obs*boltzmann_constant/100,9/4)
                if A< 1 or 0.5*np.log(A)*(1.6+np.log(A))<1:
                    thompson_coefficent_obs = 1
                    temp_obs_new = temp_obs_breakout*np.power(time/time_breakout,-2/3)*np.power(thompson_coefficent_obs/thompson_coefficent_breakout,-2)
                    print(temp_obs_new)
                    print("second loop xi is smaller than 1")
                    temp_obs = temp_obs_new
                    condition1 = 0
                    return temp_obs

                else:
                    temp_obs_new = temp_obs_breakout*np.power(time/time_breakout,-2/3)*np.power((0.5*np.log(3*np.power((density_breakout*(time_breakout/time))/(10**(-6)),-1/2)*np.power(temp_obs*boltzmann_constant/100,9/4))*(1.6+np.log(3*np.power((density_breakout*(time_breakout/time))/(10**(-6)),-1/2)*np.power(temp_obs*boltzmann_constant/100,9/4))))/thompson_coefficent_breakout,-2)
                    print("second loop" + str(temp_obs_new))
                    if ((temp_obs_new-temp_obs)**2)**(0.5)<0.001:
                        condition1 = 0
                        temp_obs = temp_obs_new
                        return temp_obs
        #            thompson_coefficent_obs = np.power(0.5*np.log(3*np.power(7*density_breakout/(10**(-6)),-1/2)*np.power(temp_breakout*boltzmann_constant/100))*(1.6+np.log(3*np.power(7*density_breakout/(10**(-6)),-1/2)*np.power(temp_breakout*boltzmann_constant/100))),-2)
                    else:
                        temp_obs = temp_obs_new

    if time>time_transition:


        if coupling_shell < 1:

            temp_obs=(temp_BB
                    *np.power(coupling_breakout,2*(1.76*poly_index+1)/(22.32*poly_index+17))
                    *np.power(time_transition/time_breakout,-(8.03*poly_index+6)/(22.32*poly_index+17))
                    *np.power(time/time_transition,-(18.48*poly_index**2+20.69*poly_index+6)/((1.19*poly_index+1)*(22.32*poly_index+17))))

            return temp_obs

        if coupling_shell > 1:

            time_1 = (time_transition
                    *np.power(coupling_breakout*np.power(time_breakout/time_transition,1/6),3*(1.19*poly_index+1)/(9.88*poly_index+5)))
            time2 = (time_transition
                    *np.power(coupling_breakout*np.power(time_breakout/time_transition,1/6),6*(1.19*poly_index+1)/(12.48*poly_index+1)))


            if time<time2 and time>time_1:
                temp_obs = (temp_BB*np.power(coupling_breakout,2)*np.power(time_transition/time_breakout,-2/3)*np.power(time_1/time_transition,-(21.27*poly_index+11)/(3*(1.19*poly_index+1)))*np.power(time/time_1,-(3*poly_index+2)/(6*(1.19*poly_index+1))))

                return temp_obs

            if time > time2:

                temp_obs=(temp_BB
                    *np.power(coupling_breakout,2*(1.76*poly_index+1)/(22.32*poly_index+17))
                    *np.power(time_transition/time_breakout,-(8.03*poly_index+6)/(22.32*poly_index+17))
                    *np.power(time/time_transition,-(18.48*poly_index**2+20.69*poly_index+6)/((1.19*poly_index+1)*(22.32*poly_index+17))))

                return temp_obs
