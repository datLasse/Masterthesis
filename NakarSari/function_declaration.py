#This file is meant to store the declarations of all used functions


import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp



c = 3e8
a = 7.5657*10**(-16)
boltzmann_constant = 8.62*10**(-5)
stefan_boltzmann_constant = 5.67*10**(-8)


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

    energy_breakout = mass_breakout*np.power(velocity_breakout,2)

    coupling_breakout = 0.2 * np.power(velocity_breakout/(10**7), 15 / 4) * np.power(density_breakout/(10**(-6)), -1 / 8)

    opticaldepth_breakout = opacity*stellar_radius*density_avg/(poly_index+1)*np.power(density_breakout/density_avg,(poly_index+1)/poly_index)


    return radius_breakout, density_breakout, mass_breakout, time_breakout, velocity_breakout, time_transition, density_avg, energy_breakout,coupling_breakout, opticaldepth_breakout



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

def luminosity(time,time_breakout,time_transition,energy_breakout,poly_index):
    if time>time_transition:

        luminosity_obs = ((energy_breakout)
        /(time_breakout)
        *np.power(time_transition/time_breakout
        ,-4/3)
        *np.power(time/time_transition
        ,-(2.28*poly_index-2)
        /(3*(1.19*poly_index+1))))

        return luminosity_obs

    else:

        luminosity_obs = ((energy_breakout)/(time_breakout))*np.power(time/time_breakout,-4/3)

        return luminosity_obs

#------------------------------------------------------------------------------------------------------

#the thermal coupling coefficent referred to as eta in the paper is calculated based on eq.10,16,18

def coupling_coefficent(mass_shell, mass_breakout, time, time_breakout, time_transition, poly_index, velocity_breakout, density_breakout, coupling_breakout):
    coupling_shell = 0

    if time == time_breakout:
        coupling_shell = coupling_breakout
    if time<time_transition:
        coupling_shell = (coupling_breakout * np.power(mass_shell / mass_breakout, -(17 * poly_index + 9)/(8 * poly_index + 8)) *
         np.power(time / time_breakout, -1/6))
    if time>time_transition:
        coupling_shell = coupling_breakout * np.power(mass_shell / mass_breakout, -(22.32 * poly_index + 17)/(8 * poly_index + 8)) * np.power(time_transition / time_breakout, -1 / 6) * np.power(time / time_transition, (42 * poly_index + 49) / (12 * (1.19 * poly_index + 1)))
    return coupling_shell

#------------------------------------------------------------------------------------------------------------------------------------

def radius(stellar_radius,time,time_transition,poly_index,velo_index):
    if time> time_transition:
        r = stellar_radius*np.power(time/time_transition,((1-velo_index)*poly_index+1)/((1+velo_index)*poly_index+1))
        return r

    else:
        return stellar_radius


#-----------------------------------------------------------------------------------------------------------------------

def optical_depth(mass_breakout,opticaldepth_breakout,radius_breakout,radius,poly_index,time,velocity_breakout,velo_index,time_transition):
    m  = mass_shell(mass_breakout,time,time_transition,poly_index,velo_index)
    if time>time_transition:
        optical_depth = opticaldepth_breakout*(m/mass_breakout)**(((1+2*velo_index)*poly_index+1)/(poly_index+1))
        return optical_depth

    else:
        optical_depth = opticaldepth_breakout
        return optical_depth

#-----------------------------------------------------------------------------------------------------------------------

def temp_blackbody(time,time_breakout,time_transition,mass_breakout,velocity_breakout,poly_index,energy_breakout,radius_breakout,opticaldepth_breakout,velo_index,stellar_radius):
    m = mass_shell(mass_breakout,time,time_transition,poly_index,velo_index)
    T_BB = ((luminosity(time,time_breakout,time_transition,energy_breakout,poly_index)*optical_depth(mass_breakout,opticaldepth_breakout,radius_breakout,radius,poly_index,time,velocity_breakout,velo_index,time_transition))/(c*a*np.power(radius(stellar_radius,time,time_transition,poly_index,velo_index),2)))**(1/4)
    return T_BB

#----------------------------------------------------------------------------------------------------------

#the coefficent for inverse thomson scattering (referred to as xi) is calculated by first finding the relevant blackbody temperature based
# on the stefan-boltzmann law and the luminosity then calculating the temperature if inverse compton scattering wouldn't make a contribution
# and checking if this is indeed the case if so the coefficent and the temperature is returned. If inverse thomspn scattering makes a non-neglible contribution
# the numercial solver scipy fsolve is used to calculate the temperature which then is used to calculate the thompson coefficent. Based on eq. 13.
#This function also returns the calculated blackbody temperature.

def thompson_coefficent(temp_obs,density_breakout,time_breakout,time):
    boltzmann_constant = 8.617*10**(-5)
    density_unit = 7*density_breakout*10**(-3)*(time_breakout/time)
    temp_unit = temp_obs*boltzmann_constant
    yMAX = 3*np.power(7*density_unit/10**(-9),-1/2)*np.power(temp_unit/100,9/4)
    thompson = 0.5*np.log(yMAX)*(1.6 + np.log(yMAX))

    if thompson < 1:
        thompson = 1
        return thompson

    else:
        return  thompson







#-----------------------------------------------------------------------------------------------------------------------------

def temp_observable(time,time_breakout,time_transition,velocity_breakout,poly_index, mass_breakout, stellar_radius, velo_index,density_breakout,energy_breakout,coupling_breakout, depth_breakout,radius_breakout,pre_temp,acc):
    mass_s = mass_shell(mass_breakout,time,time_transition,poly_index,velo_index)
#    print(mass_s)

    temp_BB_breakout = temp_blackbody(time_breakout,time_breakout,time_transition,mass_breakout,velocity_breakout,poly_index,energy_breakout,radius_breakout,depth_breakout,velo_index,stellar_radius)


    coupling_shell = coupling_coefficent(mass_s, mass_breakout, time, time_breakout, time_transition, poly_index, velocity_breakout, density_breakout,coupling_breakout)#


    if time < time_transition:

        if coupling_shell < 1:#-------------------------------------------------------------------------------------------------------------------------------------------------------

            temp_obs = temp_BB_breakout*np.power(coupling_breakout,(2*(poly_index+1)/(17*poly_index+9)))*np.power(time/time_breakout,-(2*(9*poly_index+5))/(3*(17*poly_index+9)))

            return temp_obs

        if coupling_shell > 1:#---------------------------------------------------------------------------------------------------------------------------------------------------------

            compton_breakout = thompson_coefficent(temp_BB_breakout*coupling_breakout**2,density_breakout,time_breakout,time_breakout)

            if compton_breakout == 1:
                temp_obs_breakout = temp_BB_breakout*coupling_breakout**2

            else:

                T = np.linspace(temp_BB_breakout,temp_BB_breakout*coupling_breakout**2,1000)
                xi = np.zeros(len(T))
                Txi = np.zeros(len(T))
                comp = np.zeros(len(T))
                termination = 1

                while termination == 1:

                    for i in range(len(T)):

                        xi[i] = thompson_coefficent(T[i],density_breakout,time_breakout,time_breakout)
                        Txi[i] = temp_BB_breakout*coupling_breakout**2*(xi[i])**(-2)
                        comp[i] = Txi[i] - T[i]

                        if np.abs(comp[i]) < acc:
                            temp_obs_breakout = T[i]
                            #print('done' + str(T[i]))
                            termination = 0

                        if comp[i-1] >0 and comp[i]<0 :
                            p = i


                    #fig, ax1 = plt.subplots()
                    #ax2 = ax1.twinx()
                    #ax1.plot(T,comp, label = 'T_calculated - T_coordinate',color ='red')
                    #ax2.plot(T,xi, label = 'Compton scattering correction',color = 'green')
                    #ax1.plot(T,Txi, label = 'T_calculated', color = 'blue')
                    #lines, labels = ax1.get_legend_handles_labels()
                    #lines2, labels2 = ax2.get_legend_handles_labels()

                    #ax1.set_ylabel('Temperature T_calculated [K]')
                    #ax2.set_ylabel('Compton scattering correction')
                    #ax1.set_xlabel('Temperature T_coordinate [K]')
                    #ax1.axvspan(3.1e6,3.9e6,ymin = 0.43 , ymax= 0.55, color='black',fill = 0)
                    #ax2.legend(lines + lines2, labels + labels2, loc = 'upper center')
                    #ax1.grid()
                    #plt.show()

                    T_iter = np.linspace(T[p-1],T[p],1000)
                    T = T_iter

#----------------------------------------------------------------------------------------------------------------------------------------------
            thompson_coefficent_breakout = thompson_coefficent(temp_obs_breakout,density_breakout,time_breakout,time_breakout)

            if time == int(time_breakout):

                return temp_obs_breakout


            else:
                temp_lower = temp_blackbody(time,time_breakout,time_transition,mass_breakout,velocity_breakout,poly_index,energy_breakout,radius_breakout,depth_breakout,velo_index,stellar_radius)
                temp_upper = temp_lower*coupling_shell**2
                compton = thompson_coefficent(temp_upper,density_breakout,time_breakout,time)

                if compton == 1:
                    temp = temp_upper
                    return temp

                else:

                    T = np.linspace(temp_lower,temp_upper,1000)
                    xi = np.zeros(len(T))
                    Txi = np.zeros(len(T))
                    comp = np.zeros(len(T))
                    termination = 1

                    while termination == 1:

                        for i in range(len(T)):

                            xi[i] = thompson_coefficent(T[i],density_breakout,time_breakout,time)
                            Txi[i] = temp_obs_breakout*np.power(time/time_breakout,-2/3)*np.power(xi[i]/thompson_coefficent_breakout,-2)
                            comp[i] = Txi[i] - T[i]

                            if np.abs(comp[i]) < acc:
                                temp = T[i]
                                #print('done' + str(T[i]))
                                termination = 0
                                return temp

                            if comp[i-1] >0 and comp[i]<0 :
                                p = i


                        #fig, ax1 = plt.subplots()
                        #ax2 = ax1.twinx()
                        #ax1.plot(T,comp, label = 'T_calculated - T_coordinate second loop',color ='red')
                        #ax2.plot(T,xi, label = 'Compton scattering correction',color = 'green')
                        #ax1.plot(T,Txi, label = 'T_calculated', color = 'blue')
                        #lines, labels = ax1.get_legend_handles_labels()
                        #lines2, labels2 = ax2.get_legend_handles_labels()

                        #ax1.set_ylabel('Temperature T_calculated [K]')
                        #ax2.set_ylabel('Compton scattering correction')
                        #ax1.set_xlabel('Temperature T_coordinate [K]')
                    #ax1.axvspan(3.1e6,3.9e6,ymin = 0.43 , ymax= 0.55, color='black',fill = 0)
                        #ax2.legend(lines + lines2, labels + labels2, loc = 'upper center')
                        #ax1.grid()
                        #plt.show()



                        T_iter = np.linspace(T[p-1],T[p],1000)
                        T = T_iter




    if time>time_transition:#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------


        if coupling_shell < 1:#------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            temp_obs=(temp_BB_breakout
                    *np.power(coupling_breakout,2*(1.76*poly_index+1)/(22.32*poly_index+17))
                    *np.power(time_transition/time_breakout,-(8.03*poly_index+6)/(22.32*poly_index+17))
                    *np.power(time/time_transition,-(18.48*poly_index**2+20.69*poly_index+6)/((1.19*poly_index+1)*(22.32*poly_index+17))))


            return temp_obs

        if coupling_shell > 1:#------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            time_1 = (time_transition
                    *np.power(coupling_breakout*np.power(time_breakout/time_transition,1/6),3*(1.19*poly_index+1)/(9.88*poly_index+5)))
            time2 = (time_transition
                    *np.power(coupling_breakout*np.power(time_breakout/time_transition,1/6),6*(1.19*poly_index+1)/(12.48*poly_index+1)))

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            if time<time2 and time>time_1:

                temp_obs = (temp_BB_breakout*coupling_breakout**2*np.power(time_transition/time_breakout,-2/3)*np.power(time_1/time_transition,-(21.27*poly_index+11)/(3*(1.19*poly_index+1)))*np.power(time/time_1,-(3*poly_index+2)/(6*(1.19*poly_index+1))))

                return temp_obs

#--------------------------------------------------------------------------------------------------------------------------------------------

            if time > time2:

                temp_obs=(temp_BB_breakout
                    *np.power(coupling_breakout,2*(1.76*poly_index+1)/(22.32*poly_index+17))
                    *np.power(time_transition/time_breakout,-(8.03*poly_index+6)/(22.32*poly_index+17))
                    *np.power(time/time_transition,-(18.48*poly_index**2+20.69*poly_index+6)/((1.19*poly_index+1)*(22.32*poly_index+17))))

                return temp_obs
