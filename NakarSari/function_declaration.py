#This file is meant to store the declarations of all used functions


import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp



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

    energy_breakout = mass_breakout*np.power(velocity_breakout,2)

    coupling_breakout = 0.2 * np.power(velocity_breakout/(10**7), 15 / 4) * np.power(density_breakout/(10**(-6)), -1 / 8)

    opticaldepth_breakout = opacity*mass_breakout/stellar_radius**2


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

        luminosity_obs = ((energy_breakout)
        /(time_breakout)
        *np.power(time/time_breakout
        ,-4/3))

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

#-----------------------------------------------------------------------------------------------------------------------

def temp_blackbody(stellar_radius,time,time_breakout,time_transition,mass_breakout,velocity_breakout,poly_index,energy_breakout):
    stefan_boltzmann_constant = 5.67*10**(-8)
    temp_BB = np.power(luminosity(time,time_breakout,time_transition,energy_breakout,poly_index)/(4*np.pi*stellar_radius**2*stefan_boltzmann_constant),1/4)
    return temp_BB
#There is a 4  and a pi missing in front of pi

#----------------------------------------------------------------------------------------------------------

#the coefficent for inverse thomson scattering (referred to as xi) is calculated by first finding the relevant blackbody temperature based
# on the stefan-boltzmann law and the luminosity then calculating the temperature if inverse compton scattering wouldn't make a contribution
# and checking if this is indeed the case if so the coefficent and the temperature is returned. If inverse thomspn scattering makes a non-neglible contribution
# the numercial solver scipy fsolve is used to calculate the temperature which then is used to calculate the thompson coefficent. Based on eq. 13.
#This function also returns the calculated blackbody temperature.

def thompson_coefficent(temp_obs,density_breakout,time_breakout,time):
    #print('i have been summoned')
    boltzmann_constant = 8.617*10**(-5)
    density_unit = 7*density_breakout*10**(-3)*(time_breakout/time)
    temp_unit = temp_obs*boltzmann_constant
    yMAX = 3*np.power(density_unit/10**(-9),-1/2)*np.power(temp_unit/100,9/4)
    thompson = (1/2)*np.log(yMAX)*(1.6+np.log(yMAX))

    if thompson < 1:
        thompson = 1
        return thompson

    else:
        return  thompson




def temp_break(temp_obs_breakout,density_breakout,time_breakout,coupling_breakout,temp_BB_breakout):

    return temp_obs_breakout - (coupling_breakout**2/thompson_coefficent(temp_obs_breakout,density_breakout,time_breakout,time_breakout)**2)*temp_BB_breakout

def temp_obser(temp_obs,density_breakout,time_breakout,time,coupling_shell,temp_BB_breakout,temp_obs_breakout,thompson_coefficent_breakout):

    return temp_obs -  temp_obs_breakout*np.power(time/time_breakout,-2/3)*np.power(thompson_coefficent(temp_obs,density_breakout,time_breakout,time)/thompson_coefficent_breakout,-2)





#-----------------------------------------------------------------------------------------------------------------------------

def temp_observable(time,time_breakout,time_transition,velocity_breakout,poly_index, mass_breakout, stellar_radius, velo_index,density_breakout,energy_breakout,coupling_breakout, depth_breakout,radius_breakout,pre_temp):
    mass_s = mass_shell(mass_breakout,time,time_transition,poly_index,velo_index)
#    print(mass_s)
    boltzmann_constant = 8.62*10**(-5)


    stefan_boltzmann_constant = 5.67*10**(-8)
    temp_BB_breakout = temp_blackbody(stellar_radius,time_breakout,time_breakout,time_transition,mass_breakout,velocity_breakout,poly_index,energy_breakout)
    temp_break1 = (np.power((luminosity(time_breakout,time_breakout,time_transition,energy_breakout,poly_index)*depth_breakout)
           /(radius_breakout**2*stefan_boltzmann_constant*c),
           1/4))
    #print('easy one = ' + str(temp_BB_breakout))
    #print('Nakar and Saris conoctation = ' + str(temp_break1))

    coupling_shell = coupling_coefficent(mass_s, mass_breakout, time, time_breakout, time_transition, poly_index, velocity_breakout, density_breakout,coupling_breakout)

    if time < time_transition:

        if coupling_shell < 1:#-------------------------------------------------------------------------------------------------------------------------------------------------------

            temp_obs = temp_BB_breakout*np.power(coupling_breakout,(2*(poly_index+1)/(17*poly_index+9)))*np.power(time/time_breakout,-(2*(9*poly_index+5))/(3*(17*poly_index+9)))
            return temp_obs

        if coupling_shell > 1:#---------------------------------------------------------------------------------------------------------------------------------------------------------

            temp_obs_breakout = coupling_breakout**2*temp_BB_breakout
            temp_array = np.array([temp_obs_breakout])
        #    print(temp_obs_breakout)
            itcount_breakout = 0
            terminator_breakout = 1

            while terminator_breakout == 1:

                if thompson_coefficent(temp_obs_breakout,density_breakout,time_breakout,time_breakout) == 1:

                    terminator_breakout = 0
                    thompson_coefficent_breakout = 1
                #    print("breakout thompson coefficent = 1 | " + str(itcount_breakout)+ " | "+str(temp_obs_breakout))


                else:

                    temp_iter = temp_BB_breakout*coupling_breakout**2/thompson_coefficent(temp_obs_breakout,density_breakout,time_breakout,time_breakout)**(2)

                    if np.abs(temp_iter-temp_obs_breakout) < 0.0001:

                        temp_obs_breakout = temp_iter
                        thompson_coefficent_breakout = thompson_coefficent(temp_obs_breakout,density_breakout,time_breakout,time_breakout)
                        terminator_breakout = 0
                    #    print("breakout solved by iteration | " + str(itcount_breakout)+ " | " +str(temp_obs_breakout))


                    else:
                    #print(thompson_coefficent(temp_obs_breakout,density_breakout,time_breakout,time_breakout))
                        temp_obs_breakout = temp_iter
                        temp_array = np.append(temp_array,temp_obs_breakout)

                        itcount_breakout += 1

            #x_shu = np.linspace(0,len(temp_array),len(temp_array))
            #plt.scatter(x_shu, temp_array)
            #plt.yscale('log')
            #plt.show()
            temp_obs_breakout1 = sp.fsolve(temp_break,coupling_breakout**2*temp_BB_breakout,args=(density_breakout,time_breakout,coupling_breakout,temp_BB_breakout))
            temp_obs_breakout2 = sp.root(temp_break,coupling_breakout**2*temp_BB_breakout,args=(density_breakout,time_breakout,coupling_breakout,temp_BB_breakout),method = 'anderson')
            #print('this is the fsolve one'+str(temp_obs_breakout1))
        #    print('this is the root one'+str(temp_obs_breakout2.x))
        #    print('this is mine'+str(temp_obs_breakout))

#----------------------------------------------------------------------------------------------------------------------------------------------
            temp_obs_breakout = temp_obs_breakout2.x
            thompson_coefficent_breakout = thompson_coefficent(temp_obs_breakout,density_breakout,time_breakout,time_breakout)

            if time == int(time_breakout):
                print('i have been there')
                return temp_obs_breakout


            else:
                temp_obs = pre_temp
                terminator = 1
                itcount = 0
                temp = sp.fsolve(temp_obser, temp_obs,args=(density_breakout,time_breakout,time,coupling_shell,temp_BB_breakout,temp_obs_breakout,thompson_coefficent_breakout))

                return temp




                while terminator == 1:

                    if thompson_coefficent(temp_obs,density_breakout,time_breakout,time) == thompson_coefficent(temp_obs,density_breakout,time_breakout,time):

                        terminator = 0
                        temp_obs = pre_temp
                        #    print("observable thompson coefficent = 1 | " + str(itcount)+ " | "+str(temp_obs))
                        return temp_obs

                    else:

                        temp_obs_new = temp_obs_breakout*np.power(time/time_breakout,-2/3)*np.power(thompson_coefficent(temp_obs,density_breakout,time_breakout,time)/thompson_coefficent_breakout,-2)



                        if np.abs(temp_obs_new-temp_obs) < 0.1:

                            terminator = 0
                            temp_obs = temp_obs_new
                            #        print("observable solved by iteration | " + str(itcount)+ " | " +str(temp_obs))
                            #x_shu = np.linspace(0,len(temp_array),len(temp_array))
                            #plt.scatter(x_shu, temp_array)
                            #plt.yscale('log')
                            #plt.show()
                            return temp_obs

                        else:
                            if itcount > 500:
                                terminater = 0
                                return temp_obs
                            temp_obs = temp_obs_new
                            #print('solving observable' +str(itcount))
                            temp_array = np.append(temp_array,temp_obs)
                            itcount += 1



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
