import numpy as np


def mass_shell(mass_breakout,time,time_transition,poly_index,velo_index):
    mass_shell = (np.power(time/time_transition,(2*(poly_index+1))
    /((1+velo_index)*poly_index+1))*mass_breakout)
    return mass_shell


def luminosity(time,time_breakout,time_transition,mass_breakout,velocity_breakout,poly_index):

    luminosity_obs = ((mass_breakout*velocity_breakout**2)/(time_breakout)
                    *np.power(time_transition/time_breakout,-4/3)
                    *np.power(time/time_transition,-(2.28*poly_index-2)/(3*(1.19*poly_index+1))))

    return luminosity_obs


def temp_blackbody(stellar_radius,time,time_breakout,time_transition,mass_breakout,velocity_breakout,poly_index):
    stefan_boltzmann_constant = 5.67*10**(-8)
    temp_BB = np.power(luminosity(time,time_breakout,time_transition,mass_breakout,velocity_breakout,poly_index)/(4*np.pi*stellar_radius**2*stefan_boltzmann_constant),1/4)
    return temp_BB


def coupling_coefficent(mass_shell, mass_breakout, time, time_breakout, time_transition, poly_index, velocity_breakout, density_breakout):

    coupling_breakout = 0.2 * np.power(velocity_breakout*10**-7, 15 / 4) * np.power(density_breakout*10**(6), -1 / 8)
    coupling_shell= (coupling_breakout
                    *np.power(mass_shell / mass_breakout, -(22.32 * poly_index + 17) / (8 * poly_index + 8))
                    * np.power(time_transition / time_breakout, -1 / 6)
                    * np.power(time / time_transition, (42 * poly_index + 49) / (12 * (1.19 * poly_index + 1))))
    return coupling_shell


def temp_observable(time,time_breakout,time_transition,velocity_breakout,poly_index, mass_breakout, stellar_radius,velo_index,density_breakout):
    mass_s = mass_shell(mass_breakout,time,time_transition,poly_index,velo_index)
    coupling_shell = coupling_coefficent(mass_s, mass_breakout, time, time_breakout, time_transition, poly_index, velocity_breakout, density_breakout)
    coupling_breakout = 0.2 * np.power(velocity_breakout*10**-7, 15 / 4) * np.power(density_breakout*10**(6), -1 / 8)
    temp_BB = (mass_breakout*velocity_breakout**2)/(time_breakout)


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
