import numpy as np


def breakout_values(stellar_radius,stellar_mass,poly_index,opacity,velo_index):
    c = 3*10**8
    density_avg = Mass/R**3
    radius_breakout = stellar_radius*(1-np.power((opacity*stellar_radius*density_avg)/(c*(poly_index+1)),(1)/(poly_index*(velo_index-1)-1))/1800000)
    density_breakout = density_avg*((R-radius_breakout)/R)**(poly_index)
    mass_breakout = (4*np.pi*stellar_radius**3*density_avg)/(poly_index+1)*np.power(density_breakout/density_avg,(poly_index+1)/poly_index)
    velocity_breakout = 1800*10**3*(E_51/M_15)**(1/2)*np.power(density_breakout/density_avg, -velo_index)
    time_breakout = (stellar_radius-radius_breakout)/velocity_breakout
    time_transition = stellar_radius / velocity_breakout
    return radius_breakout, density_breakout, mass_breakout, time_breakout, velocity_breakout


def mass_shell(mass_breakout,time, time_transition,poly_index):
    if time>time_transition:
        mass_shell = np.power(time/time_transition,(2*(poly_index+1))/(1.19*poly_index+1))*mass_breakout
        return mass_shell
    else:
        mass_shell = mass_breakout
        return mass_shell

def luminosity(time,time_breakout,time_transition,mass_breakout,velocity_breakout,poly_index):
    if time>time_transition:
        luminosity_obs = (mass_breakout*velocity_breakout**2)/(time_breakout)*np.power(time_transition/time_breakout,-4/3)*np.power(time/time_transition,-(2.28*poly_index-2)/(3*(1.19*poly_index+1)))
        return luminosity_obs
    else:
        luminosity_obs = (mass_breakout*velocity_breakout**2)/(time_breakout)*np.power(time/time_breakout,-4/3)
        return luminosity_obs

def thermal_coupling(mass_shell, mass_breakout, time, time_breakout, time_transition, poly_index, velocity_breakout, density_breakout):
    eta_breakout = 0, 2 * np.power(velocity_breakout / 10 ** 4, 15 / 4) * np.power(density_breakout / 10 ** (-9), -1 / 8)
    eta = 0
    if time<time_transition:
        eta = eta_breakout * np.power(mass_shell / mass_breakout, -(17 * poly_index + 9) / (8 * poly_index + 8)) * np.power(time / time_breakout, -1 / 6)
    if time>time_transition:
        eta = eta_breakout * np.power(mass_shell / mass_breakout, -(22.32 * poly_index + 17) / (8 * poly_index + 8)) * np.power(time_transition / time_breakout, -1 / 6) * np.power(time / time_transition, (42 * poly_index + 49) / (12 * (1, 19 * poly_index + 1)))
    return eta


def temp_obervble(temp_BB,thermal_coupling,inverse_thomson,time,time_breakout,time_transition):


    return temp_obs
