import numpy as np
from scipy.optimize import fsolve


def breakout_values(stellar_radius,stellar_mass,poly_index,opacity,velo_index,E_51,M_15):
    c = 3*10**8
    density_avg = stellar_mass/stellar_radius**3
    radius_breakout = stellar_radius*(1-np.power(1800*10**3*np.power(E_51/M_15,0.5)*(opacity*stellar_radius*density_avg)/(c*(poly_index+1)),(1)/(poly_index*(velo_index-1)-1))/1800000)
    density_breakout = density_avg*np.power((stellar_radius-radius_breakout)/stellar_radius,poly_index)
    mass_breakout = (4*np.pi*stellar_radius**3*density_avg)/(poly_index+1)*np.power(density_breakout/density_avg,(poly_index+1)/poly_index)
    velocity_breakout = 1800*10**3*(E_51/M_15)**(1/2)*np.power(density_breakout/density_avg, -velo_index)
    time_breakout = (stellar_radius-radius_breakout)/velocity_breakout
    time_transition = stellar_radius / velocity_breakout
    return radius_breakout, density_breakout, mass_breakout, time_breakout, velocity_breakout


def mass_shell(mass_breakout,time,time_transition,poly_index):
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

def coupling_coefficent(mass_shell, mass_breakout, time, time_breakout, time_transition, poly_index, velocity_breakout, density_breakout):
    coupling_breakout = 0, 2 * np.power(velocity_breakout / 10 ** 4, 15 / 4) * np.power(density_breakout / 10 ** (-9), -1 / 8)
    coupling_shell = 0

    if time<time_transition:
        coupling_shell = eta_breakout * np.power(mass_shell / mass_breakout, -(17 * poly_index + 9) / (8 * poly_index + 8)) * np.power(time / time_breakout, -1 / 6)
    if time>time_transition:
        coupling_shell = eta_breakout * np.power(mass_shell / mass_breakout, -(22.32 * poly_index + 17) / (8 * poly_index + 8)) * np.power(time_transition / time_breakout, -1 / 6) * np.power(time / time_transition, (42 * poly_index + 49) / (12 * (1, 19 * poly_index + 1)))
    return coupling_shell


def thompson_coefficent(luminosity_obs, stellar_radius,coupling_shell,density_shell):
    boltzmann_constant = 5.67*10**(-8)
    k_b = 8.62*10**(-5)
    temp_BB = np.power(luminosity_obs/(4*np.pi*stellar_radius**2*boltzmann_constant),1/4)
    temp_breakout = coupling_shell**2*temp_BB
    yMAX = 3*np.power(density_shell/(10**(-6)),-1/2)*np.power(temp_breakout*k_b/100)

    if yMAX<1:
        thompson_coefficent = 1
        return thompson_coefficent, temp_breakout, temp_BB
    if 1>0.5*np.log(ymax)*(1.6+np.log(ymax)):
        thompson_coefficent = 1
        return thompson_coefficent, temp_breakout, temp_BB
    else:
        eq = lambda temp : coupling_shell**2*temp_BB*np.power(0.5*np.log(3*np.power(density_shell/(10**(-6)),-1/2)*np.power(temp*k_b/100))*(1.6+np.log(3*np.power(density_shell/(10**(-6)),-1/2)*np.power(temp*k_b/100))),-2) - temp
        temp_breakout = fsolve(eq,temp_breakout)
        thompson_coefficent = np.power(0.5*np.log(3*np.power(density_shell/(10**(-6)),-1/2)*np.power(temp_breakout*k_b/100))*(1.6+np.log(3*np.power(density_shell/(10**(-6)),-1/2)*np.power(temp_breakout*k_b/100))),-2)
        return thompson_coefficent, temp_breakout, temp_BB



def temp_observable(temp_BB,coupling_shell,thompson_shell,time,time_breakout,time_transition):


    return temp_obs
