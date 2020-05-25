import numpy as np


def mass_shell(mass_breakout,ime, time_transition,poly_index):
    if time<time_transition:
        mass_shell = mass_breakout
    if time>time_transition:
        mass_shell = np.power(time/time_transition,(2*(poly_index+1))/(1.19*poly_index+1))*mass_breakout

    return mass_shell


def thermal_coupling(mass_shell, mass_breakout, time, time_breakout, time_transition, poly_index, velocity_breakout, density_breakout):
    eta_breakout = 0, 2 * np.power(velocity_breakout / 10 ** 4, 15 / 4) * np.power(density_breakout / 10 ** (-9), -1 / 8)
    eta = 0
    if time<time_transition:
        eta = eta_breakout * np.power(mass_shell / mass_breakout, -(17 * poly_index + 9) / (8 * poly_index + 8)) * np.power(time / time_breakout, -1 / 6)
    if time>time_transition:
        eta = eta_breakout * np.power(mass_shell / mass_breakout, -(22, 32 * poly_index + 17) / (8 * poly_index + 8)) * np.power(time_transition / time_breakout, -1 / 6) * np.power(time / time_transition, (42 * poly_index + 49) / (12 * (1, 19 * poly_index + 1)))
    return eta


def temp_obervble(temp_BB,thermal_coupling,inverse_thomson,time,time_breakout,time_transition):


    return temp_obs
