import numpy as np
import matplotlib.pyplot as plt

solar_M = 2e30
solar_R = 6.96e8
R = 500*solar_R
M = 15*solar_M
n = 1.5
L = 3.8e26
kappa = 0.034
mu = 0.19
E_51 = 1
M_15 = 1
c = 3e8
n_d = 10


def breakout_values(stellar_radius,stellar_mass,poly_index,opacity,velo_index,E_51,M_15):

    density_avg = stellar_mass/(stellar_radius**3)

    density_breakout = density_avg * np.power((c*(poly_index+1)/(opacity*stellar_radius*density_avg*1800e3))*np.power(M_15/E_51,0.5),poly_index/(poly_index*(1-velo_index)+1))

    radius_breakout = stellar_radius * (1 - np.power(density_breakout/density_avg,1/poly_index))

    mass_breakout = ((4*np.pi*stellar_radius**3*density_avg)/(poly_index+1))*np.power(density_breakout/density_avg,(poly_index+1)/poly_index)

    velocity_breakout = 1800e3*(E_51/M_15)**(1/2)*np.power(density_breakout/density_avg, -velo_index)

    time_breakout = (stellar_radius - radius_breakout) / velocity_breakout

    time_transition = stellar_radius / velocity_breakout


    return radius_breakout, density_breakout, mass_breakout, time_breakout, velocity_breakout, time_transition, density_avg

def luminosity(time,time_breakout,time_transition,mass_breakout,velocity_breakout,poly_index):
    if time>time_transition:

        luminosity_obs = (mass_breakout*velocity_breakout**2)/(time_breakout)*np.power(time_transition/time_breakout,-4/3)*np.power(time/time_transition,-(2.28*poly_index-2)/(3*(1.19*poly_index+1)))

        return luminosity_obs

    else:

        luminosity_obs = (mass_breakout*velocity_breakout**2)/(time_breakout)*np.power(time/time_breakout,-4/3)

        return luminosity_obs

r_0,rho_0,m_0,t_0,v_0,t_s,rho_star = breakout_values(R,M,n,kappa,mu,E_51,M_15)

print(t_s)

time = int(t_0)
luminosity_obs = np.zeros(n_d*24*3600-int(t_0))
#Mag = np.zeros(24*3600-int(t_0))
print_time = np.zeros(n_d*24*3600-int(t_0))
for time in range(int(t_0),n_d*24*3600):
    luminosity_obs[time-int(t_0)] = luminosity(time,t_0,t_s,m_0,v_0,n)
#    Mag[time-int(t_0)] = 4.77-2.5*np.log10(luminosity_obs[time]/L)
    print_time[time-int(t_0)] = time


plt.loglog(print_time,luminosity_obs)
plt.savefig("Plots\Luminosity.png")
