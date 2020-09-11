import numpy as np
import matplotlib.pyplot as plt

#---------------------------------------------------------------------------------------
c = 3e8
a = 7.56*10**(-16)
mass_solar = 2e30
radius_solar = 6.96e8

R = 500*radius_solar

M = 20 *mass_solar
M_15 = M/(15*mass_solar)

E_51 = 5

n = 1.5
kappa = 0.034

beta = 0.19

#---------------------------------------------------------------------------------------

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

def luminosity(t,t_0,t_S,E_0,n):
    if t>t_S:

        luminosity_obs = ((E_0)
        /(t_0)
        *np.power(t_S/t_0
        ,-4/3)
        *np.power(t/t_S
        ,-(2.28*n-2)
        /(3*(1.19*n+1))))

        return luminosity_obs

    else:

        luminosity_obs = ((E_0)
        /(t_0)
        *np.power(t/t_0
        ,-4/3))

        return luminosity_obs

def mass_shell(m_0,r_0,n,t,v_0,r):
    d = (R-r) + v_0*t
    d_0 = (R-r_0)
    m = m_0*(d/d_0)**(n+1)
    return m

def optical_depth(m_0,tau_0,r_0,r,n,t,v_0):
    m = mass_shell(m_0,r_0,n,t,v_0,r)
    tau = (m/m_0)*(r_0**2/r**2)
    return tau

def temp_blackbody(r,t,t_0,t_S,m_0,v_0,n,E_0,r_0,tau_0):
    m = mass_shell(m_0,r_0,n,t,v_0,r)
    T_BB = ((luminosity(t,t_0,t_S,E_0,n)*optical_depth(m_0,tau_0,r_0,r,n,t,v_0))/(c*a*r**2))**(1/4)
    return T_BB

def coupling_coefficent(eta_0,r,t,t_0,t_S,m_0,v_0,n,E_0,r_0,tau_0):
    m = mass_shell(m_0,r_0,n,t,v_0,r)
    d_0 = (R-r_0)
    rho_i = rho_0*(m/m_0)**(n/(n+1))
    rho = rho_i * (d_0/(v_0*t))
    T_BB0 = temp_blackbody(r,t_0,t_0,t_S,m_0,v_0,n,E_0,r_0,tau_0)
    T_BB = temp_blackbody(r,t,t_0,t_S,m_0,v_0,n,E_0,r_0,tau_0)
    eta = eta_0*(T_BB/T_BB0)**(7/2)*(rho/rho_0)**(-2)*(t_0/t)
    return eta

def coupling_coefficent1(eta_0,r,t,t_0,t_S,m_0,v_0,n,E_0,r_0,tau_0):
    m = mass_shell(m_0,r_0,n,t,v_0,r)
    eta = eta_0*(m/m_0)**(-(17*n+9)/(8*(n+1)))*(t/t_0)**(-1/6)
    return eta

def coupling_coefficent2(eta_0,r,t,t_0,t_S,m_0,v_0,n,E_0,r_0,tau_0):
    m = mass_shell(m_0,r_0,n,t,v_0,r)
    eta = eta_0*(m/m_0)**(-(16*beta*n+17*n+25)/(8*(n+1)))*(t/t_0)**(-1/6)
    return eta
#-------------------------------------------------------------------------------------------------
r_0,rho_0,m_0,t_0,v_0,t_S,rho_star,E_0,eta_0,tau_0 = breakout_values(R,M,n,kappa,beta,E_51,M_15)

set_r = np.linspace(r_0,R,1000)
T = np.zeros(len(set_r))
k = np.zeros(len(set_r))
m = np.zeros(len(set_r))
eta1 = np.zeros(len(set_r))
eta2 = np.zeros(len(set_r))
#fig, ax1 = plt.subplots()
#ax2 = ax1.twinx()
t = t_0 +20
print(t)

for i in range(0,len(set_r)):
    m[i] = mass_shell(m_0,r_0,n,t,v_0,set_r[i])
    T[i] = temp_blackbody(set_r[i],t,t_0,t_S,m_0,v_0,n,E_0,r_0,tau_0)
    eta1[i] = coupling_coefficent1(eta_0,set_r[i],t_0,t_0,t_S,m_0,v_0,n,E_0,r_0,tau_0)
    eta2[i] = coupling_coefficent2(eta_0,set_r[i],t_0,t_0,t_S,m_0,v_0,n,E_0,r_0,tau_0)
#    k[i] = eta[i]**2*T[i]
#ax1.plot(set_r,T, label = '$T_{BB}$',color = 'red')
plt.plot(m,eta1,label=' $\eta, Approximated Density$', color = 'red')
plt.plot(m,eta2,label=' $\eta$, Derived Density', color = 'green')
#ax1.plot(set_r,k,label='$T_{BB}\eta^2$' )


#lines, labels = ax1.get_legend_handles_labels()
#lines2, labels2 = ax2.get_legend_handles_labels()
#ax1.set_ylabel('Temperature [K] / $T_{BB}\eta^2$')
#ax2.set_ylabel('$\eta$')
#ax1.set_xlabel('Radius R [m]')
#ax2.legend(lines + lines2, labels + labels2, loc = 'center left')
plt.ylabel('$\eta$')
plt.xlabel('Mass [kg]')
#plt.yscale('log')
plt.legend()
plt.grid()
plt.show()
