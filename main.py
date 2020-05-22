import numpy as np



def eta(m,m_0,t,t_0,t_s,eta_0,n):
    eta = 0
    if t<t_s:
        eta = eta_0*np.power(m/m_0,-(17*n+9)/(8*n+8))*np.power(t/t_0,-1/6)
    if t>t_s:
        eta = eta_0*np.power(m/m_0,-(22,32*n+17)/(8*n+8))*np.power(t_s/t_0,-1/6)*np.power(t/t_s,(42*n+49)/(12*(1,19*n+1)))
    return eta


def Tobs(T_BB,eta,xi,t,t_0,t_s):


    return Tobs
