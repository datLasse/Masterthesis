import numpy as np
poly_index = 1.5
mu = 0.19
R = 800*6.9*10**(8)
Mass = 15*1.9*10**(30)
kappa = 0.034

def breakout_point(R,Mass,poly_index,kappa,mu):
    condition = 1
    r = 800*6*10**(8)
    while condition==1:
        c = 3*10**8
        RHO = Mass/R**3
        rho = RHO*((R-r)/R)**(poly_index)
        v = (1800*10**3)*(rho/RHO)**(-mu)
        tau = (kappa*R*RHO)/(poly_index+1)*np.power(rho/RHO,(poly_index+1)/poly_index)
        if tau == c/v:
            print("here we are boys",r)
            condition = 0
            return r
        if r>R:
            print("it happend you donkey")
            condition = 0
            return 0
        else:
            r=r+1
    return 0

r = breakout_point(R,Mass,poly_index,kappa,mu)
