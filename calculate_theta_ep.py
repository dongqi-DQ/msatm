from load_constants import load_constants
import numpy as np

def calculate_theta_ep(T_in,r_in,p_in,c=load_constants('bolton')):
    """
     Calculate pseudo-equivelant potential temperature
     assuming no ice
     Use Bolton formula, taken from Emanuel (1994), p 132
     This assumes no ice!
    """
    # stupid copy for python reasons
    T,r,p = T_in.copy(),r_in.copy(),p_in.copy()

    evap = r*p/(c.eps + r) 
    Tstar = 2840/(3.5*np.log(T) - np.log(evap/100) - 4.805 ) + 55 
    theta_ep = T*(100000.0/p)**(0.2854*(1-0.28*r))*np.exp( r * (1+0.81*r) * (3376.0/Tstar -2.54) )

    return theta_ep,Tstar
