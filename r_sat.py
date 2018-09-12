from load_constants import load_constants
import numpy as np
from e_sat import e_sat

def r_sat(T_in,p_in,c=load_constants('default')):
    """
     Function to calculate saturation mixing ratio
      
     Calling format:
    
     qs = q_sat(T,p,[type,ice,deltaT)
    
     qs = saturation mixing ratio (kg/kg)
     T = temperature (K)
     p = pressure (Pa)
    
     See e_sat(T) function for info about optional arguments
    """
    # copy objects
    T,p = T_in.copy(),p_in.copy()

    es,junk = e_sat(T,c=c);

    rs = c.eps * es/(p-es);

    return rs

