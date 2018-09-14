from msatm.calculate_frac_ice import calculate_frac_ice
from msatm.load_constants import load_constants
import numpy as np
import pdb

def simple_fallout(T_in,cond_in,rl_in,ri_in,c=load_constants('default'),gamma=0):
    """
     A very simple precipitation parameterization for use in parcel ascent calculations.
     Assumes a fraction gamma of the condensation is reomved from the parcel
    
     inputs:
    		T	temperature (K)
    		cond	amount of condensation (kg condesate / kg dry air)
    		rl	liquid water mixing ratio (kg/kg)
    		ri	ice water mixing ratio (kg/kg)
    		gamma   fraction to precipitate (0-1)
    """
    # stupid copy for python reasons
    T,cond,rl,ri = T_in.copy(),cond_in.copy(),rl_in.copy(),ri_in.copy()

    # Break the condensation into ice and liquid
    fliq,fice,junk = calculate_frac_ice(T,c=c)
 
    # Remove a fraction gamma
    rl_new = np.maximum(0,rl - gamma*(cond)*fliq)
    ri_new = np.maximum(0,ri - gamma*(cond)*fice)

    return rl_new,ri_new
