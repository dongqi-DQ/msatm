from msatm.load_constants import load_constants
from msatm.e_sat import e_sat
from msatm.calculate_frac_ice import calculate_frac_ice
import numpy as np
import pdb

def saturation_adjustment(p_in,T_in,r_t_in,c=load_constants('default')):
    """
     Function to calculate the breakdown of water into vapor, liquid and solid 
     based on different microphysical assumptions
    
      [r,rl,ri,[rs]] = saturation_adjustment(p,T,r_t [type,ice,deltaT])
    
     Outputs: 
     r = water vapor mixing ratio (kg/kg)
     rl = liquid water mixing ratio (kg/kg)
     ri = solid water mixing ratio (kg/kg)
     rs = saturation mixing ratio (kg/kg)
    
     Inputs:
     p = pressure (Pa)
     T = temperature (K)
     r_t = total water mixing ratio (kg/kg)
     deltaT = mixed phase range (K) (default = 40)
    """
    # stupid copy for python reasons
    p,T,r_t = p_in.copy(),T_in.copy(),r_t_in.copy()


    
    es,junk = e_sat(T,c=c)
    rs = c.eps*(es/(p-es))
    r = np.minimum(r_t,rs)
    r[r<0] = 0
    fliq,fice,junk = calculate_frac_ice(T,c=c)
    

    rl = fliq*(r_t-r)
    ri = fice*(r_t-r)

    rl[np.abs(r_t-r)<1e-8] = 0
    ri[np.abs(r_t-r)<1e-8] = 0
    
    #if nargout > 3 varargout{1} = rs end
    varargout = rs

    return r,rl,ri,varargout 
