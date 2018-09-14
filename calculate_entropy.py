from msatm.load_constants import load_constants
import numpy as np
from msatm.r_sat import r_sat
from msatm.saturation_adjustment import saturation_adjustment

def calculate_entropy(T_in,p_in,rt_in,c=load_constants('default'),**kwargs):
    """
     Function to calculate the entropy

     s = calculate_entropy(T,p,rt [,rl,ri,type,ice,deltaT])

     Calculates entropy based on pressure (p Pa) Temperature (T K) and total water mixing ratio (rt kg/kg))

     Optional arguments

     rl = liquid mixing ratio (kg/kg)
     ri = ice mixing ratio (kg/kg)

     If these are not given, the liquid and solid are divided according
     to the saturation adjustment scheme. See saturation_adjustment.m.

     If rt is replaced by 'sat', the saturation entropy is given
    """
    # stupid copy for python objects
    T,p,rt = T_in.copy(),p_in.copy(),rt_in.copy()

    rv,rl,ri,junk = saturation_adjustment(p,T,rt,c=c)
    
    condensed_input = 0

    # if saturation entropy
    if type(rt) == str:
        rv = r_sat(T,p,c=c)
        rl = np.zeros(rv.shape)
        ri = np.zeros(rv.shape)
        condensed_input = 0
    else:
        rv,rl,ri,junk = saturation_adjustment(p,T,rt,c=c)
    
    if 'rl' in kwargs.keys():
        rl = kwargs['rl']
        condensed_input = 1
    if 'ri' in kwargs.keys():
        ri = kwargs['ri']
        condensed_input = 1

    if condensed_input == 1:
        rv = rt - rl - ri

    e = p*rv/(c.eps+rv)
    pd = p-e

    # calculate specific entropy (dry, vapour, liquid, ice) 
    s_d = c.cp   * np.log(T/c.T0) - c.Rd* np.log(pd/c.p00)
    s_v = c.cpv  * np.log(T/c.T0) - c.Rv* np.log(e/c.e0) + c.Lv0/c.T0
    s_l = c.cpl  * np.log(T/c.T0) - c.Lv0/c.T0         + c.Lv0/c.T0
    s_i = c.cpi  * np.log(T/c.T0) - c.Ls0/c.T0         + c.Lv0/c.T0
    
   
    # Kludge to deal with zero vapour pressure
    s_v[rv<1e-10] = 0 
    
    s = s_d + rv*s_v + rl*s_l + ri*s_i

    
    
    return s    
