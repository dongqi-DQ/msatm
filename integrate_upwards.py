from calculate_dTdp_adiabatic import calculate_dTdp_adiabatic
from simple_fallout import simple_fallout
from saturation_adjustment import saturation_adjustment
from load_constants import load_constants
import pdb

def integrate_upwards(Tk_in,rtk_in,pk_in,dp_in,c=load_constants('default'),gamma=0):
    """
     Integrate the equation for enthalpy upwards. 
     1) 2nd order Runge-Kutta while conserving total water
     2) precipitation fall out by fraction gamma
    
     If gamma == 1, then the calculation of dT/dp is done assuming no water,
     and no precipitation fallout is required.
    """
    # stupid copy for python reasons
    Tk,rtk,pk,dp = Tk_in.copy(),rtk_in.copy(),pk_in.copy(),dp_in.copy()
    
    # calculate predictor
    dTdp,[rvk,junk] = calculate_dTdp_adiabatic(Tk,rtk,pk,c=c,gamma=gamma)
    
    T1 = Tk + dTdp*dp/2
    p1 = (pk+dp/2)

    rv1,rl1,ri1,junk = saturation_adjustment(p1,T1,rtk,c=c)
    nonvapour = rl1 + ri1
    if nonvapour.sum() > 0:
        rl1,ri1 = simple_fallout(T1,rvk-rv1,rl1,ri1,c=c,gamma=gamma)
   
    rt1 = rv1+rl1+ri1
    
    # calculate corrector 
    dTdp,[junk1,junk2] = calculate_dTdp_adiabatic(T1,rt1,p1,c=c,gamma=gamma)
    
    Tkp1 = Tk + dTdp*dp
    pkp1 = pk+dp
    rvkp1,rlkp1,rikp1,junk = saturation_adjustment(pkp1,Tkp1,rtk,c=c)
    nonvapour = rlkp1 + rikp1
    if nonvapour.sum() > 0:
        rlkp1,rikp1 = simple_fallout(Tkp1,rvk-rvkp1,rlkp1,rikp1,c=c,gamma=gamma)
    
    rtkp1 = rvkp1+rlkp1+rikp1
    
    return Tkp1,rtkp1,rvkp1,rlkp1,rikp1    
    

