from calculate_dTdp_adiabatic import calculate_dTdp_adiabatic
from simple_fallout import simple_fallout
from saturation_adjustment import saturation_adjustment
def integrate_upwards(Tk,rtk,pk,dp,gamma,modeltype,ice,deltaT):
    """
    # Integrate the equation for enthalpy upwards. 
    # 1) 2nd order Runge-Kutta while conserving total water
    # 2) precipitation fall out by fraction gamma
    #
    # If gamma == 1, then the calculation of dT/dp is done assuming no water,
    # and no precipitation fallout is required.
    """
    # calculate predictor
    dTdp,rvk,junk = calculate_dTdp_adiabatic(Tk,rtk,pk,gamma,modeltype,ice,deltaT)
    
    T1 = Tk + dTdp*dp/2
    p1 = (pk+dp/2)

    rv1,rl1,ri1,junk = saturation_adjustment(p1,T1,rtk,[modeltype,ice,deltaT])
    if rl1+ri1 > 0:
        rl1,ri1 = simple_fallout(T1,rvk-rv1,rl1,ri1,gamma,modeltype,ice,deltaT)
    
    rt1 = rv1+rl1+ri1
    
    # calculate corrector 
    dTdp,junk1,junk2 = calculate_dTdp_adiabatic(T1,rt1,p1,gamma,modeltype,ice,deltaT)
    
    Tkp1 = Tk + dTdp*dp
    pkp1 = pk+dp
    rvkp1,rlkp1,rikp1,junk = saturation_adjustment(pkp1,Tkp1,rtk,modeltype,ice,deltaT)
    if rlkp1+rikp1 > 0:
       rlkp1,rikp1 = simple_fallout(Tkp1,rvk-rvkp1,rlkp1,rikp1,gamma,modeltype,ice,deltaT)
    
    rtkp1 = rvkp1+rlkp1+rikp1
    
    return Tkp1,rtkp1,rvkp1,rlkp1,rikp1    
    

