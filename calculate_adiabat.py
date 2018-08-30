def calculate_adiabat(Ts,rts,ps,p,varargin):
    """
    # Function to calculate an adiabatic parcel ascent with various thermodynamic assumptions
    #
    # [T,rv,rl,ri,T_rho] = calculate_adiabat(Ts,rts,ps,p,[gamma,type,ice,deltaT])
    #
    # INPUTS:
    #
    # Ts = parcel temperature (K)
    # rts = parcel total mixing ratio (kg/kg)
    # ps = Parcel pressure (Pa)
    # p = pressure levels to integrate parcel (Pa)
    #
    # gamma = fraction of condensate removed from parcel
    # type, ice, deltaT are microphysical parameters
    # See e_sat.m for details.
    #
    # OUTPUTS:
    #
    # T = temperature of parcel at all pressure levels (K)
    # rv = water vapor mixing raio of parcel (kg/kg)
    # rl = liquid water mixing ratio of parcel (kg/kg)
    # ri = solid water mixing ratio of parcel (kg/kg)
    # T_rho = density temperature of parcel (K)
    #
    # p must either be a vector or of size [Np size(Ts)] and should be sorted largest to smallest
    # Ts, ps, and rts may be up to 3 dimensions
    # The outputs are of size [Np size(Ts)]
    #
    # Parcel is lifted from ps to the top of the pressure matrix
    # all values of p greater than  ps are filled with Nans
    #
    # This is helpful to vectorize the parcel ascent calculation over a 
    # matrix of different columns which may have different surface pressures
    #
    #
    # The adiabat is calculated in three steps.
    #
    # 1) The LCL pressure is calculated using the method of Romps (2017).
    #    This method provides an analytic exact solution for the LCL in terms
    #    of the Lambert W function. The MATLAB library for the W function must
    #    be loaded, so this part of the calculation can be slow the first time
    #    the adiabat calculator is called. 
    #
    #    If the LCL lies in the mixed-phase range, the calculation is no longer
    #    analytic, and an iteration is required. See calculate_LCL.m for more
    #    details
    #   
    # 2) For levels below the LCL, the parcel is assumed to conserve its
    #    potential temperature theta = T*(p/p_LCL)^(Rm/cpm) and its mixing
    #    ratio r. This gives the temperature by simple rearrangement
    #
    # 3) Above the LCL, the parcel temperature is calculated by by integrating 
    #    the equation for the lapse rate:
    #
    #                       dT/dp = f(T,rt,p) 
    #
    #    See the function calculate_dTdp_adiabatic.m for more details
    #
    #    The ascent is performed by calculating dT/dp assuming no fallout and
    #    then removing a fraction gamma of the water condensed during that step
    #    isothermally. In the case of gamma == 1, dT/dp is calculated assuming 
    #    no condensed water for pseudoadiabatic ascent.
    #
    #
    #
    # The accuracy of these calculations may be checked using the script
    # test_script.m. For the case of no ice and the default thermodynamics,
    # an exact solution for a reversible parcel ascent is conservation of
    # entropy. This script is accurate to within 0.5 K for dp = 50 hPa, and to
    # within 0.1 K for dp = 20 hPa.
    # 
    """




    # Read in the optional arguments
    # These modify the default characteristics set in load_constants.m

    c = atm.load_constants(varargin[1:])
    gamma = c.gamma
    #if nargin >= 5 gamma  = varargin{1}  end

    # Arrangement of the columns
    xygrid = Ts.shape

    # Make pressure interval matrix
    if len(p.shape) == 1: 	# If the pressure is a vector
        p = p[:]
        p = np.tile(p,(*xygrid,1)).T
      
    dp = np.diff(p,n=1,axis=0)

    # persistent big_dp
    # if max(abs(dp(:))) >= 5000 & isempty(big_dp)
    #    warning(['calculate_adiabat.m: You have pressure intervals greater than or equal to 50 hPa. ' ...
    #             'It is recommened that you subsample the pressure matrix for better accuracy of the integration.'])
    #    big_dp = 1
    # end

    # Initialize profiles
    T  = np.zeros((p.shape[0],*xygrid))*np.nan #nan([size(p,1) xygrid ]) # temperature
    rv = np.zeros((p.shape[0],*xygrid))*np.nan #nan([size(p,1) xygrid ]) # vapor
    rl = np.zeros((p.shape[0],*xygrid))*np.nan #nan([size(p,1) xygrid ]) # liquid
    ri = np.zeros((p.shape[0],*xygrid))*np.nan #nan([size(p,1) xygrid ]) # solid
    rt = np.zeros((p.shape[0],*xygrid))*np.nan #nan([size(p,1) xygrid ]) # total water

    # Initialize surface values
    [rvs,junk1,junk2,junk3] = saturation_adjustment(ps,Ts,rts,[c.modeltype,c.ice,c.deltaT])
    Rm = c.Rd + c.Rv*rvs
    cpm = c.cp + c.cpv*rvs

    if gamma ==1:
        rts = rvs 
   
# Sugata is up to here!


[T_LCL,p_LCL] = atm.calculate_LCL(Ts,rts,ps,varargin{2:end})

Im_prev = zeros(size(Ts))
for k = 1:size(p,1)

    # Pressure at this level
    pk = squeeze(p(k,:,:,:))

    # Find regions where the surface pressure is larger than the current level p(k)
    Ia = ps>=reshape(p(k,:,:,:),xygrid)
    
    # Find regions where the LCL pressure is smaller than the current level p(k)
    Id = Ia & p_LCL<=reshape(p(k,:,:,:),xygrid)
    
    # Find the temperature in the regions below the LCL
    T(k,Id) = T_LCL(Id).*(pk(Id)./p_LCL(Id)).^(Rm(Id)./cpm(Id))
    rt(k,Id) = rts(Id)
    rv(k,Id) = rvs(Id)
    rl(k,Id) = 0
    ri(k,Id) = 0
    
    # Find regions where the LCL pressure is smaller than the current level p(k)
    Im = Ia & p_LCL>reshape(p(k,:,:,:),xygrid)
    
   
    # Find regions where the LCL is between levels (k-1) and k
    Is = Im & ~Im_prev
 
    ## Integrate upwards (towards lower pressure) from the surface to the current model level for the regions where the current level p(k) is just above the surface

    # Calculate the dp
    if sum(Is(:))>0
      dps = reshape(p(k,:,:,:),xygrid) - p_LCL
      [T(k,Is),rt(k,Is),rv(k,Is),rl(k,Is),ri(k,Is)] = integrate_upwards(T_LCL(Is),rts(Is),p_LCL(Is),dps(Is),gamma,c.type,c.ice,c.deltaT)
    end
    
    
    ## Now integrate from level k to level k+1
    if k < size(p,1)
       [T(k+1,Im),rt(k+1,Im),rv(k+1,Im),rl(k+1,Im),ri(k+1,Im)] = integrate_upwards(T(k,Im),rt(k,Im),p(k,Im),dp(k,Im),gamma,c.type,c.ice,c.deltaT)
    end
        
    Im_prev = Im

end



T_rho = T.*(1 + rv./c.eps)./(1+rv+rl+ri)




end
    return T,rv,rl,ri,T_rho 




function [Tkp1,rtkp1,rvkp1,rlkp1,rikp1] = integrate_upwards(Tk,rtk,pk,dp,gamma,type,ice,deltaT)
#
# Integrate the equation for enthalpy upwards. 
# 1) 2nd order Runge-Kutta while conserving total water
# 2) precipitation fall out by fraction gamma
#
# If gamma == 1, then the calculation of dT/dp is done assuming no water,
# and no precipitation fallout is required.
#
    # calculate predictor
    [dTdp,rvk] = atm.calculate_dTdp_adiabatic(Tk,rtk,pk,gamma,type,ice,deltaT)
    
    T1 = Tk + dTdp.*dp./2
    p1 = (pk+dp./2)
    [rv1,rl1,ri1] = atm.saturation_adjustment(p1,T1,rtk,type,ice,deltaT)
    if rl1+ri1 > 0
       [rl1,ri1] = atm.simple_fallout(T1,rvk-rv1,rl1,ri1,gamma,type,ice,deltaT)
    end
    rt1 = rv1+rl1+ri1
    
    # calculate corrector 
    dTdp = atm.calculate_dTdp_adiabatic(T1,rt1,p1,gamma,type,ice,deltaT)
    
    Tkp1 = Tk + dTdp.*dp
    pkp1 = pk+dp
    [rvkp1,rlkp1,rikp1] = atm.saturation_adjustment(pkp1,Tkp1,rtk,type,ice,deltaT)
    if rlkp1+rikp1 > 0
       [rlkp1,rikp1] = atm.simple_fallout(Tkp1,rvk-rvkp1,rlkp1,rikp1,gamma,type,ice,deltaT)
    end
    rtkp1 = rvkp1+rlkp1+rikp1
    
    
end
    

