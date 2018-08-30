from load_constants import load_constants
from e_sat import e_sat
from calculate_frac_ice import calculate_frac_ice

def saturation_adjustment(p,T,r_t,varargin):
    """
    # Function to calculate the breakdown of water into vapor, liquid and solid 
    # based on different microphysical assumptions
    #
    #  [r,rl,ri,[rs]] = saturation_adjustment(p,T,r_t [type,ice,deltaT])
    #
    # Outputs: 
    # r = water vapor mixing ratio (kg/kg)
    # rl = liquid water mixing ratio (kg/kg)
    # ri = solid water mixing ratio (kg/kg)
    # rs = saturation mixing ratio (kg/kg)
    #
    # Inputs:
    # p = pressure (Pa)
    # T = temperature (K)
    # r_t = total water mixing ratio (kg/kg)
    # deltaT = mixed phase range (K) (default = 40)
    """

    c = load_constants(varargin)
    es = e_sat(T,varargin)
    rs = c.eps*(es/(p-es))
    r = np.min(r_t,rs)
    r(r<0) = 0
    [fliq,fice,junk] = calculate_frac_ice(T,varargin)
    

    rl = fliq*(r_t-r)
    ri = fice*(r_t-r)

    rl(np.abs(r_t-r)<1e-8) = 0
    ri(np.abs(r_t-r)<1e-8) = 0
    
    #if nargout > 3 varargout{1} = rs end
    varargout = [rs]

    return r,rl,ri,varargout 
