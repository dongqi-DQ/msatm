from calculate_frac_ice import calculate_frac_ice

def simple_fallout(T,cond,rl,ri,gamma,varargin):
    """
    # A very simple precipitation parameterization for use in parcel ascent calculations.
    # Assumes a fraction gamma of the condensation is reomved from the parcel
    #
    # inputs:
    #		T	temperature (K)
    #		cond	amount of condensation (kg condesate / kg dry air)
    #		rl	liquid water mixing ratio (kg/kg)
    #		ri	ice water mixing ratio (kg/kg)
    #		gamma   fraction to precipitate (0-1)
    """

    # Break the condensation into ice and liquid
    fliq,fice,junk = calculate_frac_ice(T,varargin)
 
    # Remove a fraction gamma
      rl_new = max(0,rl - gamma*(cond)*fliq)
      ri_new = max(0,ri - gamma*(cond)*fice)

    return rl_new,ri_new
