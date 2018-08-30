from load_constants import load_constants
from e_sat import e_sat
import numpy as np

def calculate_LCL(T,r,p,varargin):
    """
    #
    # Function to calculate the lifted condensation level (LCL) for a parcel.
    #
    #
    # [T_LCL,p_LCL,[T_LCLl,p_LCLl,T_LCLi,p_LCLi]] = calculate_LCL(T,r,p,[type,ice,deltaT])
    #
    #
    # Finds the temperature (T_LCL) and pressure (p_LCL) at which a parcel
    # with a given temperature (T), mixing ratio (r) and pressure (p) becomes
    # saturated.
    #
    # Phase changes are assumed to have a mixed phase range, with consitent
    # alterations to the saturation vapour pressure. These options are set by
    # the optional input arguments (see e_sat.m for details).
    #
    # Optional output arguments give the lifted condensation level (LCL liquid
    # water) and the lifted deposition level (LDL for ice).
    #
    # The parcel is assumed to conserve its potential temperature and its
    # mixing ratio as it is lifted. These assumptions allow the derivation of
    # an analytic expression for the pressure of the LCL in terms of Lambert W
    # functions. See Romps (2017) for details.
    #
    #
    # The MATLAB library for the Lambert W function appears to take some time to
    # load on the first call to the W function within a session. This means
    # that calling calculate_LCL will be rather slow on the first call, but
    # will speed up for subsequent calls. This should be taken into account
    # when testing functions that use calculate_LCL (for example
    # calculate_adiabat) for speed.
    #
    # If the LCL is within the mixed-phase range, the solution for the LCL is
    # no longer analytic, but an iteration must be performed. This has the
    # potenial to slow down the calculation further, but testing indicates that
    # only a few iterations are normally required.
    #
    """
    # Load constants
    c = load_constants(varargin)

    # Vapor dependent variables
    cpm = c.cp + r*c.cpv
    Rm = c.Rd + r*c.Rv

    # Vapor pressures
    e = p*r/(r+c.eps)
    es,esl,esi = e_sat(T,varargin)

    # Relative humidities
    RHl = e/esl
    RHi = e/esi

    # Calculate lifted condensation temperature
    ac = cpm/Rm + (c.cpl-c.cpv)/c.Rv
    bc = - ( c.Lv0 - (c.cpv-c.cpl)*c.T0 )/(c.Rv*T)
    cc = bc/ac

    # Calculate LCL (see Romps 2017)
    T_LCLl = cc*lambertw(-1,RHl**(1/ac)*cc*np.exp(cc))**(-1)*T


    # Calculate lifted deposition temperature
    ac = cpm/Rm + (c.cpi-c.cpv)/c.Rv
    bc = - ( c.Ls0 - (c.cpv-c.cpi)*c.T0 )/(c.Rv*T)
    cc = bc/ac

    # Calculate LDL (see Romps 2017)
    T_LCLi = cc*lambertw(-1,RHi**(1/ac)*cc*np.exp(cc))**(-1)*T

    ## Set the temperature of lifted saturation level
    # First set LCL to condensation level
    T_LCL = T_LCLl.copy()

    if c.ice > 0:

        # Below homogenous freezing temperature set T_LCL to deposition level
        T_LCL(T_LCLl < c.T0-c.deltaT)   = T_LCLi(T_LCLl < c.T0-c.deltaT)

        # For Temperatre of LCL is in mixed-phase range
        I = (T_LCLl > c.T0-c.deltaT) & (T_LCLl < c.T0)

        # Need to iterate to find correct value
        myiter = 0
        dT = np.ones(T_LCL.shape)

        while max(abs(dT))>0.02:
            [es,esl,esi] = e_sat(T_LCL[I])
            RHi_LCL = es/esi
            T_LCL_new = cc[I]*lambertw(-1,(RHi[I]/RHi_LCL)**(1/ac[I])*cc[I]*np.exp(cc[I]))**(-1)*T[I]

            dT = T_LCL_new-T_LCL[I]
            T_LCL[I] = T_LCL_new
            myiter += 1

    else:
        T_LCLi = T_LCLl

    p_LCL = p*(T_LCL/T)**(cpm/Rm)
    p_LCLl = p*(T_LCLl/T)**(cpm/Rm)
    p_LCLi = p*(T_LCLi/T)**(cpm/Rm)       
    
    varargout = [T_LCLl,p_LCLl,T_LCLi,p_LCLi]
    
    return T_LCL,p_LCL,varargout
