from load_constants import load_constants
import numpy as np

def e_sat(T,varargin):
    """
    # Function ot calculate saturation vapor pressure
    #
    # [es[,esl,esi]] = e_sat(T[,type,ice,deltaT])
    #
    # es = saturation vapor pressure (Pa)
    # T = temperature (K)
    #
    # type = {'default','bolton','teten','sam'}
    # ice =  {[0] , [1] }
    # deltaT = mixed-phase temperature ranges
    """
    c = load_constants(varargin)

    ## Calculate saturation vapor pressure over liquid and solid
    if c.modeltype = 'default':
        # Thermodynamically consistent definition of saturation curves
        # assuming that the heat capacities of all consituents are
        # independent of temperature
        # i.e. integral of Clausius-Clapeyron equation with constant heat
        # capacities. See Romps (2008) for details.

        esl = c.e0*(T/c.T0).^((c.cpv-c.cpl)/c.Rv)*
             np.exp( ( c.Lv0 - c.T0*(c.cpv-c.cpl) )/c.Rv * ( 1/c.T0 - 1/T ) )

        esi = c.e0*(T/c.T0).^((c.cpv-c.cpi)/c.Rv)*
             np.exp( ( c.Ls0 - c.T0*(c.cpv-c.cpi) )/c.Rv * ( 1/c.T0 - 1/T ) )
    
    elif c.modeltype = 'bolton':
        # Bolton's formulas (Bolton, 1980) 
        # Used in George Bryan's CM1
        esl = 611.2*np.exp( 17.67      * ( T  - 273.15 ) / ( T  - 29.65 ) )
        esi = 611.2*np.exp( 21.8745584 * ( T  - 273.15 ) / ( T  - 7.66  ) )

    elif c.modeltype = 'teten':
        # Teten's formulas
        # Used in ECMWF
        # coefficients in Tetens saturation vapor pressure es = es0 * np.exp(a3 * (T-T0)/(T-a4))      
        es0       = 611.21  # saturation vapor pressure at T0 (Pa)
        T0        = 273.16  # (K)
        a3l       = 17.502  # liquid water (Buck 1981)
        a4l       = 32.19   # (K)
        a3i       = 22.587  # ice (Alduchov and Eskridge 1996)
        a4i       = -0.7    # (K)
        # saturation vapor pressure over liquid and ice
        esl       = es0 * np.exp(a3l * (T - T0)/(T - a4l))
        esi       = es0 * np.exp(a3i * (T - T0)/(T - a4i))

    elif c.modeltype = 'sam':
        # Formulation in SAM model
        # Probably not a good idea to use this without ice

        a0 = 6.11239921
        a1 = 0.443987641
        a2 = 0.142986287e-1
        a3 = 0.264847430e-3
        a4 = 0.302950461e-5
        a5 = 0.206739458e-7
        a6 = 0.640689451e-10
        a7 = -0.952447341e-13
        a8 = -0.976195544e-15
        dt = max(-80.,T-273.16)

        esl = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))

        esi = esl

        a0 = 6.11147274
        a1 =  0.503160820
        a2 = 0.188439774e-1
        a3 =  0.420895665e-3
        a4 =  0.615021634e-5
        a5 = 0.602588177e-7
        a6 = 0.385852041e-9
        a7 = 0.146898966e-11
        a8 = 0.252751365e-14
        dt = T[T>185]-273.16

        esi(T>185) = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))

        dt = max(-100.,T[T<=185]-273.16)
        esi[T<=185] = 0.00763685 + dt*(0.000151069+dt*7.48215e-07)


        esl = esl*100
        esi = esi*100

    elif c.modeltype == 'fms':
        # Simple thermodynamics with no ice
        # Ice should not be turned on if you use this option!

        e0 = 610.78 
        esl = c.e0 * np.exp(c.Lv0/c.Rv*(1.0/c.T0 - 1.0/T))
        esi = esl

    else:
        raise ValueError('Unknown thermodynamics type. Cannot calculate saturation vapor pressure')

    ## Calculate the mixture of liquid and solid
    fliq,fice,junk = calculate_frac_ice(T,varargin)


    es = esl*(fliq) + esi*fice

#     if nargout > 1 varargout{1} = esl end
#     if nargout > 2 varargout{2} = esi end

    varargout = [esl,esi]    

    return es,varargout





       
