import numpy as np

def calculate_frac_ice(T,varargin):
    """
    # Calculate the fraction of liquid and ice for a saturation adjustment
    #
    # [fliq,fice,[dfliqdT,dficedT]] = calculate_frac_ice(T,[type,ice,deltaT])
    #
    """
    c = load_constants(varargin)

    if c.ice == 0: 
        fliq = np.ones(T.shape)
        fice = np.zeros(T.shape)
        dfliqdT = 0
        dficedT = 0
    else:
       fice = -( T-(c.T0) )./(c.deltaT)
       fice(fice<0) = 0
       fice(fice>1) = 1
       fliq = 1-fice
       
    #if nargout > 2:
    dfliqdT = np.zeros(T.shape)
    if c.ice != 0:
      dfliqdT[(T>c.T0-c.deltaT) & (T<=c.T0)] = 1./c.deltaT

    varargout = [dfliqdT,-dfliqdT]

    return fliq,fice,varargout
