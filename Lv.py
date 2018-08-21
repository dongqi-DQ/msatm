from load_constants import load_constants 
def Lv(T,varargin):
    """
    Calculate the latent heats of vaporization 
    and freezing at given Temperature
    """
    c = load_constants(varargin)
    Lcond = c.Lv0 + (c.cpv - c.cpl )*(T-c.T0)
    Ldep = c.Ls0 + (c.cpv - c.cpi )*(T-c.T0)
    Lfrz = Ldep-Lcond

    return Lcond,Lfrz,Ldep


