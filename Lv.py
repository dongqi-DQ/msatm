from msatm.load_constants import load_constants 

def Lv(T,c=load_constants('default')):
    """
    Calculate the latent heats of vaporization 
    and freezing at given Temperature
    """
    Lcond = c.Lv0 + (c.cpv - c.cpl )*(T-c.T0)
    Ldep = c.Ls0 + (c.cpv - c.cpi )*(T-c.T0)
    Lfrz = Ldep-Lcond

    return Lcond,Lfrz,Ldep


