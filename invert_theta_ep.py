from msatm.load_constants import load_constants
import numpy as np
from msatm.r_sat import r_sat
from msatm.calculate_theta_ep import calculate_theta_ep
import pdb

def invert_theta_ep(theta_ep_in,r_in,Tstar_in,p_in,c=load_constants('bolton')):
    """
    Invert the pseudo-equivelant potetial temperature
    using iteraction procedure. Should converge within a few iterations
    Use Bolton formula, taken from Emanuel (1994), p 132
    """
    
    # stupid copy for python reasons
    theta_ep,r,Tstar,p = theta_ep_in.copy(),r_in.copy(),Tstar_in.copy(),p_in.copy() 
    
    c.ice = 0
    r_out = r.copy() 
    # No condensation: simply calculate T assuming Tstar and r are unchanged.
    T = theta_ep/ (  (c.p00/p)**(0.2854*(1-0.28*r))*np.exp( r * (1+0.81*r) * (3376.0/Tstar -2.54) )   )  

    # Find out where the super saturation is
    rsat = r_sat(T,p,c=c)
    saturated = rsat<r
    if saturated.sum()>0:
        Tsat = T[saturated]
        psat = p[saturated]
        rsat = rsat[saturated]
        Tstarsat = Tstar[saturated]

        pisat = (c.p00/psat)**(0.2854)

        # Begin iteration to fix supersaturation
        i_iter = 0
        dT = 1
        while dT > 0.01: # Assume we are converged if error in theta_ep < 0.01 K.
            # Iteration number
            i_iter = i_iter+1

            # Calculate error in theta_ep
            theta_ep_shoot,Tstarsat = calculate_theta_ep(Tsat,rsat,psat)
            dtheta_ep = theta_ep[saturated] - theta_ep_shoot

            # Calculate partial derivatives (simply)
            dthepdT = (  pisat*np.exp( rsat * (3376./Tstarsat -2.54) )   )  
            dthepdr = 3376./Tstarsat*Tsat*(  pisat*np.exp( rsat *(3376./Tstarsat -2.54) )   )  
            drdT = rsat*c.Lv0/(c.Rv*Tstarsat**2)

            dTdthep = 1./(dthepdT + dthepdr*drdT)

            # Find new Temperature in saturated regions
            Tnew = Tsat + 0.8*dtheta_ep*dTdthep


            # Calculate difference
            dT = np.max(np.abs(dtheta_ep))

            if i_iter > 12: 
                Tnew[np.abs(dtheta_ep)>100] = np.nan
            Tsat = Tnew
            rsat = r_sat(Tsat,psat,c=c)

            if i_iter > 12: 
                print(['iter: ' +str(i_iter)+ ', error: ' +str(dT) +' K'])

        T[saturated] = Tsat
        r_out[saturated] = rsat

    return T,r_out

