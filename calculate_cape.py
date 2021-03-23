# using marty's code, calculate CAPE and CIN
from msatm.calculate_adiabat import calculate_adiabat
import xarray as xr
import numpy as np
from msatm.saturation_adjustment import saturation_adjustment
from msatm.load_constants import load_constants
from scipy.interpolate import griddata,interp1d
import pdb
import matplotlib.pyplot as plt
from msatm.regridder import regrid_1d

c = load_constants('default')
# open a test case
f = xr.open_dataset('/g/data/ua8/Martin/va_analysis/syn599/CPOL_large-scale_forcing.nc')
p = f['lev'].data*100.
T = f['T'][80:100,:].data
r = f['r'][80:100,:].data/1000.
Ts = f['T_srf'][80:100].data + 273.15
rs = f['r_srf'][80:100].data/1000.
ps = np.ones(rs.shape)*p[0]

# get the env mix of phases
rv1,rl1,ri1,rs1 = saturation_adjustment(p,T,r,c=c)
# calculate the environmental density temperature
T_rho = T*(1 + rv1/c.eps)/(1+rv1+rl1+ri1)

#######################################################
# calculate the moist adiabat
# first calculate moist adiabat on high res p (5hPa?)
ptop = p[-1]
nsteps = int((ps[0]-ptop)/500.)
highp = np.linspace(ps[0],ptop,nsteps,endpoint=True)

# reversible
c.ice=1
mygam=0
T_r, rv, rl_r,ri_r,T_rho_r = calculate_adiabat(Ts,rs,ps,highp,gamma=mygam,c=c)
# # pseudoadiabatic
# c.ice=1
# mygam=1
# T_p, rv, rl_p,ri_p,T_rho_p = calculate_adiabat(Ts,rs,ps,p,gamma=mygam,c=c)
# # partially reversible
# c.ice=1
# mygam=0.5
# T_h, rv, rl_h,ri_h,T_rho_h = calculate_adiabat(Ts,rs,ps,p,gamma=mygam,c=c)

# resample T_r,rv,rl_r,T_rho_r from highp to coarse p
parcel_T_rho = regrid_1d(T_rho_r,highp,p,0) 

#######################################################
# calculate cape as the integral of Rd*(T_parcel - T_env) dln(p)
coarse_b = c.Rd*(parcel_T_rho - T_rho.T)
fine_b = regrid_1d(coarse_b,p,highp,0)

# multiply by ln(p)
highp_arr  = np.tile(np.expand_dims(highp,axis=1),(fine_b.shape[1]))
buoyancy = fine_b * np.log(highp_arr)

pos_b = buoyancy[...]
pos_b[pos_b<0] = 0
cape = pos_b.sum(axis=0)

neg_b = buoyancy[...]
neg_b[neg_b>0] = 0
cin = neg_b.sum(axis=0)




