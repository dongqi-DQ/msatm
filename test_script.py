# 7/9/2018 sugsnarsey@gmail.com
# perform some basic tests to check python version
# of marty singh's +atm module.

from load_constants import load_constants
import numpy as np
from e_sat import e_sat
import matplotlib.pyplot as plt
from load_constants import load_constants
from calculate_adiabat import calculate_adiabat
import pdb
from calculate_entropy import calculate_entropy
from saturation_adjustment import saturation_adjustment
from calculate_theta_ep import calculate_theta_ep
from invert_theta_ep import invert_theta_ep
from scipy.optimize import fsolve

####################################################################
# # test sat vapour pressure
# 
# T = np.linspace(201,300,100)
# 
# e_sat_cc,[esl,esi] = e_sat(T,c=load_constants('default'))
# e_sat_b,junk = e_sat(T,c=load_constants('bolton'))
# e_sat_t,junk = e_sat(T,c=load_constants('teten'))
# 
# fig = plt.figure(figsize=(10,5))
# ax1 = fig.add_subplot(121)
# ax1.plot(T,e_sat_cc,'k',label='e_sat')
# ax1.plot(T,esl,'b:',label='over liquid')
# ax1.plot(T,esi,'r:',label='over ice')
# 
# ax1.set_yscale('log')
# 
# ax1.set_ylabel('vapour pressure (Pa)')
# ax1.set_xlabel('Temperature (K)')
# ax1.legend()
# 
# ax2 = fig.add_subplot(122)
# ax2.plot(T,e_sat_b/e_sat_cc,'r',label='bolton')
# ax2.plot(T,e_sat_t/e_sat_cc,'g',label='teten')
# 
# ax2.set_ylabel('vapour pressure ratio')
# ax2.set_xlabel('Temperature (K)')
# 
# ax2.legend()
# 
# fig.tight_layout()
# 
# fig.subplots_adjust(bottom=0.2)
# 
# ax3 = fig.add_axes([0.2,0.05,0.6,0.1])
# ax3.axis('off')
# caption = 'Fig. 1 (Left) saturation vapor pressure uisng constant c_p method (black) and saturation vapor pressure over liquid (blue) and ice (red). \n (Right) Ratio of saturation vapor pressure calculated using Bolton (red) and Teten (green) methods compared to the default method (const. c_p).'
# 
# ax3.text(0.5,0.5,caption, fontsize=10, style='oblique', ha='center', va='top', wrap=True)
# 
# 
# fig.show()
# pdb.set_trace()
# 
# fig.savefig('test_e_sat.pdf')
####################################################################
# ##test for adiabats with a dry adiabatic layer as well
# 
# Ts1000 = np.array([300.])
# rs1000 = np.array([0.02])
# ps1000 = np.array([100000])
# p1000 = np.linspace(100000,10000,101)
# myc = load_constants('default')
# 
# T_r1000,rv1000,rl_r1000,ri_r1000,T_rho_r1000 = calculate_adiabat(Ts1000,rs1000,ps1000,p1000,gamma=0,c=myc)
# T_p1000,rv1000,rl_r1000,ri_r1000,T_rho_p1000 = calculate_adiabat(Ts1000,rs1000,ps1000,p1000,gamma=1,c=myc)
# T_h1000,rv1000,rl_r1000,ri_r1000,T_rho_h1000 = calculate_adiabat(Ts1000,rs1000,ps1000,p1000,gamma=0.5,c=myc)
# 
# 
# fig = plt.figure(figsize=(10,10))
# ax = fig.add_subplot(111)
# ax.plot(T_r1000,p1000/100.,'k',label='reversible')
# ax.plot(T_p1000,p1000/100.,'r',label='pseudoadiabatic')
# ax.plot(T_h1000,p1000/100.,'g',label='gamma=0.5')
# 
# ax.set_xlabel('Temperature (K)')
# ax.set_ylabel('Pressure (hPa)')
# ax.invert_yaxis()
# ax.legend()
# 
# fig.tight_layout()
# 
# fig.subplots_adjust(bottom=0.2)
# 
# ax3 = fig.add_axes([0.2,0.05,0.6,0.1])
# ax3.axis('off')
# 
# caption = 'Fig. 1: Temperature for adiabatic parcel ascents assuming reversible (black) and pseudoadiabatic (red) thermodynamics with a mixed-phase range between 233.15 and 273.15 K.. Adiabat initialized with (T,r,p) = (%.2f K, %.2f kg/kg, %.2f hPa)' % (Ts1000,rs1000,p1000[0]/100)
# 
# ax3.text(0.5,0.5,caption, fontsize=10, style='oblique', ha='center', va='top', wrap=True)
# 
# fig.show()
# pdb.set_trace()
# fig.savefig('test_adiabat1.pdf')

####################################################################
# # compare cases against conservation of entropy - no ice (conserve equivalent potential temperature). 
# 
# Ts = np.array([298.15])
# rs = np.array([0.02])
# ps = np.array([95000])
# p = np.linspace(95000,10000,96)
# c = load_constants('default')
# 
# # no ice reversible
# c.ice=0
# mygam=0
# T_noice_r, rv, rl_noice_r,ri_noice_r,T_rho_noice_r = calculate_adiabat(Ts,rs,ps,p,gamma=mygam,c=c)
# # no ice pseudoadiabatic
# c.ice=0
# mygam=1
# T_noice_p, rv, rl_noice_p,ri_noice_p,T_rho_noice_p = calculate_adiabat(Ts,rs,ps,p,gamma=mygam,c=c)
# 
# # reference case (no ice, conserve entropy)
# c = load_constants('default')
# c.ice = 0
# 
# ss = calculate_entropy(Ts,ps,rs,c=c)
# T_entropy = np.zeros(p.shape)
# T_entropy[0] = Ts[...]
# 
# for kk in range(1,len(p)):
#     T_entropy[kk] = fsolve( lambda x: ss - calculate_entropy(x,np.array([p[kk]]),rs,c=c),T_entropy[kk-1])
# 
# rv,rl,ri,junk = saturation_adjustment(p,T_entropy,rs*np.ones(p.shape),c=c)
# T_rho_entropy = T_entropy*(1 + rv/c.eps)/(1+rv+rl+ri)
# 
# theta_ep,Tstar = calculate_theta_ep(Ts,rs,ps)
# 
# T_psu_entropy = np.zeros(p.shape)
# rv_psu_entropy = np.zeros(p.shape)
# T_psu_entropy[0] = Ts[...]
# rv_psu_entropy[0] = rs[...]
# 
# for kk in range(1,len(p),1):
#     T_psu_entropy[kk],rv_psu_entropy[kk] = invert_theta_ep(theta_ep,rs,Tstar,np.array([p[kk]])) 
# 
# T_rho_psu_entropy = T_psu_entropy*(1 + rv_psu_entropy/c.eps)/(1+rv_psu_entropy)
# 
# fig = plt.figure(figsize=(15,8))
# ax1 = fig.add_subplot(121)
# ax1.plot(T_noice_r - T_entropy,p/100.,'k',label='reversible')
# ax1.plot(T_noice_p - T_entropy,p/100.,'r',label='pseudoadiabatic')
# ax1.plot(T_psu_entropy - T_entropy,p/100.,'b',label='Bolton-pseudo')
# ax1.set_xlabel('Temperature (K)')
# ax1.set_ylabel('Pressure (hPa')
# ax1.invert_yaxis()
# ax1.legend()
# ax2 = fig.add_subplot(122)
# ax2.plot(T_rho_noice_r - T_rho_entropy,p/100.,'k',label='reversible')
# ax2.plot(T_rho_noice_p - T_rho_entropy,p/100.,'r',label='pseudoadiabatic')
# ax2.plot(T_rho_psu_entropy - T_rho_entropy,p/100.,'b',label='Bolton-pseudo')
# ax2.set_xlabel('Density temperature (K)')
# ax2.invert_yaxis()
# 
# fig.subplots_adjust(bottom=0.2)
# 
# ax3 = fig.add_axes([0.2,0.05,0.6,0.1])
# ax3.axis('off')
# 
# caption = 'Fig. 3: Temperature (left) and density temperature (right) for adiabatic parcel ascents assuming reversible (black) and pseudoadiabatic (red) thermodynamics and assuming no ice formation. The profiles are plotted as an anomaly from a control ascent that is calculated based on conservation of entropy, assuming no ice formation (an exact solution for the reversible case). Approximate pseudoadiabatic ascent calculated by assuming conservation of pseudo-equivalent potential temperature as defined by Bolton (1980) is shown in blue. Adiabat initialized with (T,r,p) = (%.2f K, %.2f kg/kg, %.2f hPa)' % (Ts,rs,p[0]/100.)
# 
# ax3.text(0.5,0.5,caption, fontsize=10, style='oblique', ha='center', va='top', wrap=True)
# 
# fig.show()
# pdb.set_trace()
# fig.savefig('test_adiabat2.pdf')

####################################################################
# ## again, compare a few cases against reference case
# Ts = np.array([298.15])
# rs = np.array([0.02])
# ps = np.array([95000])
# p = np.linspace(95000,10000,96)
# c = load_constants('default')
# 
# # reversible
# c.ice=1
# mygam=0
# T_r, rv, rl_r,ri_r,T_rho_r = calculate_adiabat(Ts,rs,ps,p,gamma=mygam,c=c)
# # pseudoadiabatic
# c.ice=1
# mygam=1
# T_p, rv, rl_p,ri_p,T_rho_p = calculate_adiabat(Ts,rs,ps,p,gamma=mygam,c=c)
# # partially reversible
# c.ice=1
# mygam=0.5
# T_h, rv, rl_h,ri_h,T_rho_h = calculate_adiabat(Ts,rs,ps,p,gamma=mygam,c=c)
# 
# # reference case (no ice, conserve entropy)
# c = load_constants('default')
# c.ice = 0
# 
# ss = calculate_entropy(Ts,ps,rs,c=c)
# T_entropy = np.zeros(p.shape)
# T_entropy[0] = Ts[...]
# 
# for kk in range(1,len(p)):
#     T_entropy[kk] = fsolve( lambda x: ss - calculate_entropy(x,np.array([p[kk]]),rs,c=c),T_entropy[kk-1])
# 
# rv,rl,ri,junk = saturation_adjustment(p,T_entropy,rs*np.ones(p.shape),c=c)
# T_rho_entropy = T_entropy*(1 + rv/c.eps)/(1+rv+rl+ri)
# 
# 
# fig = plt.figure(figsize=(20,8))
# ax1 = fig.add_subplot(131)
# ax1.plot(T_r - T_entropy,p/100.,'k',label='reversible')
# ax1.plot(T_p - T_entropy,p/100.,'r',label='pseudoadiabatic')
# ax1.plot(T_h - T_entropy,p/100.,'g',label='gamma=0.5')
# ax1.set_xlabel('Temperature (K)')
# ax1.set_ylabel('Pressure (hPa')
# ax1.invert_yaxis()
# ax1.legend()
# 
# 
# ax2 = fig.add_subplot(132)
# ax2.plot(T_rho_r - T_rho_entropy,p/100.,'k',label='reversible')
# ax2.plot(T_rho_p - T_rho_entropy,p/100.,'r',label='pseudoadiabatic')
# ax2.plot(T_rho_h - T_rho_entropy,p/100.,'g',label='gamma=0.5')
# ax2.set_xlabel('Density temperature (K)')
# ax2.invert_yaxis()
# 
# 
# ax3 = fig.add_subplot(133)
# ax3.plot(rl_r,p/100.,'k',linestyle=':',label='liquid')
# ax3.plot(ri_r,p/100.,'k',linestyle='--',label='ice')
# ax3.plot(ri_r+rl_r,p/100.,'k',linestyle='-',label='total')
# 
# ax3.plot(rl_p,p/100.,'r',linestyle=':')
# ax3.plot(ri_p,p/100.,'r',linestyle='--')
# ax3.plot(ri_p+rl_p,p/100.,'r',linestyle='-')
# 
# ax3.plot(rl_h,p/100.,'g',linestyle=':')
# ax3.plot(ri_h,p/100.,'g',linestyle='--')
# ax3.plot(ri_h+rl_h,p/100.,'g',linestyle='-')
# 
# ax3.set_xlabel('mixing ratio (kg/kg)')
# ax3.invert_yaxis()
# ax3.legend()
# 
# fig.subplots_adjust(bottom=0.2)
# 
# ax4 = fig.add_axes([0.2,0.05,0.6,0.1])
# ax4.axis('off')
# 
# caption = 'Fig. 4: Temperature (left), density temperature (middle), and condensate mixing ratio (right) for parcel ascents including ice and using reversible (blue) pseudo-adiabatic (red) and intermediate thermodynamics. The intermediate case is calculated by assuming half of the condensate created at a given level is precipitated out. The temperature profiles are plotted as ananomaly from a control parcel ascent calculated assuming conservation of entropy with no ice. Adiabat initialized with (T,r,p) = (%.2f K, %.2f kg/kg, %.2f hPa)' % (Ts,rs,p[0]/100.)
# 
# ax4.text(0.5,0.5,caption, fontsize=10, style='oblique', ha='center', va='top', wrap=True)
# 
# fig.show()
# pdb.set_trace()
# fig.savefig('test_adiabat3.pdf')


###################################################################
# test array calcs

Ts = np.linspace(1,1.1,12).reshape((4,3))*290.
rs = np.ones((4,3))*0.02
ps = np.ones((4,3))*95000.
p = np.linspace(95000,10000,96)
c = load_constants('default')
c.ice = 1
gamma = 0

T,rv,rl,ri,T_rho = calculate_adiabat(Ts,rs,ps,p,gamma=gamma,c=c)

for x in range(4):
    for y in range(3):
        plt.plot(T[:,x,y],p,label='Ts = '+str(np.round(Ts[x,y],2)))

plt.gca().invert_yaxis()
plt.legend()
plt.show()


