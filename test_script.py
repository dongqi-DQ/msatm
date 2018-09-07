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
# test for adiabats

Ts1000 = np.array([300.])
rs1000 = np.array([0.02])
ps1000 = np.array([100000])
p1000 = np.linspace(100000,10000,101)
myc = load_constants('default')

T_r1000,rv1000,rl_r1000,ri_r1000,T_rho_r1000 = calculate_adiabat(Ts1000,rs1000,ps1000,p1000,gamma=0,c=myc)
T_p1000,rv1000,rl_r1000,ri_r1000,T_rho_p1000 = calculate_adiabat(Ts1000,rs1000,ps1000,p1000,gamma=1,c=myc)
T_h1000,rv1000,rl_r1000,ri_r1000,T_rho_h1000 = calculate_adiabat(Ts1000,rs1000,ps1000,p1000,gamma=0.5,c=myc)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
ax.plot(T_r1000,p1000/100.,'k',label='reversible')
ax.plot(T_p1000,p1000/100.,'r',label='pseudoadiabatic')
ax.plot(T_h1000,p1000/100.,'g',label='gamma=0.5')

ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Pressure (hPa)')
ax.invert_yaxis()
ax.legend()

fig.tight_layout()

fig.subplots_adjust(bottom=0.2)

ax3 = fig.add_axes([0.2,0.05,0.6,0.1])
ax3.axis('off')

caption = 'Fig. 1: Temperature for adiabatic parcel ascents assuming reversible (black) and pseudoadiabatic (red) thermodynamics with a mixed-phase range between 233.15 and 273.15 K.. Adiabat initialized with (T,r,p) = (%.2f K, %.2f kg/kg, %.2f hPa)' % (Ts1000,rs1000,p1000[0]/100)

ax3.text(0.5,0.5,caption, fontsize=10, style='oblique', ha='center', va='top', wrap=True)

fig.show()
pdb.set_trace()
fig.savefig('test_adiabat1.pdf')
