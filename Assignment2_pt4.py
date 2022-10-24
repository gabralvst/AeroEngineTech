import numpy as np
from matplotlib import pyplot as plt

gross_th= 15167 #N
mdot_a = 23.81  #kg/s
mdot_f = 0.4267 #kg/s
tot_comp_ratio = 5.5
comb_eff = 1
comp_eff = 0.89
turb_eff = 0.8
nozz_eff = 1
gamma_air = 1.4
gamma_gas = 1.3
gas_const = 287 #J/kg/K
LHV = 43 * 10**6 #J
cp_a = 1000
cp_g = 1150
R = 287
h = 10668 #m
Mach = 0.78
TIT = 1150 #k
t_amb = 219
p_amb = 23911 #Pa
rho_amb = 0.3796
kh = 0.161
A_turbine = 0.5103267031956857 * 0.07679508070244695  #m2
A_nozzle = 0.0768 #m2
beta_tt = 8.10682594846814 #compressor pressure ratio
t_t1 = 245.64792 #K
p_t1 = 35738.733374570875 #Pa
reaction = 0.5
flow_coeff = 0.6
work_coeff = 0.4
u_max = 450 #m/s
Nstages = 8
RPM = 13800 #rot/min
V1 = 150.38
V2 = 206.68
mdot_engine = 12.947804546087896

t_t2 = t_t1 * beta_tt**((gamma_air-1)/gamma_air)
work = cp_a*(t_t2 - t_t1)
work_per_stage = work / Nstages
tip_speed = np.sqrt(work_per_stage / work_coeff)
deltaT = work_per_stage / cp_a
t_s1 = t_t1 - ((gamma_air - 1) / 2) * (V1 ** 2 / (gamma_air * R))

T_lst = [t_t1]
P_lst = [p_t1]
Ps_lst = []
Ts_lst = [t_s1]
M_lst = []
T_ratio_lst = []
Ts_ratio_lst = []
Psratio_lst = []
beta_lst = []
delta_s_lst = []
h_lst = []
x_axis = np.arange(0,Nstages*2,1) #delete *2 to plot pressure ratio per stage
for i in range(Nstages):
    Tt = T_lst[i] + deltaT
    T_lst.append(Tt)
    Ts = Tt - ((gamma_air - 1) / 2) * (V1 ** 2 / (gamma_air * R))
    Ts_lst.append(Ts)
    M = V1 / (np.sqrt(gamma_air * R * Ts_lst[i]))
    M_lst.append(M)
    Ttratio = Tt/T_lst[i]
    T_ratio_lst.append(Ttratio)
#   T_ratio_lst.append(Ttratio)
    Tsratio = Ts / Ts_lst[i]
    Ts_ratio_lst.append(Tsratio)
    beta = Ttratio**((gamma_air*comb_eff)/(gamma_air-1))
    beta_lst.append(beta)
    beta_lst.append(beta)
    Pt = P_lst[i] * beta_lst[i]
    P_lst.append(Pt)
    Ps = P_lst[i]*(1+(gamma_air-1)/2*M**2)**(-gamma_air/(gamma_air-1))
    Ps_lst.append(Ps)
    Psratio = Ps / P_lst[i]
    Psratio_lst.append(Psratio)
    Psratio_lst.append(Psratio)
    delta_s = cp_a * np.log(Ts_ratio_lst[i]) - R * np.log(Psratio_lst[i]) #Enthalpy
    delta_s_lst.append(delta_s)
    delta_s_lst.append(delta_s)
    h = cp_a * deltaT
    print(Tsratio)

 #   print(Ts_ratio_lst)

mean_radius = tip_speed / ((2*np.pi*RPM)/60)
blade_height = mdot_engine / (rho_amb * V1 * 2 * np.pi * mean_radius)




#plt.plot(x_axis[1:], T_ratio_lst[1:])
#plt.plot(x_axis[1:], beta_lst[1:])
#plt.plot(x_axis[1:], Psratio_lst[1:])
plt.plot(x_axis[1:], delta_s_lst[1:])
plt.show()