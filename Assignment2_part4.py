import numpy as np
from matplotlib import pyplot as plt

gross_th= 15167 #N
mdot_a = 23.81  #kg/s   #check mdot values again
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
beta_tt = 6.989635476612321 #compressor pressure ratio
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
mdot_engine = 11.824418267533746

t_t2 = t_t1 * beta_tt**((gamma_air-1)/gamma_air)
work = cp_a*(t_t2 - t_t1)
work_per_stage = work / Nstages
tip_speed = np.sqrt(work_per_stage / work_coeff)
deltaT = work_per_stage / cp_a
t_s1 = t_t1 - ((gamma_air - 1) / 2) * (V1 ** 2 / (gamma_air * R))
h1 = cp_a * (t_t1 - t_amb)
deltah = cp_a * deltaT


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
h_lst = [h1]
rho0 = p_amb/(R*t_amb)
density = [rho0]
A_0= mdot_a /(rho0*V1)
# Area = []
# Area_list = [A_0]
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
    #delta_s_lst.append(delta_s)
    h = h_lst[i] + deltah
    h_lst.append(h)
    print(h_lst)
    # rho = Pt / (Tt*R)
    # density.append(rho)
    # A_rotor = mdot_a /(density[i]*V1)
    # Area.append(A_rotor)
    # A_stator = mdot_a / (density[i]*V2)
    # Area.append(A_stator)
    
    # A2=(density[i]*V1*Area_list[i])/(density[i+1]*V2)
    # Area_list.append(A2)
    # print(Tsratio)


    
# plt.plot([0,2,4,6,8,10,12,14,16],Area_list,'b')
# plt.plot(range(16), Area,'r')
# plt.show()
# print(Area)
#   print(Ts_ratio_lst)

chord = 0.03
rotor_space = chord * 0.1
stator_space = chord * 0.2
Acn = (mdot_a * (cp_a * t_t2) ** 0.5) / (1.281 * P_lst[-1])
Ac = 1 / ((p_t1/P_lst[-1]) ** (1 - (gamma_air - 1)/ (comp_eff * gamma_air))) * Acn
print(Ac, Acn)
Area_lst = [Ac]
for i in range(Nstages):
    A = Ac - ((i+1)*(Ac-Acn) / Nstages)
    A_rotor_exit = (Area_lst[-1] + A)/2
    Area_lst.append(A_rotor_exit)
    Area_lst.append(A)


mean_radius = tip_speed / ((2*np.pi*RPM)/60)
# blade_height = mdot_engine / (rho_amb * V1 * 2 * np.pi * mean_radius)
Blade_H = Area_lst / (2 * np.pi * mean_radius)
tip_rad = mean_radius + Blade_H
# print(tip_rad)

step = [*range(Nstages)]
rotor_x_coor1 = [i * (2*chord + rotor_space + stator_space) for i in step]
rotor_x_coor2 = [chord + i for i in rotor_x_coor1]
rotor_x_coor = rotor_x_coor1 * 2 + rotor_x_coor2 * 2
rotor_x_coor.sort()


# print(rotor_x_coor)
# xlst = [*range(17)]
# plt.plot(xlst, tip_rad)
# plt.show()


#plt.plot(x_axis[1:], T_ratio_lst[1:])
#plt.plot(x_axis[1:], beta_lst[1:])
#plt.plot(x_axis[1:], Psratio_lst[1:])
#plt.plot(x_axis[1:], delta_s_lst[1:])
#plt.show()