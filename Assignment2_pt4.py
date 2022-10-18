import numpy as np


gross_th= 15167 #N
mdot_a = 23.81  #kg/s
mdot_f = 0.4267 #kg/s
tot_comp_ratio = 5.5
comb_eff = 1
comp_eff = 1
turb_eff = 1
nozz_eff = 1
gamma_air = 1.4
gamma_gas = 1.3
gas_const = 287 #J/kg/K
LHV = 43 * 10**6 #J
cp_a = 1000
cp_g = 1150
h = 10668 #m
Mach = 0.78
TIT = 1150 #k
t_amb = 219
p_amb = 23911 #Pa
kh = 0.152
A_turbine = 0.037   #m2
A_nozzle = 0.079  #m2
beta_tt = 8.10682594846814 #compressor pressure ratio
t_t1 = 245.64792 #K
p_t1 = 35738.733374570875 #Pa
reaction = 0.5
flow_coeff = 0.6
work_coeff = 0.4
u_max = 450 #m/s
Nstages = 10



t_t2 = t_t1 * beta_tt**((gamma_air-1)/gamma_air)
work = cp_a*(t_t2 - t_t1)
work_per_stage = work / Nstages
tip_speed = np.sqrt(work_per_stage / work_coeff)
deltaT = work_per_stage / cp_a
t_tout1 = t_t1 + deltaT
t_tratio1 = t_tout1 / t_t1
beta_stage1 = (t_tratio1)**((gamma_air*comb_eff)/(gamma_air-1))
print(work)
print(tip_speed)
print(t_tratio1)
print(beta_stage1)