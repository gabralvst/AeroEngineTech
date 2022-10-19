'Assignment 2 - Deliverable 3'
import numpy as np
gross_th= 15167 #N
mdot_a = 23.81  #kg/s
mdot_f = 0.4267 #kg/s
tot_comp_ratio = 5.5
comb_eff = 1
comp_eff = 0.8
turb_eff = 0.8
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

v0 = Mach*np.sqrt(gamma_air * gas_const * t_amb)
def total_T_amb(T, gamma, M):
    total_T_amb = T *(1+ (( gamma - 1) /2) * M **2)
    return total_T_amb

def total_p_amb(p, gamma, M): #in Pa
    total_p_amb = p * ((1+ (( gamma - 1) /2) * M **2)) ** (gamma/(gamma-1))
    return total_p_amb

def isentropic_flow(t0, p2, p0, gamma_air):
    isentropic_flow = t0 * (p2 / p0)**((gamma_air-1)/gamma_air)
    return isentropic_flow

def T_current(T_previous, isentropic_eff, p_current, p_previous, gamma): #in K
    T_current = T_previous * (1 + 1/isentropic_eff * ((p_current/p_previous)**((gamma-1)/gamma) -1))
    return T_current


t_t_amb = total_T_amb(t_amb, gamma_air, Mach)
p_t_amb = total_p_amb(p_amb, gamma_air, Mach)
p_t1 = p_t_amb #assuming no pressure change in the inlet
t_t1 = t_t_amb
print(t_t1)
print(p_t1)
t_t2 = (1 + cp_g/cp_a * kh * TIT/t_t1) * t_t1
comp_t_ratio = (1 + cp_g/cp_a * kh * (TIT/t_t1))
comp_p_ratio = (1 + cp_g/cp_a * kh * (TIT/t_t1)) ** (gamma_air * comp_eff/ (gamma_air-1))
p_t2 = p_t1 * comp_p_ratio
p_t3 = p_t2 #assuming no pressure loses in CC
turb_p_ratio = 1 / ((A_turbine/A_nozzle)**(2*gamma_gas / (2*gamma_gas - turb_eff*(gamma_gas-1)))) #assuming choked turbine and nozzle
turb_t_ratio = 1 / (turb_p_ratio)**(turb_eff * (gamma_gas-1) / gamma_gas)
t_t4 = TIT * turb_t_ratio
p_t4 = p_t3 / turb_p_ratio
choked_core = p_t4/p_t1
choked_core_lim = (1-(1/nozz_eff)*((gamma_gas-1)/(gamma_gas+1)))**(-gamma_gas/(gamma_gas-1))
mdot_engine = (1.389 * p_t3 * A_turbine) / (cp_g * TIT)**0.5
if choked_core > choked_core_lim:
    print('The core nozzle is choked', choked_core, choked_core_lim)
    t_jet = t_t4 * 2 / (gamma_gas + 1)
    v_jet = np.sqrt(gamma_gas * gas_const * t_jet) # [m/s]
    p_5 = p_t4 / choked_core_lim
    F_core_gross = mdot_engine * (v_jet) + A_nozzle * (p_5 - p_amb) #N
    print('Jet Temperature =', t_jet )
    print('Jet Velocity =', v_jet)
    print('Gross Thrust =', F_core_gross)
else:
    print('The core nozzle is NOT choked', choked_core)

print('Compressor T ratio', comp_t_ratio, 'Compressor P ratio', comp_p_ratio)
print('Turbine T ratio', turb_t_ratio, 'Turbine P ratio', turb_p_ratio)
print(mdot_engine)