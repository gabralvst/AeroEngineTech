from matplotlib import pyplot as plt
import numpy as np
import math
import pyromat as pm


#Engine parameters
inlet_p_ratio = 0.98
mdot_a = 173 #kg/s
BPR = 12
FPR = 1.35
LPCPR = 1.7
HPCPR = 12.5
T4 = 1405
fan_isentropic_eff = 0.90
LPC_isentropic_eff = 0.92
HPC_isentropic_eff = LPC_isentropic_eff
LPT_isentropic_eff = 0.90
HPT_isentropic_eff = LPT_isentropic_eff
mech_eff = 0.99
comb_eff = 0.995 #Combustor
comb_p_ratio = 0.96
nozz_eff = 0.98 #Convergent
T_amb = 218.8 #K
p_amb = 23842 #Pa
rho_amb = 0.3796
gas_const = 287 #J/kg/K
LHV = 43 * 10**6 #MJ/kg
cp_a = 1000 #J/kg/K
cp_g = 1150 #J/kg/K
gamma_a = 1.4
gamma_g = 1.33
r_fan = 0.5 #meters
W_generator = 4.76 * 10 **6 #W

#Flight Conditions
M = 0.78
h = 10668 #m, altitude

#Definitions of required EQs
def total_T_amb(T, gamma, M):
    total_T_amb = T *(1+ (( gamma - 1) /2) * M **2)
    return total_T_amb

def total_p_amb(p, gamma, M): #in Pa
    total_p_amb = p * ((1+ (( gamma - 1) /2) * M **2)) ** (gamma/(gamma-1))
    return total_p_amb

def T_current(T_previous, isentropic_eff, p_current, p_previous, gamma): #in K
    T_current = T_previous * (1 + 1/isentropic_eff * ((p_current/p_previous)**((gamma-1)/gamma) -1))
    return T_current

def work(mdot_current, cp_a, t_current, t_previous): #in MW
    W = mdot_current * cp_a * (t_current - t_previous)
    return W

def ff_required(mdot_previous, cp_g, t_current, t_previous, comb_eff, LHV):
    ff_required = (mdot_previous * cp_g * (t_current - t_previous)) / (comb_eff * LHV)
    return ff_required

def total_p_45(total_p_previous, isentropic_eff, t_current, t_previous, gamma_g):
    total_p_45 = total_p_previous * (1 - (1 / isentropic_eff) * (1 - (t_current / t_previous))) ** (gamma_g / (gamma_g - 1))
    return total_p_45

#------------------ PART 1 ---------------------------------
#Computations
#Total upstream values
v0 = M*np.sqrt(gamma_a * gas_const * T_amb)
total_T_upstream = total_T_amb(T_amb, gamma_a, M)
total_p_upstream = total_p_amb(p_amb, gamma_a, M)
print('pt_a=',total_p_upstream)
print('Tt_a=',total_T_upstream)

#Station 2
total_t2 = total_T_upstream
total_p2 = total_p_upstream * inlet_p_ratio
print('pt_2=',total_p2)
print('Tt_2=',total_t2)
#Mass Flows
A_fan = 0.65 * np.pi * r_fan **2
mdot_FANS = rho_amb * A_fan * v0
mdot_core = mdot_FANS / (BPR + 1) #Change index for BPR for Part 2
mdot_bypass = mdot_FANS * BPR / (BPR + 1)
print('mdot fans=', mdot_FANS)
print('mdot_core=',mdot_core)
print('mdot_BP',mdot_bypass)

#Station 21 (FAN)
mdot_21 = mdot_core
total_p21 = total_p2 * FPR
total_t21 = T_current(total_t2, fan_isentropic_eff, total_p21, total_p2, gamma_a)


W_fan = work(mdot_FANS, cp_a, total_t21, total_t2)
number_fans = W_generator / W_fan
print('pt_21=',total_p21)
print('Tt_21=',total_t21)
print('W_fan=',W_fan)
print('Fans required=', number_fans)