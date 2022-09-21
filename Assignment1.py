from matplotlib import pyplot as plt
import numpy as np
import math


#Engine parameters
inlet_p_ratio = 0.98
mdot_a = 173 #kg/s
BPR = np.array([12, 9, 11, 13]) #12 for Part 1, rest for Part 2
FPR = np.array([1.4, 1.5]) #1.4 for Part 1, both needed for Part 2
LPCPR = 1.7
HPCPR = np.array([12.5, 11, 13, 15]) #12.5 for Part 1, rest for Part 2
T4 = np.array([1400, 1500, 1600]) #K
fan_isentropic_eff = 0.9
LPC_isentropic_eff = 0.92
HPC_isentropic_eff = LPC_isentropic_eff
LPT_isentropic_eff = 0.9
HPT_isentropic_eff = LPT_isentropic_eff
mech_eff = 0.99
comb_eff = 0.995 #Combustor
comb_p_ratio = 0.96
nozz_eff = 0.98 #Convergent
T_amb = 218.8 #K
p_amb = 23842 #Pa
gas_const = 287 #J/kg/K
LHV = 43 * 10**6 #MJ/kg
cp_a = 1000 #J/kg/K
cp_g = 1150 #J/kg/K
gamma_a = 1.4
gamma_g = 1.33

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

#Computations
#Total upstream values
v0 = M*np.sqrt(gamma_a * gas_const * T_amb)
total_T_upstream = total_T_amb(T_amb, gamma_a, M)
total_p_upstream = total_p_amb(p_amb, gamma_a, M)

#Station 2
total_t2 = total_T_upstream
total_p2 = total_p_upstream * inlet_p_ratio

#Mass Flows
mdot_core = mdot_a / (BPR[0] + 1) #Change index for BPR for Part 2
mdot_bypass = mdot_a * BPR[0] / (BPR[0] + 1)

#Station 21
mdot_21 = mdot_core
total_p21 = total_p2 * FPR[0]
total_t21 = T_current(total_t2, fan_isentropic_eff, total_p21, total_p2, gamma_a)
W_fan = work(mdot_21, cp_a, total_t21, total_t2)

#Station 13
mdot_13 = mdot_bypass
total_t13 = total_t21
total_p13 = total_p21

#Station 25
total_p25 = total_p21 * LPCPR
total_t25 = T_current(total_t21, LPC_isentropic_eff, total_p25, total_p21, gamma_a)
mdot_25 = mdot_21
W_lpc = work(mdot_25, cp_a, total_t25, total_t21) #for LPT

#Station 3
total_p3 = total_p25 * HPCPR[0]
total_t3 = T_current(total_t25, HPC_isentropic_eff, total_p3, total_p25, gamma_a)
mdot_3 = mdot_25
W_hpc = work(mdot_3, cp_a, total_t3, total_t25) #for HPT

#Station 4
total_p4 = total_p3 * comb_p_ratio
mdot_f = ff_required(mdot_3, cp_g, T4[0], total_t3, comb_eff, LHV)
mdot_4 = mdot_3 + mdot_f

#Station 45
W_hpt = W_hpc / mech_eff #as HPT drives the HPC
total_t45 = T4[0] - (W_hpt / (mdot_4 * cp_g))
total_p45 = total_p_45(total_p4, HPT_isentropic_eff, total_t45, T4[0], gamma_g)
mdot_45 = mdot_4

#Station 5
W_lpt = (W_lpc + W_fan) / mech_eff  #as LPT drives LPC and Fan
total_t5 = total_t45 - (W_lpt / (mdot_45 * cp_g))
total_p5 = total_p_45(total_p45, LPT_isentropic_eff, total_t5, total_t45, gamma_g)
mdot_5 = mdot_45

