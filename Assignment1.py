import matplotlib as plt
import numpy as np
import math


#Engine parameters
inlet_p_ratio = 0.98
mdot_a = 173 #kg/s
BPR = np.array([12, 9, 11, 13]) #12 for Part 1, rest for Part 2
FPR = np.array([1.4, 1.5]) #1.4 for Part 1, both needed for Part 2
LPCPR = 1.7
HPCPR = np.array([12.5, 11, 13, 15]) #12.5 for Part 1, rest for Part 2
T4 = 1400 #K
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
LHV = 43 #MJ/kg
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

def total_p_amb(p, gamma, M):
    total_p_amb = p * ((1+ (( gamma - 1) /2) * M **2)) ** (gamma/(gamma-1))
    return total_p_amb

def T_current(T_previous, isentropic_eff, p_current, p_previous, gamma):
    T_current = T_previous * (1 + 1/isentropic_eff * ((p_current/p_previous)**((gamma-1)/gamma) -1))
    return T_current






#Computations
#Total upstream values
v0 = M*np.sqrt(gamma_a * gas_const * T_amb)
total_T_upstream = total_T_amb(T_amb, gamma_a, M)
total_p_upstream = total_p_amb(p_amb, gamma_a, M)

#Inlet
total_T2 = total_T_upstream
total_p2 = total_p_upstream * inlet_p_ratio

#Mass Flows
mdot_core = mdot_a / (BPR[0] + 1) #Change index for BPR for Part 2
mdot_21 = mdot_core
mdot_bypass = mdot_a * BPR[0] / (BPR[0] + 1)
mdot_13 = mdot_bypass

#Fan
