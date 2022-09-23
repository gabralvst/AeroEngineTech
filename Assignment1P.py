from matplotlib import pyplot as plt
import numpy as np
import math


#Engine parameters
inlet_p_ratio = 0.98
mdot_a = 173 #kg/s
BPR = 12 #12 for Part 1, Part 2 (9,11,13)
FPR = 1.4 #1.4 for Part 1, Part 2 (1.4, 1.5)
LPCPR = 1.7
HPCPR = 15 #12.5 for Part 1, Part 2 (11,13,15)
T4 = 1400 #K (1400,1500,1600)
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
mdot_core = mdot_a / (BPR + 1) #Change index for BPR for Part 2
mdot_bypass = mdot_a * BPR / (BPR + 1)
print('mdot_core=',mdot_core)
print('mdot_BP',mdot_bypass)

#Station 21
mdot_21 = mdot_core
total_p21 = total_p2 * FPR
total_t21 = T_current(total_t2, fan_isentropic_eff, total_p21, total_p2, gamma_a)
W_fan = work(mdot_a, cp_a, total_t21, total_t2)
print('pt_21=',total_p21)
print('Tt_21=',total_t21)
print('W_fan=',W_fan)

#Station 13
mdot_13 = mdot_bypass
total_t13 = total_t21
total_p13 = total_p21

#Station 25
total_p25 = total_p21 * LPCPR
total_t25 = T_current(total_t21, LPC_isentropic_eff, total_p25, total_p21, gamma_a)
mdot_25 = mdot_21
W_lpc = work(mdot_25, cp_a, total_t25, total_t21) #for LPT
print('pt_25=',total_p25)
print('Tt_25=',total_t25)
print('W_lpc=',W_lpc)

#Station 3
total_p3 = total_p25 * HPCPR
total_t3 = T_current(total_t25, HPC_isentropic_eff, total_p3, total_p25, gamma_a)
mdot_3 = mdot_25
W_hpc = work(mdot_3, cp_a, total_t3, total_t25) #for HPT
print('pt_3=',total_p3)
print('Tt_3=',total_t3)
print('W_hpc=',W_hpc)

#Station 4
total_p4 = total_p3 * comb_p_ratio
mdot_f = ff_required(mdot_3, cp_g, T4, total_t3, comb_eff, LHV)
mdot_4 = mdot_3 + mdot_f
print('pt_4=',total_p4)
print('Tt_4=',T4)
print('mdot_f=',mdot_f)
print('mdot_4=',mdot_4)

#Station 45
W_hpt = W_hpc / mech_eff #as HPT drives the HPC
total_t45 = T4 - (W_hpt / (mdot_4 * cp_g))
total_p45 = total_p_45(total_p4, HPT_isentropic_eff, total_t45, T4, gamma_g)
mdot_45 = mdot_4
print('pt_45=',total_p45)
print('Tt_45=',total_t45)

#Station 5
W_lpt = (W_lpc + W_fan) / mech_eff  #as LPT drives LPC and Fan
total_t5 = total_t45 - (W_lpt / (mdot_45 * cp_g))
total_p5 = total_p_45(total_p45, LPT_isentropic_eff, total_t5, total_t45, gamma_g)
mdot_5 = mdot_45
print('pt5=', total_p5)
print('Tt5=', total_t5)
print('mdot5=', mdot_5)

# Nozzle core, station 5,7,8
total_p7 = total_p5
total_t7 = total_t5
mdot_8 = mdot_5
chocked_core = total_p7/p_amb
chocked_limit_core = (1-(1/nozz_eff)*((gamma_g-1)/(gamma_g+1)))**(-gamma_g/(gamma_g-1))
print('chocked_limit_core', chocked_limit_core)
if chocked_core > chocked_limit_core:
    print('The core nozzle is chocked', chocked_core)
    T8 = total_t7 * 2/(gamma_g+1)
    p8 = total_p7/1.876
    v8 = np.sqrt(gamma_g * gas_const * T8)          #[m/s]
    rho8 = p8 / (gas_const * T8)                    #[kg/m3]
    A8 = mdot_8 / (rho8 * v8)                       #[m2]
    F_core = mdot_8 * (v8 - v0) + A8 * (p8-p_amb)   # [N]
    print('p8=', p8)
    print('T8=', T8)
    print('v8=', v8)
    print('rho8=',rho8)
    print('A8=',A8)
    print('F_core=', F_core)

else:
    p8 = p_amb
    T8 = T_current(total_t7,nozz_eff,p8,total_p7,gamma_g)
    v8 = np.sqrt(2*cp_g*(total_t7-T8))
    F_core = mdot_8 * (v8-v0)

# Nozzle bypass
total_p16 = total_p13
total_t16 = total_t13
mdot_18 = mdot_13
chocked_BP = total_p16/p_amb
chocked_limit_BP = (1-(1/nozz_eff)*((gamma_a-1)/(gamma_a+1)))**(-gamma_a/(gamma_a-1))
print('chocked_limit_BP', chocked_limit_BP)
if chocked_BP > chocked_limit_BP:
    print('The bypass nozzle is chocked', chocked_BP)
    T18 = total_t16 * 2/(gamma_a+1)
    p18 = total_p16/1.92
    v18 = np.sqrt(gamma_a * gas_const * T18)
    rho18 = p18 / (gas_const * T18)
    A18 = mdot_18 / (rho18 * v18)
    F_bp = mdot_18 * (v18 - v0) + A18 * (p18-p_amb)
    print('p18=', p18)
    print('T18=', T18)
    print('F_bp=', F_bp)
else:
    p18 = p_amb
    T18 = T_current(total_t16,nozz_eff,p18,total_p16,gamma_a)
    v18 = np.sqrt(2*cp_g*(total_t16-T18))
    F_core = mdot_18 * (v18-v0)
# total
F_N = F_core + F_bp
TSFC = mdot_f*10**3/ (F_N * 10**(-3)) #[g/(kN s)]
OPR = total_p3/total_p2
print('F_N=',F_N, 'N')
print('TSFC=',TSFC)
print('OPR=', OPR)


#------------- PART 2 ----------------------------
# This number come form running the program after
# manually changing  the constants the values

#Thust and SFC for FPR 1.4 and 1.5 for BPR 9, 11 and 13
FPR = np.array([1.4, 1.5])
Fn9 = np.array([21.10925036827095, 22.049937961140808])
Fn11 = np.array([18.7085814455392, 19.236106036439596])
Fn13 = np.array([16.76661704214361, 15.788601923841086])
SFC9 = np.array([15.614270308483599,14.628977695105219])
SFC11 = np.array([14.681566594592049,13.974070861020905])
SFC13 = np.array([14.041743118447886,14.593169905006821])

# Thrust and SFC varying HPCPR (11,13,15) with BPR 12 and FPR 1.4
HPCPR = np.array([11,13,15])
Fn_HPCPR = np.array([17.837785291126256,17.59537487891046,17.315404687820876])
SFC_HPCPR = np.array([14.730982270860046,14.244883239554586,13.848258769740251])
OPR = np.array([26.18,30.93, 35.7])

# Thrust and SFC varying Turbine inlet temp (1400,1500,1600) with BPR 12 and FPR 1.4 HPCPR 12.5
TIT = np.array([1400, 1500, 1600])
Fn_TIT = np.array([17.66001982889069,19.328605761285995, 20.65996997155316])
SFC_TIT = np.array([14.35687652123111,14.968069266906191,15.734829050430928])

plt.subplot(3,2,1)
plt.plot(FPR,Fn9, 'r', label='BPR=9')
plt.plot(FPR,Fn11, 'b', label='BPR=11')
plt.plot(FPR,Fn13, 'g', label='BPR=13')
plt.xlabel("FPR [-]")
plt.ylabel("Thrust [kN]")
plt.title('Thrust vs FPR for different BPR')
plt.legend()
plt.grid()

plt.subplot(3,2,2)
plt.plot(FPR,SFC9, 'r', label='BPR=9')
plt.plot(FPR,SFC11, 'b', label='BPR=11')
plt.plot(FPR,SFC13, 'g', label='BPR=13')
plt.xlabel("FPR [-]")
plt.ylabel("FSFC [g/(kN s)]")
plt.title('FSFC vs FPR for different BPR')
plt.legend()
plt.grid()

plt.subplot(3,2,3)
plt.plot(OPR,Fn_HPCPR )
plt.xlabel("OPR [-]")
plt.ylabel("Thrust [kN]")
plt.title('Thrust vs OPR for BPR=12 and FPR=1.4')
plt.grid()

plt.subplot(3,2,4)
plt.plot(OPR, SFC_HPCPR)
plt.xlabel("OPR [-]")
plt.ylabel("FSFC [g/(kN s)]")
plt.title('FSFC vs OPR for BPR=12 and FPR=1.4')
plt.grid()



plt.subplot(3,2,5)
plt.plot(TIT,Fn_TIT )
plt.xlabel("TIT [K]")
plt.ylabel("Thrust [kN]")
plt.title('Thrust vs TIT for BPR=12 and FPR=1.4')
plt.grid()

plt.subplot(3,2,6)
plt.plot(TIT,SFC_TIT )
plt.xlabel("TIT [K]")
plt.ylabel("FSFC [g/(kN s)]")
plt.title('FSFC vs TIT for BPR=12 and FPR=1.4')
plt.grid()
plt.show()