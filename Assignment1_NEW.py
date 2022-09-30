from matplotlib import pyplot as plt
import numpy as np
import math
import pyromat as pm


#Engine parameters
inlet_p_ratio = 0.98
mdot_a = 173 #kg/s
BPR = 9# 12 for Part 1, Part 2 (9,11,13)
FPR = 1.4 # 1.4 for Part 1, Part 2 (1.4, 1.5)
LPCPR = 1.7
HPCPR = 12.5# 12.5 for Part 1, Part 2 (11,13,15)
T4 = 1400 #K 1400 for part 1, Part 2  (1400,1500,1600)
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
print('chocked_limit_core = ', chocked_limit_core)
if chocked_core > chocked_limit_core:
    print('The core nozzle is chocked', chocked_core)
    T8 = total_t7 * 2/(gamma_g+1)
    p8 = total_p7/chocked_limit_core
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
    print('The core nozzle is NOT chocked', chocked_core)
    p8 = p_amb
    T8 = total_t7-(total_t7*nozz_eff*(1-(p_amb/total_p7)**((gamma_g-1)/gamma_g)))
    v8 = np.sqrt(2*cp_g*(total_t7-T8))
    F_core = mdot_8 * (v8-v0)
    print('p8=', p8)
    print('T8=', T8)
    print('v8=', v8)
    print('F_core=', F_core)

# Nozzle bypass
total_p16 = total_p13
total_t16 = total_t13
mdot_18 = mdot_13
chocked_BP = total_p16/p_amb
chocked_limit_BP = (1-(1/nozz_eff)*((gamma_a-1)/(gamma_a+1)))**(-gamma_a/(gamma_a-1))
print('chocked_limit_BP=', chocked_limit_BP)
if chocked_BP > chocked_limit_BP:
    print('The bypass nozzle is chocked', chocked_BP)
    T18 = total_t16 * 2/(gamma_a+1)
    p18 = total_p16/chocked_limit_BP
    v18 = np.sqrt(gamma_a * gas_const * T18)
    rho18 = p18 / (gas_const * T18)
    A18 = mdot_18 / (rho18 * v18)
    F_bp = mdot_18 * (v18 - v0) + A18 * (p18-p_amb)
    print('p18=', p18)
    print('T18=', T18)
    print('F_bp=', F_bp)
else:
    print('The bypass nozzle is NOT chocked', chocked_BP)
    p18 = p_amb
    T18 = total_t16 - (total_t16 * nozz_eff * (1 - (p_amb / total_p16) ** ((gamma_a - 1) / gamma_a)))
    v18 = np.sqrt(2*cp_a*(total_t16-T18))
    F_bp = mdot_18 * (v18-v0)
    print('p18=', p18)
    print('T18=', T18)
    print('F_bp=', F_bp)

# total
M8 = v8/np.sqrt(gamma_g*gas_const*T8)
M18 = v18/np.sqrt(gamma_a*gas_const*T18)
total_t8 = total_T_amb(T8,gamma_g,M8)
total_t18 = total_T_amb(T18,gamma_a,M18)
total_p8 = total_p_amb(p8,gamma_g,M8)
total_p18 = total_p_amb(p18,gamma_a,M18)
print('Tt8',total_t8)
print('pt8', total_p8)
F_N = (F_core + F_bp)*10**(-3) # kN
TSFC = mdot_f*10**3 /F_N #[g/(kN s)]
OPR = inlet_p_ratio*FPR*LPCPR*HPCPR
v_j_eff_core = (F_core / mdot_8) +  v0 #If choked, v_j becomes v_j_eff
v_j_eff_bp = (F_bp / mdot_18) + v0
K_e_core = 0.5 * mdot_8 * (v_j_eff_core**2 - v0**2)
K_e_bp = 0.5 * mdot_18 * (v_j_eff_bp**2 - v0**2)
print('F_N=',F_N, 'kN')
print('TSFC=',TSFC)
print('OPR=', OPR)


#------------------T-S ---------------
pm.config['unit_pressure']='Pa'
pm.config['unit_temperature']='K'
pm.config['unit_energy']='J'
air = pm.get('ig.air')

s1 = air.s(T_amb,p_amb)
s2 = air.s(total_t2,total_p2)
s21 = air.s (total_t21, total_p21)
s25 = air.s (total_t25, total_p25)
s3 = air.s (total_t3, total_p3)
s4 = air.s (T4, total_p4)
s45 = air.s (total_t45, total_p45)
s5 = air.s (total_t5, total_p5)
s8 = air.s (T8,p8)
s18 = air.s (T18,p18)


T34 = np.linspace(total_t3,T4)     # This plots the curved line between 3 and 4
p34 = np.linspace(total_p3,total_p4)
s34 = np.zeros(np.shape(T34))
for i in range(len(T34)):
    s34[i] = air.s(T34[i],p34[i])
    i = i+1
plt.plot(s34, T34,'b',linewidth=1.5)

T08 = np.linspace(T_amb,T8)    # This plots the curved line between 1 and 8
p08 = np.linspace(p_amb,p8)
s08 = np.zeros(np.shape(T08))
for i in range(len(T08)):
    s08[i] = air.s(T08[i],p08[i])
    i = i+1
plt.plot(s08, T08,'b',linewidth=1.5)


plt.plot([s1,s2,s21,s25,s3],[T_amb,total_t2,total_t21,total_t25,total_t3],'b-o',linewidth=1.5)
plt.plot ([s4,s45,s5,s8], [T4,total_t45,total_t5,T8,],'b-o',linewidth=1.5)
plt.xlabel("s[J/(K kg)]")
plt.ylabel("T [K]")
plt.title('T-S Diagram')
plt.grid()
plt.show()

plt.plot ([s1,s2,s21,s18,s1], [T_amb,total_t2,total_t21,T18,T_amb],'g-o',linewidth=1.5)
plt.xlabel("s[J/(K kg)]")
plt.ylabel("T [K]")
plt.title('T-S Diagram')
plt.grid()
plt.show()



#------------- PART 2 ----------------------------
# This number come from running the program after
# manually changing  the values
# Thrust and SFC for BPR 9, 11 and 13
BPR = np.array ([9,11,13])
F_BPR_1_4 = np.array([21.105759284706633,18.705033402814526,16.64172821384842])
SFC_BPR_1_4 = np.array([15.616853050080687,14.684351450647688,14.147120205656035])
F_BPR_1_5 = np.array([22.046691016805507, 19.103040824486065,15.73834106092552])
SFC_BPR_1_5 = np.array ([14.631132189683212,14.071409432301765,14.63977362958362])

#Thust and SFC for FPR 1.4 and 1.5 for BPR 9, 11 and 13
FPR = np.array([1.4, 1.5])
Fn9 = np.array([21.105759284706633, 22.046691016805507])
Fn11 = np.array([18.705033402814526, 19.103040824486065])
Fn13 = np.array([16.64172821384842, 15.73834106092552]) # for BPR 13 and FPR 1.5 the BP thrust is negative

SFC9 = np.array([15.616853050080687,14.631132189683212])
SFC11 = np.array([14.684351450647688,14.071409432301765])
SFC13 = np.array([14.147120205656035,14.63977362958362]) #

# Thrust and SFC varying HPCPR (11,13,15) with BPR 12 and FPR 1.4
HPCPR = np.array([11,13,15])
Fn_HPCPR = np.array([17.834219852263917,17.591809404438415,17.311841253316317])
SFC_HPCPR = np.array([14.733927306701482,14.247770365397088,13.851109267408646])
OPR = np.array([25.65,30.32, 34.99])

# Thrust and SFC varying Turbine inlet temp (1400,1500,1600) with BPR 12 and FPR 1.4 HPCPR 12.5
TIT = np.array([1400, 1500, 1600])
Fn_TIT = np.array([17.656454142426373,19.325009704876194, 20.65635392983446])
SFC_TIT = np.array([14.359775864432676,14.970854570626686,15.73758354420443])

plt.subplot(4,2,1)
plt.plot(BPR,F_BPR_1_4, label = 'FPR=1.4')
plt.plot(BPR,F_BPR_1_5, label = 'FPR=1.5')
plt.xlabel("BPR [-]")
plt.ylabel("Thrust [N]")
plt.title('Thrust vs BPR')
plt.grid()
plt.legend()

plt.subplot(4,2,2)
plt.plot(BPR,SFC_BPR_1_4, label = 'FPR=1.4')
plt.plot(BPR,SFC_BPR_1_5, label = 'FPR=1.5')
plt.xlabel("BPR [-]")
plt.ylabel("FSFC [g/(kN s)]")
plt.title('FSFC vs BPR')
plt.legend()
plt.grid()

plt.subplot(4,2,3)
plt.plot(FPR,Fn9, 'r', label='BPR=9')
plt.plot(FPR,Fn11, 'b', label='BPR=11')
plt.plot(FPR,Fn13, 'g', label='BPR=13')
plt.xlabel("FPR [-]")
plt.ylabel("Thrust [kN]")
plt.title('Thrust vs FPR for different BPR')
plt.legend(loc='upper left', ncol=2)
plt.grid()

plt.subplot(4,2,4)
plt.plot(FPR,SFC9, 'r', label='BPR=9')
plt.plot(FPR,SFC11, 'b', label='BPR=11')
plt.plot(FPR,SFC13, 'g', label='BPR=13')
plt.xlabel("FPR [-]")
plt.ylabel("TSFC [g/(kN s)]")
plt.title('TSFC vs FPR for different BPR')
plt.legend()
plt.grid()

plt.subplot(4,2,5)
plt.plot(OPR,Fn_HPCPR )
plt.xlabel("OPR [-]")
plt.ylabel("Thrust [kN]")
plt.title('Thrust vs OPR for BPR=12 and FPR=1.4')
plt.grid()

plt.subplot(4,2,6)
plt.plot(OPR, SFC_HPCPR)
plt.xlabel("OPR [-]")
plt.ylabel("TSFC [g/(kN s)]")
plt.title('TSFC vs OPR for BPR=12 and FPR=1.4')
plt.grid()

plt.subplot(4,2,7)
plt.plot(TIT,Fn_TIT )
plt.xlabel("TIT [K]")
plt.ylabel("Thrust [kN]")
plt.title('Thrust vs TIT for BPR=12 and FPR=1.4')
plt.grid()

plt.subplot(4,2,8)
plt.plot(TIT,SFC_TIT )
plt.xlabel("TIT [K]")
plt.ylabel("TSFC [g/(kN s)]")
plt.title('TSFC vs TIT for BPR=12 and FPR=1.4')
plt.grid()
plt.show()

# #--------------PART 3 -------------
#Thermodynamic Efficiency
total_tg = T4 - ((mdot_core * cp_a * (total_t3 - total_t2)) / (mdot_4 * cp_g))
total_pg = total_p_45(total_p4, LPT_isentropic_eff, total_tg, T4, gamma_g)
total_t8dash = total_tg / ((total_pg/total_p_upstream) ** (1 - 1/gamma_g))
W_gg = mdot_45 * cp_g * (total_tg - total_t8dash)
thdy_eff = W_gg / (mdot_3 * cp_a * (T4 - total_t3))
print('Ttg',total_tg)
# print(total_pg)
# print(total_t8dash)
print('The thermodynamic efficiency is=',thdy_eff)

#GG Efficiency
gg_eff = ((K_e_core + K_e_bp)) / W_gg
print('The gas generation efficiency is=',gg_eff)

#Propulsive Efficiency
W_thr_core = mdot_8 * (v_j_eff_core - v0) * v0
W_thr_bp = mdot_18 * (v_j_eff_bp - v0) * v0
prop_eff = (W_thr_core + W_thr_bp) / (K_e_core + K_e_bp)
print('The propulsive efficiency is=', prop_eff)

#Thermal Efficiency
thermal_eff = ((K_e_bp + K_e_core)) / (mdot_f * LHV) #wrong for now, maybe units??
print('The thermal efficiency is =', thermal_eff)


#Overall efficiency
total_eff = prop_eff * thermal_eff
print('The overall engine efficiency is =', total_eff)


