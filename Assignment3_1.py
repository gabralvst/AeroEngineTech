from matplotlib import pyplot as plt
import numpy as np
import math
import pyromat as pm


#Engine parameters
inlet_p_ratio = 0.98
mdot_a = 173 #kg/s
BPR = 12
FPR = 1.3
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
Thrust_original = 17635

#new eff 2030
motor_eff = 0.95
cable_eff = 0.995
rectfier_eff = 0.99

#new eff 2040
motor_eff2 = 0.97
cable_eff2 = 0.998
rectfier_eff2 = 0.995

#Flight Conditions
M = 0.78
h = 10668 #m, altitude
M_fan = 0.55

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
v_fan = M_fan*np.sqrt(gamma_a * gas_const * T_amb)
total_T_upstream = total_T_amb(T_amb, gamma_a, M_fan)
total_p_upstream = total_p_amb(p_amb, gamma_a, M_fan)
print('pt_a=',total_p_upstream)
print('Tt_a=',total_T_upstream)

#Station 2
total_t2 = total_T_upstream
total_p2 = total_p_upstream * inlet_p_ratio
print('pt_2=',total_p2)
print('Tt_2=',total_t2)
#Mass Flows
# A_fan = 0.65 * np.pi * r_fan **2
A_fan = np.pi * (r_fan **2 - (0.35 * r_fan)**2)
v_fan = M_fan*np.sqrt(gamma_a * gas_const * T_amb)

mdot_FANS = rho_amb * A_fan * v_fan
mdot_core = 13.31

# mdot_bypass = mdot_FANS * BPR / (BPR + 1)
print('mdot fans=', mdot_FANS)
print('mdot_core=',mdot_core)
# print('mdot_BP',mdot_bypass)

#Station 21 (FAN)
# mdot_21 = mdot_core
total_p21 = total_p2 * FPR
total_t21 = T_current(total_t2, fan_isentropic_eff, total_p21, total_p2, gamma_a)

total_p16 = total_p21
total_t16 = total_t21

W_fan = work(mdot_FANS, cp_a, total_t21, total_t2)
# Thrust_fan = W_fan/v0

chocked_BP = total_p16/p_amb
chocked_limit_BP = (1-(1/nozz_eff)*((gamma_a-1)/(gamma_a+1)))**(-gamma_a/(gamma_a-1))
print('chocked_limit_BP=', chocked_limit_BP)
if chocked_BP > chocked_limit_BP:
    print('The bypass nozzle is chocked', chocked_BP)
    T18 = total_t16 * 2/(gamma_a+1)
    p18 = total_p16/chocked_limit_BP
    v18 = np.sqrt(gamma_a * gas_const * T18)
    rho18 = p18 / (gas_const * T18)
    A18 = mdot_FANS / (rho18 * v18)
    F_bp = mdot_FANS * (v18 - v_fan) + A18 * (p18-p_amb)
    print('p18=', p18)
    print('T18=', T18)
    print('F_bp=', F_bp)
else:
    print('The bypass nozzle is NOT chocked', chocked_BP)
    p18 = p_amb
    T18 = total_t16 - (total_t16 * nozz_eff * (1 - (p_amb / total_p16) ** ((gamma_a - 1) / gamma_a)))
    v18 = np.sqrt(2*cp_a*(total_t16-T18))
    F_bp = mdot_FANS * (v18-v_fan)
    print('p18=', p18)
    print('T18=', T18)
    print('F_bp=', F_bp)

# number_fans = W_generator / W_fan
number_fans_per_wing_original = Thrust_original / F_bp
Thrust_total = math.ceil(number_fans_per_wing_original) * F_bp
print('pt_21=',total_p21)
print('Tt_21=',total_t21)
print('W_fan=',W_fan)
print('Fans required=', number_fans_per_wing_original)
print('Thrust fan', F_bp)
print('Total Thrust', Thrust_total)
W_mot = W_fan/mech_eff

# Iterate FPR for optimized generator power
number_fans_per_wing = number_fans_per_wing_original
# for a in range(0,10):
while (number_fans_per_wing > math.floor(number_fans_per_wing_original)):
    print(FPR)
    FPR = FPR+0.005

    total_T_upstream = total_T_amb(T_amb, gamma_a, M_fan)
    total_p_upstream = total_p_amb(p_amb, gamma_a, M_fan)
    # print('pt_a=', total_p_upstream)
    # print('Tt_a=', total_T_upstream)

    # Station 2
    total_t2 = total_T_upstream
    total_p2 = total_p_upstream * inlet_p_ratio
    # print('pt_2=', total_p2)
    # print('Tt_2=', total_t2)
    # Mass Flows
    # A_fan = 0.65 * np.pi * r_fan **2
    A_fan = np.pi * (r_fan ** 2 - (0.35 * r_fan) ** 2)
    v_fan = M_fan * np.sqrt(gamma_a * gas_const * T_amb)

    mdot_FANS = rho_amb * A_fan * v_fan
    mdot_core = 13.31

    # mdot_bypass = mdot_FANS * BPR / (BPR + 1)
    # print('mdot fans=', mdot_FANS)
    # print('mdot_core=', mdot_core)
    # print('mdot_BP',mdot_bypass)

    # Station 21 (FAN)
    # mdot_21 = mdot_core
    total_p21 = total_p2 * FPR
    total_t21 = T_current(total_t2, fan_isentropic_eff, total_p21, total_p2, gamma_a)

    total_p16 = total_p21
    total_t16 = total_t21

    W_fan = work(mdot_FANS, cp_a, total_t21, total_t2)
    # Thrust_fan = W_fan/v0

    chocked_BP = total_p16 / p_amb
    chocked_limit_BP = (1 - (1 / nozz_eff) * ((gamma_a - 1) / (gamma_a + 1))) ** (-gamma_a / (gamma_a - 1))
    # print('chocked_limit_BP=', chocked_limit_BP)
    if chocked_BP > chocked_limit_BP:
        # print('The bypass nozzle is chocked', chocked_BP)
        T18 = total_t16 * 2 / (gamma_a + 1)
        p18 = total_p16 / chocked_limit_BP
        v18 = np.sqrt(gamma_a * gas_const * T18)
        rho18 = p18 / (gas_const * T18)
        A18 = mdot_FANS / (rho18 * v18)
        F_bp = mdot_FANS * (v18 - v_fan) + A18 * (p18 - p_amb)
        # print('p18=', p18)
        # print('T18=', T18)
        # print('F_bp=', F_bp)
    else:
        print('The bypass nozzle is NOT chocked', chocked_BP)
        p18 = p_amb
        T18 = total_t16 - (total_t16 * nozz_eff * (1 - (p_amb / total_p16) ** ((gamma_a - 1) / gamma_a)))
        v18 = np.sqrt(2 * cp_a * (total_t16 - T18))
        F_bp = mdot_FANS * (v18 - v_fan)
        # print('p18=', p18)
        # print('T18=', T18)
        # print('F_bp=', F_bp)

    # number_fans = W_generator / W_fan
    number_fans_per_wing = Thrust_original / F_bp
    Thrust_total = math.ceil(number_fans_per_wing) * F_bp
    # print('pt_21=', total_p21)
    # print('Tt_21=', total_t21)
    # print('W_fan=', W_fan)
    # print('Fans required=', number_fans_per_wing)
    # print('Thrust fan', F_bp)
    # print('Total Thrust', Thrust_total)
    W_mot = W_fan / mech_eff

print('W_fan=', W_fan)
print('Thrust fan', F_bp)
print('Fans required=', number_fans_per_wing)
Thrust_total2 = math.ceil(number_fans_per_wing) * F_bp
print('Total Thrust', Thrust_total)

W_generator2 = (math.ceil(number_fans_per_wing)*W_mot)/ rectfier_eff/cable_eff/motor_eff

powertrain_eff1 = mech_eff * motor_eff * rectfier_eff * cable_eff
powertrain_eff2 = mech_eff * motor_eff2 * rectfier_eff2 * cable_eff2
print('Required generator Power', W_generator2)
print('Motor Power', W_mot)
print('Ppwertrain eff', powertrain_eff1, powertrain_eff2)


# Part 2 Cycle calculations
F_core = 500        #Value to start iteration
# for a in range(0,30):
while F_core>5:
    T4 = T4 - 0.1
    # print(T4)

    #Total upstream values core
    v0_c = M*np.sqrt(gamma_a * gas_const * T_amb)
    total_T_upstream_c = total_T_amb(T_amb, gamma_a, M)
    total_p_upstream_c = total_p_amb(p_amb, gamma_a, M)
    # print('pt_a=',total_p_upstream_c)
    # print('Tt_a=',total_T_upstream_c)

    #Station 21 core
    mdot_21c = mdot_core
    total_p21 = total_p_upstream_c
    total_t21 = total_T_upstream_c
    # print('pt_21=',total_p21)
    # print('Tt_21=',total_t21)


    #Station 25
    total_p25 = total_p21 * LPCPR
    total_t25 = T_current(total_t21, LPC_isentropic_eff, total_p25, total_p21, gamma_a)
    mdot_25 = mdot_21c
    W_lpc = work(mdot_25, cp_a, total_t25, total_t21) #for LPT
    # print('pt_25=',total_p25)
    # print('Tt_25=',total_t25)
    # print('W_lpc=',W_lpc)

    #Station 3
    total_p3 = total_p25 * HPCPR
    total_t3 = T_current(total_t25, HPC_isentropic_eff, total_p3, total_p25, gamma_a)
    mdot_3 = mdot_25
    W_hpc = work(mdot_3, cp_a, total_t3, total_t25) #for HPT
    # print('pt_3=',total_p3)
    # print('Tt_3=',total_t3)
    # print('W_hpc=',W_hpc)

    #Station 4
    total_p4 = total_p3 * comb_p_ratio
    mdot_f = ff_required(mdot_3, cp_g, T4, total_t3, comb_eff, LHV)
    mdot_4 = mdot_3 + mdot_f
    # print('pt_4=',total_p4)
    # print('Tt_4=',T4)
    # print('mdot_f=',mdot_f)
    # print('mdot_4=',mdot_4)

    #Station 45
    W_hpt = W_hpc / mech_eff #as HPT drives the HPC
    total_t45 = T4 - (W_hpt / (mdot_4 * cp_g))
    total_p45 = total_p_45(total_p4, HPT_isentropic_eff, total_t45, T4, gamma_g)
    mdot_45 = mdot_4
    # print('pt_45=',total_p45)
    # print('Tt_45=',total_t45)

    #Station 5
    W_lpt = (W_lpc) / mech_eff + W_generator2 #as LPT drives LPC and Fan
    total_t5 = total_t45 - (W_lpt / (mdot_45 * cp_g))
    total_p5 = total_p_45(total_p45, LPT_isentropic_eff, total_t5, total_t45, gamma_g)
    mdot_5 = mdot_45
    # print('pt5=', total_p5)
    # print('Tt5=', total_t5)
    # print('mdot5=', mdot_5)

    # Nozzle core, station 5,7,8
    total_p7 = total_p5
    total_t7 = total_t5
    mdot_8 = mdot_5
    chocked_core = total_p7/p_amb
    chocked_limit_core = (1-(1/nozz_eff)*((gamma_g-1)/(gamma_g+1)))**(-gamma_g/(gamma_g-1))
    # print('chocked_limit_core = ', chocked_limit_core)
    if chocked_core > chocked_limit_core:
        # print('The core nozzle is chocked', chocked_core)
        T8 = total_t7 * 2/(gamma_g+1)
        p8 = total_p7/chocked_limit_core
        v8 = np.sqrt(gamma_g * gas_const * T8)          #[m/s]
        rho8 = p8 / (gas_const * T8)                    #[kg/m3]
        A8 = mdot_8 / (rho8 * v8)                       #[m2]
        F_core = mdot_8 * (v8 - v0) + A8 * (p8-p_amb)   # [N]
        # print('p8=', p8)
        # print('T8=', T8)
        # print('v8=', v8)
        # print('rho8=',rho8)
        # print('A8=',A8)
        # print('F_core=', F_core)

    else:
        # print('The core nozzle is NOT chocked', chocked_core)
        p8 = p_amb
        T8 = total_t7-(total_t7*nozz_eff*(1-(p_amb/total_p7)**((gamma_g-1)/gamma_g)))
        v8 = np.sqrt(2*cp_g*(total_t7-T8))
        F_core = mdot_8 * (v8-v0)
        # print('p8=', p8)
        # print('T8=', T8)
        # print('v8=', v8)
        # print('F_core=', F_core)

print('p8=', p8)
print('T8=', T8)
print('v8=', v8)
print('F_core=', F_core)
print('TIT', T4)

# Fuel consumption
mdot_fuel = mdot_3 * cp_g * (T4 - total_t3) / comb_eff / LHV
mdot_fuel_original = mdot_3 * cp_g * (1405 - total_t3) / comb_eff / LHV
mdot_fuel_reduction = mdot_fuel_original - mdot_fuel
print('mdot_fuel reduction', mdot_fuel_reduction, 'kg/s', mdot_fuel_reduction/mdot_fuel_original*100, '%')

l_cable = 35.8/2    #wingspan/2
m_motor2030 = W_generator2/ rectfier_eff / 10e3
m_cable2030 = W_generator2/rectfier_eff/1e6 * 2 * l_cable
m_motor2040 = W_generator2/rectfier_eff/20e3
m_cable2040 = W_generator2/rectfier_eff/1e6 * l_cable
print('motor mass 2030', m_motor2030, 'cable mass 2030', m_cable2030)
print('motor mass 2040', m_motor2040, 'cable mass 2040', m_cable2040)

m_PowerTrain2030 = round(number_fans_per_wing) * m_motor2030 * 2 + m_cable2030
m_PowerTrain2040 = round(number_fans_per_wing) * m_motor2040 * 2 + m_cable2040
print('Powertrain mass 2030', m_PowerTrain2030)
print('Powertrain mass 2040', m_PowerTrain2040)