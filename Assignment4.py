import numpy as np
import matplotlib.pyplot as plt

Cp_air = 1000
Cp_gas = 1150
LHV = 43e6 #J/kg

# Fuel composition
C = 11
H = 22
combustor_stages = 18
O2_mole = C + H/4
air_mass = O2_mole * (79/21 * 28 + 32)
fuel_mass = C*12 + H
ratio_stoi = air_mass / fuel_mass

def m_fuel(mdot_air,Tt3, Tt4):
    m_f = mdot_air * (Cp_air * (Tt3-298) - Cp_gas * (Tt4-298)) / ((Tt4-298) * Cp_gas - LHV)
    return m_f

def eq_ratio(mdot_air_2,Tt3_2, Tt4_2):
    ratio = ratio_stoi / (mdot_air_2 / m_fuel(mdot_air_2, Tt3_2, Tt4_2))
    return ratio

def eq_ratio2(mdot_air_2,Tt3_2, Tt4_2, mdot_fuel):
    ratio = ratio_stoi / (mdot_air_2 / mdot_fuel)
    return ratio

def residence_time(mdot_gas, combustorVOL):
    time = combustorVOL / mdot_gas
    return time

def abs_heat_density(heat_release, combustorVOL):
    abs_heat_dens = heat_release / combustorVOL
    return abs_heat_dens

def norm_heat_density(heat_release, combustorVOL, air_pressure):
    norm_heat_dens = heat_release / combustorVOL / (air_pressure / 100000)
    return norm_heat_dens

mdot_air_lst = [14.5, 29.5, 37.5]
P3_lst = [1100e3, 1950e3, 3350e3]
T3_lst = [620, 765, 850]
T4_lst = [1325, 1575, 1820]
X_lst = [0.16, 0.6, 1]
mdot_fuel_lst = []
eq_ratio_lst = []
eq_ratio_RQLlst = []
T_adiabatic_lst = []
mdot_gas_lst = []
residence_lst=[] #s
combustorVOL = 0.012 * combustor_stages #m3
heat_released_lst = []
norm_heat_density_lst = []
abs_heat_density_lst = []
for i in range(0,3):
    mdot_fuel_lst.append(m_fuel(mdot_air_lst[i], T3_lst[i], T4_lst[i]))
    eq_ratio_lst.append(eq_ratio(mdot_air_lst[i], T3_lst[i], T4_lst[i]))
    for i2 in range(0,3):
        eq_ratio_RQLlst.append(eq_ratio2(mdot_air_lst[i] * X_lst[i2], T3_lst[i], T4_lst[i], mdot_fuel_lst[i]))
    mdot_gas = mdot_fuel_lst[i] + mdot_air_lst[i]
    mdot_gas_lst.append(mdot_gas)
    residence = residence_time(mdot_gas_lst[i], combustorVOL)
    residence_lst.append(residence)
    heat_released = LHV * mdot_fuel_lst[i] / 1000 #kW
    heat_released_lst.append(heat_released)
    abs_heat = abs_heat_density(heat_released_lst[i], combustorVOL)
    abs_heat_density_lst.append(abs_heat)
    norm_heat = norm_heat_density(heat_released_lst[i], combustorVOL, P3_lst[i])
    norm_heat_density_lst.append(norm_heat)
    print('heat released [kW]', heat_released_lst)
    print('mass flow of gas', mdot_gas_lst)
    print('residence times', residence_lst)
    print('absolute heat density[kW/m3]', abs_heat_density_lst)
    print('normalized heat density[kW/m3/bar]', norm_heat_density_lst)
print('fuel flow rate [kg/s]', mdot_fuel_lst)
print('overall eq. ratio', eq_ratio_lst)
id1 = 0
id2 = 3
for a in range(0,3):
    print('eq. ratio red dot', [a+1], eq_ratio_RQLlst[id1:id2])
    id1 = id1 + 3
    id2 = id2 + 3

# Enthalpy of formation Tref 1100K    [J/mol]
kCal_to_j = 4184
DeltaH_f_fuel = -66.8e3
DeltaH_f_CO2 = -394e3
DeltaH_f_H2O = -244.5e3

# Cp at 1100K [J/mol]
Cp_fuel = 0.136718 * kCal_to_j
Cp_CO2 = 55.396
Cp_H2O = 42.44
Cp_N2 = 33
Cp_O2 = 35.29

# For stoichiometric conditions
Hproducts_1 = DeltaH_f_CO2 * C + DeltaH_f_H2O * H/2
Hproducts_2 = Cp_CO2 * C + Cp_H2O * H/2 + O2_mole * 79/21 * Cp_N2

def T_adiabatic(equivalence):
    if equivalence>1:
        Hreactants = 1/ equivalence * DeltaH_f_fuel
        Hprod = 1/ equivalence * Hproducts_1
        Hprod2 = 1 / equivalence * Hproducts_2 + (1 - 1/equivalence) * (Cp_fuel)

    if equivalence<1:
        Hreactants = DeltaH_f_fuel
        Hprod = Hproducts_1
        Hprod2 = Hproducts_2 + (1/equivalence - 1) * O2_mole * (79/21 * Cp_N2 + Cp_O2)

    T = (Hreactants - Hprod) / (Hprod2) + 298
    return T

for a in eq_ratio_RQLlst:
    T_adiabatic_lst.append(T_adiabatic(a))

id1 = 0
id2 = 3
for a in range(0,3):
    print('T adiabatic red dot', [a+1], T_adiabatic_lst[id1:id2])
    id1 = id1 + 3
    id2 = id2 + 3

