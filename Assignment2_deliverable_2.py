'Assignment 2'
import numpy as np
import matplotlib.pyplot as plt

gross_th= 15167 #N
mdot_a = 23.81  #kg/s
mdot_f = 0.4267 #kg/s
tot_comp_ratio = 5.5
t0 = 288 #k
p0 = 100000 #Pa
comb_eff = 1
nozz_eff = 1
cp_a = 1000 #J/kg/K
cp_g = 1150 #J/kg/K     #may be different
gamma_a = 1.4
gamma_g = 1.3
gas_const = 287 #J/kg/K

kh = 0.16095
A_nozzle = 0.076795  #m2
A_turbine = A_nozzle * 0.5103267031956857
Tt1 = 288   #K
Pt1 = 100000    #Pa
Tt3 = 900   #K
nc = 0.89  #compressor eff
nt = 0.89  #turbine eff

# kh = 1 - (A_turbine/A_nozzle) ** (2*(gamma_g-1)/ (gamma_g+1))
print('kh', kh)

Tt2 = (1 + cp_g/cp_a * kh * Tt3/Tt1) * Tt1
comp_T_ratio = (1 + cp_g/cp_a * kh * Tt3/Tt1)
Pt2 = (1 + cp_g/cp_a * kh * Tt3/Tt1) ** (gamma_a * nc/ (gamma_a-1)) * Pt1   #attention to comp eff
comp_P_ratio = (1 + cp_g/cp_a * kh * Tt3/Tt1) ** (gamma_a * nc/ (gamma_a-1))
Pt3 = Pt2 #no losses in combustion chamber

Pt4 = (A_turbine/A_nozzle)**(2*gamma_g / (2*gamma_g - nt*(gamma_g-1))) * Pt3
turb_P_ratio = (A_turbine/A_nozzle)**(2*gamma_g / (2*gamma_g - nt*(gamma_g-1)))
Tt4 = (Pt4/Pt3) ** (nt * (gamma_g-1) / gamma_g) * Tt3
turb_T_ratio = (Pt4/Pt3) ** (nt * (gamma_g-1) / gamma_g)
chocked_core = Pt4/Pt1
chocked_limit_core = (1-(1/nozz_eff)*((gamma_g-1)/(gamma_g+1)))**(-gamma_g/(gamma_g-1))

if chocked_core > chocked_limit_core:
    print('The core nozzle is chocked', chocked_core, chocked_limit_core)
else:
    print('The core nozzle is NOT chocked', chocked_core)

mdot_turbine = (1.389 * Pt3 * A_turbine) / (cp_g * Tt3)**0.5
mdot_nozzle = (1.389 * Pt4 * A_nozzle) / (cp_g * Tt4)**0.5

if chocked_core > chocked_limit_core:
    T5 = Tt4 * 2/(gamma_g+1)
    v5 = np.sqrt(gamma_g * gas_const * T5)
    v5_theirs = (cp_g * T5 * (gamma_g-1)) ** 0.5
    P4 = Pt4 / (1 + (gamma_a-1)/2 * 1)**(gamma_g/ (gamma_g-1))
    P5 = Pt4/chocked_limit_core
    F_core = mdot_turbine * (v5) + A_nozzle * (P5-Pt1)
else:
    P5 = Pt1
    T5 = Tt4 - (Tt4 * nozz_eff * (1 - (Pt1 / Pt4) ** ((gamma_g - 1) / gamma_g)))
    v5 = np.sqrt(2 * cp_g * (Tt4 - T5))
    F_core = mdot_turbine * (v5)


print('mass flow engine', mdot_turbine, mdot_nozzle)

print('Compressor T ratio', comp_T_ratio, 'Compressor P ratio', comp_P_ratio)
print('Turbine T ratio', 1/turb_T_ratio, 'Turbine P ratio', 1/turb_P_ratio)
print('T5', T5, 'v5', v5, 'Thrust', F_core)

