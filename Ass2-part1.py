import numpy as np

T0 = 288 #K
p0 = 100000 #Pa
comp_PR = 5.5
gamma_a = 1.4
gamma_g = 1.3
cp_gas =1150 # j/kg k
cp_air = 1000
eff_comp = 0.89
eff_cc = 1
eff_turb = 0.89
eff_noz = 1
LHV_f = 43*10**6 # j/kg
mdot_air = 23.81 #kg/s
mdot_f = 0.4267 
R = 287

pt2 = p0
Tt2 = T0

pt3 = comp_PR * pt2
Tt3 = Tt2 * (1 + (1/eff_comp)*((pt3/pt2)**((gamma_a-1)/gamma_a)-1))

pt4 = pt3
TIT = (mdot_air * cp_air*Tt3 + mdot_f * LHV_f*eff_cc)/(cp_gas*(mdot_air + mdot_f))

Tt5 = TIT-((cp_air*mdot_air*(Tt3-Tt2))/(cp_gas*(mdot_air+mdot_f)))
pt5 = pt4 * (1 - (1 / eff_turb) * (1 - (Tt5 / TIT))) ** (gamma_g / (gamma_g - 1))

print('pt3=',pt3)
print('Tt3=',Tt3)
print('pt4=',pt4)
print('TIT=',TIT)
print('pt5=',pt5)
print('Tt5=',Tt5)
print('pt5/pt1=',pt5/p0)

kh=1-Tt5/TIT
Area_ratio = (1-kh)**((gamma_g+1)/(2*gamma_g-2))
print('kh=',kh)
print('At/An=',Area_ratio)

choked_core_lim = (1-(1/eff_noz)*((gamma_g-1)/(gamma_g+1)))**(-gamma_g/(gamma_g-1))
T8 = Tt5 * 2 / (gamma_g + 1)
v8 = np.sqrt(gamma_g * R * T8) # [m/s]
p8 = pt5 / choked_core_lim
print('p8=',p8)
An = ((mdot_air+mdot_f)*np.sqrt(cp_gas*Tt5))/(1.389*pt5)
print('An=',An)
F = (mdot_air+mdot_f)*v8 + An*(p8-p0)
print('F=',F)







