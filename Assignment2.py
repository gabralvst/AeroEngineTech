'Assignment 2'

gross_th= 15167 #N
mdot_a = 23.81  #kg/s
mdot_f = 0.4267 #kg/s
tot_comp_ratio = 5.5
t0 = 288 #k
p0 = 100000 #Pa
comb_eff = 1
nozz_eff = 1
gamma_air = 1.4
gamma_gas = 1.33
LHV = 43 * 10**6 #J
cp_a = 1000
cp_g = 1150
p_t2 = p0 * tot_comp_ratio


def isentropic_flow(t0, p2, p0, gamma_air):
    isentropic_flow = t0 * (p2 / p0)**((gamma_air-1)/gamma_air)
    return isentropic_flow

def T_current(T_previous, isentropic_eff, p_current, p_previous, gamma): #in K
    T_current = T_previous * (1 + 1/isentropic_eff * ((p_current/p_previous)**((gamma-1)/gamma) -1))
    return T_current

t_t2isentropic = isentropic_flow(t0, p_t2, p0, gamma_air)
t_t2 = T_current(t0, 0.92, p_t2, p0, gamma_air)
print(t_t2isentropic)
print(t_t2)