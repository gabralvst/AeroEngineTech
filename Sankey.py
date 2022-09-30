from matplotlib import pyplot as plt
import numpy as np
from matplotlib.sankey import Sankey

comb_eff = 0.995
thdy_eff = 0.62935
gg_eff = 0.88987
prop_eff = 0.77297

loss_comb = 100 * (1 - comb_eff)
loss_thdy = (100 - loss_comb) * (1 - thdy_eff)
loss_gg = (100 - loss_comb - loss_thdy) * (1 - gg_eff)
loss_prop = (100 - loss_comb - loss_thdy - loss_gg) * (1 - prop_eff)
rest = 100 - loss_comb - loss_thdy - loss_gg - loss_prop
losses = np.array([loss_comb, loss_thdy, loss_gg, loss_prop])
print(losses)

Sankey(flows=[1, -0.005, -0.3688, -0.06896, -0.1265, -0.43],
       labels=['', 'Combustion Loss', 'Heat Loss', 'Heat Loss', 'Kinetic Energy', 'Thrust Power'],
       orientations=[0, 1, 1, 1, 1, 0]).finish()

plt.tight_layout()
plt.show()



