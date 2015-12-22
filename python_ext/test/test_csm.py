import numpy as np
import py_cosmo_mad as csm
import matplotlib.pyplot as plt

r_arr=np.linspace(0,200,1000)

csm.background_set(0.315,0.685,0.049,-1,0,0.67,2.725);
csm.set_linear_pk("EH",-3,3,0.01,0.95,0.8)
xi_arr=np.array([csm.xi2p_L(r,2,2,"Gauss","Gauss") for r in r_arr])

plt.xlabel("$r\\,({\\rm Mpc}/h)$",fontsize=20)
plt.ylabel("$r^2\\,\\xi(r)$",fontsize=20)
plt.plot(r_arr,r_arr*r_arr*xi_arr)
plt.show()
