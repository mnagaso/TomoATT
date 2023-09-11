# %% [markdown]
# # notebook for create init and true test model

# %%
import numpy as np
import math
import os

# study region
R_earth = 6371.0

rr1=R_earth - 50
rr2=R_earth + 10
tt1=(30.0)
tt2=(32.0)
pp1=(30.0)
pp2=(32.0)

try:
    os.mkdir("models")
except:
    print("dir models exists")


# %%
# build the initial model
import h5py



# model
n_rtp = [61,61,61]
dr = (rr2-rr1)/(n_rtp[0]-1)
dt = (tt2-tt1)/(n_rtp[1]-1)
dp = (pp2-pp1)/(n_rtp[2]-1)
rr = np.array([rr1 + x*dr for x in range(n_rtp[0])])
tt = np.array([tt1 + x*dt for x in range(n_rtp[1])])
pp = np.array([pp1 + x*dp for x in range(n_rtp[2])])

vel_init    = np.zeros(n_rtp)   # velocity
xi_init     = np.zeros(n_rtp)   # xi and eta are azimuthal anisotropy
eta_init    = np.zeros(n_rtp)
zeta_init   = np.zeros(n_rtp)


for ir in range(n_rtp[0]):
    for it in range(n_rtp[1]):
        for ip in range(n_rtp[2]):

            dep = R_earth - rr[ir]
            if (dep < 0):
                vel_init[ir,it,ip]  = 6.0
            elif (dep >= 0  and dep < 40):
                vel_init[ir,it,ip]  = 6.0 + dep/40*2.0
            else:
                vel_init[ir,it,ip]  = 8.0
            

            xi_init[ir,it,ip]   = 0.0
            eta_init[ir,it,ip]  = 0.0
            zeta_init[ir,it,ip] = 0.0


# write out in hdf5 format

fout_init = h5py.File('models/model_init_N%d_%d_%d.h5'%(n_rtp[0],n_rtp[1],n_rtp[2]), 'w')
# write out the arrays eta_init, xi_init, zeta_init, fun_init, a_init, b_init, c_init, f_init
fout_init.create_dataset('eta', data=eta_init)
fout_init.create_dataset('xi', data=xi_init)
fout_init.create_dataset('zeta', data=zeta_init)
fout_init.create_dataset('vel', data=vel_init)

fout_init.close()





