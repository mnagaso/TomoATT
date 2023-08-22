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
# read initial model
import h5py

model_fname = "models/model_init_N61_61_61.h5"   # 检测板反演的最终模型
fmodel  = h5py.File(model_fname, 'r')

vel_init    = np.array(fmodel['vel'])
xi_init     = np.array(fmodel['xi'])
eta_init    = np.array(fmodel['eta'])

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

eta_ckb     = np.zeros(n_rtp)
xi_ckb      = np.zeros(n_rtp)
zeta_ckb    = np.zeros(n_rtp)
vel_ckb     = np.zeros(n_rtp)

# perturbation amplitude
vel_pert = 0.06
ani_pert = 0.04

for ir in range(n_rtp[0]):
    for it in range(n_rtp[1]):
        for ip in range(n_rtp[2]):

            dep = R_earth - rr[ir]

            if (tt[it] >= 30.5 and tt[it] <= 31.5 and pp[ip] >= 30.5 and pp[ip] <= 31.5 and dep >= 0  and dep <= 40):
                sigma = math.sin(math.pi*(tt[it]-30.5)/(0.5)) \
                      * math.sin(math.pi*(pp[ip]-30.5)/(0.5))  \
                      * math.sin(math.pi*(dep)/(20))
            else:
                sigma = 0.0

            if sigma < 0:
                psi = 60.0/180.0*math.pi
            elif sigma > 0:
                psi = 150.0/180.0*math.pi
            else:
                psi = 0.0

            vel_ckb[ir,it,ip]   = vel_init[ir,it,ip] * (1.0 + vel_pert * sigma)
            xi_ckb[ir,it,ip]    = ani_pert * abs(sigma) * math.cos(2*psi) 
            eta_ckb[ir,it,ip]   = ani_pert * abs(sigma) * math.sin(2*psi) 
            zeta_ckb[ir,it,ip]  = 0.0

# write out in hdf5 format
fout_ckb = h5py.File('models/model_ckb_N%d_%d_%d.h5'%(n_rtp[0],n_rtp[1],n_rtp[2]), 'w')
# write out the arrays eta_init, xi_init, zeta_init, fun_init, a_init, b_init, c_init, f_init
fout_ckb.create_dataset('eta', data=eta_ckb)
fout_ckb.create_dataset('xi', data=xi_ckb)
fout_ckb.create_dataset('zeta', data=zeta_ckb)
fout_ckb.create_dataset('vel', data=vel_ckb)

fout_ckb.close()





