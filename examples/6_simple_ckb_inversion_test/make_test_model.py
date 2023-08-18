# %% [markdown]
# # notebook for create init and true test model

# %%
import numpy as np
import math

# grid
R_earth = 6371.0

rr1=6371 - 50
rr2=6371 + 10
tt1=(30.0)
tt2=(32.0)
pp1=(30.0)
pp2=(32.0)



# %%
# build the models
import h5py
import os

try:
    os.mkdir("models")
except:
    print("dir models exists")

try:
    os.mkdir("OUTPUT_FILES")
except:
    print("dir models exists")

vel_pert = 0.06
ani_pert = 0.04

Ngrid = 61
Nt_ckb = 2
Np_ckb = 2
Nr_ckb = 2

# model
n_rtp = [Ngrid,Ngrid,Ngrid]
dr = (rr2-rr1)/(n_rtp[0]-1)
dt = (tt2-tt1)/(n_rtp[1]-1)
dp = (pp2-pp1)/(n_rtp[2]-1)
rr = np.array([rr1 + x*dr for x in range(n_rtp[0])])
tt = np.array([tt1 + x*dt for x in range(n_rtp[1])])
pp = np.array([pp1 + x*dp for x in range(n_rtp[2])])

eta_init    = np.zeros(n_rtp)
xi_init     = np.zeros(n_rtp)
zeta_init   = np.zeros(n_rtp)
vel_init    = np.zeros(n_rtp)

eta_ckb     = np.zeros(n_rtp)
xi_ckb      = np.zeros(n_rtp)
zeta_ckb    = np.zeros(n_rtp)
vel_ckb     = np.zeros(n_rtp)

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


            if (tt[it] >= 30.5 and tt[it] <= 31.5 and pp[ip] >= 30.5 and pp[ip] <= 31.5 and dep >= 0  and dep <= 40):
                sigma = math.sin(math.pi*(tt[it]-30.5)/(1/Nt_ckb)) \
                      * math.sin(math.pi*(pp[ip]-30.5)/(1/Np_ckb))  \
                      * math.sin(math.pi*(dep)/(40/Nr_ckb))
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

            # print(vel_init[ir,it,ip],vel_ckb[ir,it,ip])

# write out in hdf5 format

fout_init = h5py.File('models/model_init_N%d_%d_%d.h5'%(n_rtp[0],n_rtp[1],n_rtp[2]), 'w')
# write out the arrays eta_init, xi_init, zeta_init, fun_init, a_init, b_init, c_init, f_init
fout_init.create_dataset('eta', data=eta_init)
fout_init.create_dataset('xi', data=xi_init)
fout_init.create_dataset('zeta', data=zeta_init)
fout_init.create_dataset('vel', data=vel_init)

fout_init.close()

fout_ckb = h5py.File('models/model_ckb_N%d_%d_%d.h5'%(n_rtp[0],n_rtp[1],n_rtp[2]), 'w')
# write out the arrays eta_init, xi_init, zeta_init, fun_init, a_init, b_init, c_init, f_init
fout_ckb.create_dataset('eta', data=eta_ckb)
fout_ckb.create_dataset('xi', data=xi_ckb)
fout_ckb.create_dataset('zeta', data=zeta_ckb)
fout_ckb.create_dataset('vel', data=vel_ckb)

fout_ckb.close()



# %%
class SRC:
    def __init__(self):
        self.id = 0
        self.year    = 1998
        self.month   = 1
        self.day     = 1
        self.hour    = 0
        self.minute  = 0
        self.sec     = 0
        self.lat     = 0 / math.pi * 180
        self.lon     = 0 / math.pi * 180
        self.dep     = 0
        self.mag     = 3.0
        self.name    = "s0"
    
class REC:
    def __init__(self):
        self.id = 0
        self.lat     = 0 / math.pi * 180
        self.lon     = 0 / math.pi * 180
        self.ele     = 0
        self.name    = "r0"

def cal_azimuth(lat1, lon1, lat2, lon2):
    lat1_rad = lat1 * math.pi / 180
    lon1_rad = lon1 * math.pi / 180
    lat2_rad = lat2 * math.pi / 180
    lon2_rad = lon2 * math.pi / 180

    y = math.sin(lon2_rad - lon1_rad) * math.cos(lat2_rad)
    x = math.cos(lat1_rad) * math.sin(lat2_rad) - math.sin(lat1_rad) * math.cos(lat2_rad) * math.cos(lon2_rad - lon1_rad)
    brng = math.atan2(y, x) * 180 / math.pi
    if((lat1-lat2)**2+(lon1-lon2)**2<0.0001):
        return 0
    return float((brng + 360.0) % 360.0)

# %%
# generate src:

import random
seed_value = 42
random.seed(seed_value)

srcs = []
Nsrcs = 600

for i in range(Nsrcs):
    src = SRC()
    src.lat     = random.uniform(30.2, 31.8)
    src.lon     = random.uniform(30.2, 31.8)
    src.dep     = random.uniform(2, 35)
    src.id      = i
    src.name    = "s%d"%(i)
    srcs.append(src)

# generate rec:
recs = []
Nrecs = 25
count = 0
for i in range(5):
    for j in range(5):
        rec = REC()
        rec.lat     = 30.2 + i/4*1.6 + random.uniform(-0.1, 0.1)
        rec.lon     = 30.2 + j/4*1.6 + random.uniform(-0.1, 0.1)
        rec.ele     = 300
        rec.id      = count
        rec.name    = "r%d"%(count)
        recs.append(rec)
        count = count + 1




# %%
# -------------------- plot stations and receivers -------------------------
import pygmt
from pygmt.clib import Session
vel_per     = 4     

# 开始画图
with pygmt.clib.Session() as session:
    session.call_module('gmtset', 'FONT 20p')
fig = pygmt.Figure()
pygmt.config(IO_SEGMENT_MARKER="<<<")


projection  = "M11c"
frame       =   ["xa2","ya2","NsWe"]
stlat=30.0
edlat=32.0
stlon=30.0
edlon=32.0

# ----------- base map ----------------
fig.shift_origin(xshift=5,yshift = 5)
fig.basemap(
    frame=frame,            
    projection=projection,  
    region=[stlon,edlon,stlat,edlat],      
)

pygmt.makecpt(cmap="jet", series=[0, 50], background = True)
x = []
y = []
z = []
for src in srcs:
    x.append(src.lon)
    y.append(src.lat)
    z.append(src.dep)
fig.plot(x = x, y = y, size = [0.4]*len(x), style = "a", fill = z, cmap = True, pen = "0.5p,black")
fig.colorbar(frame=["x+lDepth", "y+lkm"])  

x = []
y = []
for rec in recs:
    x.append(rec.lon)
    y.append(rec.lat)
fig.plot(x = x, y = y, size = [1]*len(x), style = "t", fill = 'black', pen = "1p,white")


fig.savefig('img/src_rec.jpg')

# %%
# generate src_rec_file
fname = 'src_rec_config.dat'
doc_src_rec = open(fname,'w')

accept_thd_abs = 0.3
accept_thd_dd = 0.9

azimuth_thd = 15


# loop all sources
for isrc in range(len(srcs)):
    src = srcs[isrc]
    # data info
    data_info = []

    # case 1. abs
    phase = "P"
    for rec in recs:
        accept = random.uniform(0, 1)
        if (accept > accept_thd_abs):
            data_info.append('%8d %8d %6s %9.4f %9.4f %9.1f %s %8.3f\n'%(src.id,rec.id,rec.name,rec.lat,rec.lon,rec.ele,phase,0.0))
    
    # case 2. cs (common source difference)
    phase = "P,cs"
    for i in range(len(recs)-1):
        rec1 = recs[i]
        for j in range(i+1,len(recs)):
            rec2 = recs[j]

            azi1 = cal_azimuth(rec1.lat, rec1.lon, src.lat, src.lon)
            azi2 = cal_azimuth(rec2.lat, rec2.lon, src.lat, src.lon)
            azi_gap = min(abs(azi1-azi2),360 - abs(azi1-azi2))

            accept = random.uniform(0, 1)
            if (accept > accept_thd_dd and azi_gap < azimuth_thd):
                data_info.append('%8d %8d %6s %9.4f %9.4f %9.1f %8d %6s %9.4f %9.4f %9.1f %s %8.3f\n'%
                                 (src.id,rec1.id,rec1.name,rec1.lat,rec1.lon,rec1.ele,rec2.id,rec2.name,rec2.lat,rec2.lon,rec2.ele,phase,0.0))
    
    # case 3. cr (common receiver difference)
    phase = "P,cr"
    for i in range(len(recs)-1):
        rec = recs[i]
        for j in range(isrc+1,len(srcs)):
            src2 = srcs[j]

            azi1 = cal_azimuth(rec.lat, rec.lon, src.lat, src.lon)
            azi2 = cal_azimuth(rec.lat, rec.lon, src2.lat, src2.lon)
            azi_gap = min(abs(azi1-azi2),360 - abs(azi1-azi2))

            accept = random.uniform(0, 1)
            if (accept > accept_thd_dd and azi_gap < azimuth_thd):
                data_info.append('%8d %8d %6s %9.4f %9.4f %9.1f %8d %6s %9.4f %9.4f %9.4f %s %8.3f\n'%
                                 (src.id,rec.id,rec.name,rec.lat,rec.lon,rec.ele,src2.id,src2.name,src2.lat,src2.lon,src2.dep,phase,0.0))
    
    Ndata = len(data_info)

    # write event info 
    doc_src_rec.write('%8d %8d %2d %2d %2d %2d %5.2f %9.4f %9.4f %8.4f %5.2f %5d %s\n'%
                 (src.id,src.year,src.month,src.day,src.hour,src.minute,src.sec,src.lat,src.lon,src.dep,src.mag,Ndata,src.name))
    # write data info
    for info in data_info:
        doc_src_rec.write(info)

doc_src_rec.close()


