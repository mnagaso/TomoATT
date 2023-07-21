# %% [markdown]
# # notebook for create init and true test model

# %%
import numpy as np
import math

# grid
R_earth = 6371.0

rr1=6371 - 40
rr2=6371 + 10
tt1=(0.0)/180*math.pi
tt2=(1.0)/180*math.pi
pp1=(0.0)/180*math.pi
pp2=(1.0)/180*math.pi

# source
r0=(rr1+rr2)/2.0+0.01
t0=(tt1+tt2)/2.0+0.0001
p0=(pp1+pp2)/2.0+0.0001

z0 = r0 * math.sin(t0)
x0 = r0 * math.cos(t0) * math.cos(p0)
y0 = r0 * math.cos(t0) * math.sin(p0)



# %%
import math
def cal_time(src, rec, c0):
    r1  = R_earth - src.dep
    t1  = src.lat/180*math.pi
    p1  = src.lon/180*math.pi
    x1  = r1 * math.cos(t1) * math.cos(p1)
    y1  = r1 * math.cos(t1) * math.sin(p1)
    z1  = r1 * math.sin(t1)


    r2  = R_earth + rec.ele/1000
    t2  = rec.lat/180*math.pi
    p2  = rec.lon/180*math.pi
    x2  = r2 * math.cos(t2) * math.cos(p2)
    y2  = r2 * math.cos(t2) * math.sin(p2)
    z2  = r2 * math.sin(t2)


    time = math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)/c0
    return time


# %%
def linear_velocity_time(x,y,z):
    c0=7.0
    gx=-0.002*x0/math.sqrt(x0**2+y0**2+z0**2)
    gy=-0.002*y0/math.sqrt(x0**2+y0**2+z0**2)
    gz=-0.002*z0/math.sqrt(x0**2+y0**2+z0**2)
    vel = c0 + gx*(x-x0) +gy*(y-y0) + gz*(z-z0)

    normg = math.sqrt(gx**2+gy**2+gz**2)
    dis = math.sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
    arcz = 1+0.5*(1.0/vel)*(1.0/c0)*normg**2*dis**2
    arccosh = math.log(arcz+math.sqrt(arcz**2-1))

    time = arccosh/normg

    return [vel,time]

constant_v = 5.0

# %%
import h5py
import os

try:
    os.mkdir("models")
except:
    print("dir models exists")

for Ngrid in [61]:
    # model
    n_rtp = [Ngrid,Ngrid,Ngrid]
    dr = (rr2-rr1)/(n_rtp[0]-1)
    dt = (tt2-tt1)/(n_rtp[1]-1)
    dp = (pp2-pp1)/(n_rtp[2]-1)
    rr = np.array([rr1 + x*dr for x in range(n_rtp[0])])
    tt = np.array([tt1 + x*dt for x in range(n_rtp[1])])
    pp = np.array([pp1 + x*dp for x in range(n_rtp[2])])

    eta_init = np.zeros(n_rtp)
    xi_init  = np.zeros(n_rtp)
    zeta_init = np.zeros(n_rtp)
    fun_init = np.zeros(n_rtp)
    vel_init = np.zeros(n_rtp)

    for ir in range(n_rtp[0]):
        for it in range(n_rtp[1]):
            for ip in range(n_rtp[2]):

                z = rr[ir] * math.sin(tt[it])
                x = rr[ir] * math.cos(tt[it]) * math.cos(pp[ip])
                y = rr[ir] * math.cos(tt[it]) * math.sin(pp[ip])

                # [vel,time] = linear_velocity_time(x,y,z)
                vel = constant_v
                vel_init[ir,it,ip] = vel

    # write out in hdf5 format

    fout_init = h5py.File('models/model_N%d_%d_%d.h5'%(n_rtp[0],n_rtp[1],n_rtp[2]), 'w')

    # write out the arrays eta_init, xi_init, zeta_init, fun_init, a_init, b_init, c_init, f_init
    fout_init.create_dataset('eta', data=eta_init)
    fout_init.create_dataset('xi', data=xi_init)
    fout_init.create_dataset('zeta', data=zeta_init)
    fout_init.create_dataset('vel', data=vel_init)

    fout_init.close()


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

# %%
# generate src:
srcs = []
Nsrcs = 4
for i in range(Nsrcs):
    src = SRC()
    src.lat     = i/(Nsrcs-1)*0.8 + 0.1
    src.lon     = 0.5
    src.dep     = 30
    src.id      = i
    src.name    = "s%d"%(i)
    srcs.append(src)

# generate rec:
recs = []
Nrecs = 6
for i in range(Nrecs):
    rec = REC()
    rec.lat     = i/(Nrecs-1)*0.7 + 0.15
    rec.lon     = 0.5
    rec.ele     = 300
    rec.id      = i
    rec.name    = "r%d"%(i)
    recs.append(rec)


fname = 'src_rec_obs.dat'
doc_src_rec = open(fname,'w')

# event 1
#   earthquake:         s0
#   data1 abs:          r0    
#   data2 cs_dif:       r1  r2              
#   data3 cr_dif:       r3  s1

src = srcs[0]
Ndata = 3
doc_src_rec.write('%8d %8d %2d %2d %2d %2d %5.2f %9.4f %9.4f %8.4f %5.2f %5d %s\n'%
                 (src.id,src.year,src.month,src.day,src.hour,src.minute,src.sec,src.lat,src.lon,src.dep,src.mag,Ndata,src.name))

# data 1
rec = recs[0]
phase = "P"
time  = cal_time(src,rec,constant_v)
noise = 1.0
doc_src_rec.write('%8d %8d %6s %9.4f %9.4f %9.1f %s %8.3f\n'%(src.id,rec.id,rec.name,rec.lat,rec.lon,rec.ele,phase,time+noise))

# data 2
rec1 = recs[1]
rec2 = recs[2]
time1 = cal_time(src,rec1,constant_v)
time2 = cal_time(src,rec2,constant_v)
phase = "P,cs"
noise = 1.0
doc_src_rec.write('%8d %8d %6s %9.4f %9.4f %9.1f %8d %6s %9.4f %9.4f %9.1f %s %8.3f\n'%(src.id,rec1.id,rec1.name,rec1.lat,rec1.lon,rec1.ele,rec2.id,rec2.name,rec2.lat,rec2.lon,rec2.ele,phase,time1-time2+noise))

# data 3
rec = recs[3]
src2 = srcs[1]
time1 = cal_time(src,rec,constant_v)
time2 = cal_time(src2,rec,constant_v)
phase = "P,cr"
noise = 1.0
doc_src_rec.write('%8d %8d %6s %9.4f %9.4f %9.1f %8d %6s %9.4f %9.4f %9.4f %s %8.3f\n'%(src.id,rec.id,rec.name,rec.lat,rec.lon,rec.ele,src2.id,src2.name,src2.lat,src2.lon,src2.dep,phase,time1-time2+noise))


# event 2
#   earthquake:         s2  
#   data1 cs_dif:       r2  r4  
#   data2 cr_dif:       r3  s1            
#   data3 cr_dif:       r4  s1
#   data4 cr_dif:       r2  s3

src = srcs[2]
Ndata = 4
doc_src_rec.write('%8d %8d %2d %2d %2d %2d %5.2f %9.4f %9.4f %8.4f %5.2f %5d %s\n'%
                 (src.id,src.year,src.month,src.day,src.hour,src.minute,src.sec,src.lat,src.lon,src.dep,src.mag,Ndata,src.name))

# data 1
rec1 = recs[2]
rec2 = recs[4]
time1 = cal_time(src,rec1,constant_v)
time2 = cal_time(src,rec2,constant_v)
phase = "P,cs"
noise = 1.0
doc_src_rec.write('%8d %8d %6s %9.4f %9.4f %9.1f %8d %6s %9.4f %9.4f %9.1f %s %8.3f\n'%(src.id,rec1.id,rec1.name,rec1.lat,rec1.lon,rec1.ele,rec2.id,rec2.name,rec2.lat,rec2.lon,rec2.ele,phase,time1-time2+noise))

# data 2
rec = recs[3]
src2 = srcs[1]
time1 = cal_time(src,rec,constant_v)
time2 = cal_time(src2,rec,constant_v)
phase = "P,cr"
noise = 1.0
doc_src_rec.write('%8d %8d %6s %9.4f %9.4f %9.1f %8d %6s %9.4f %9.4f %9.4f %s %8.3f\n'%(src.id,rec.id,rec.name,rec.lat,rec.lon,rec.ele,src2.id,src2.name,src2.lat,src2.lon,src2.dep,phase,time1-time2+noise))

# data 3
rec = recs[4]
src2 = srcs[1]
time1 = cal_time(src,rec,constant_v)
time2 = cal_time(src2,rec,constant_v)
phase = "P,cr"
noise = 1.0
doc_src_rec.write('%8d %8d %6s %9.4f %9.4f %9.1f %8d %6s %9.4f %9.4f %9.4f %s %8.3f\n'%(src.id,rec.id,rec.name,rec.lat,rec.lon,rec.ele,src2.id,src2.name,src2.lat,src2.lon,src2.dep,phase,time1-time2+noise))

# data 4
rec = recs[2]
src2 = srcs[3]
time1 = cal_time(src,rec,constant_v)
time2 = cal_time(src2,rec,constant_v)
phase = "P,cr"
noise = 1.0
doc_src_rec.write('%8d %8d %6s %9.4f %9.4f %9.1f %8d %6s %9.4f %9.4f %9.4f %s %8.3f\n'%(src.id,rec.id,rec.name,rec.lat,rec.lon,rec.ele,src2.id,src2.name,src2.lat,src2.lon,src2.dep,phase,time1-time2+noise))


# event 3
#   earthquake:         s3  
#   data1 abs:          r0    
#   data2 abs:          r2              
#   data3 cs_dif:       r1  r2
#   data4 cr_dif:       r3  s0

src = srcs[3]
Ndata = 4
doc_src_rec.write('%8d %8d %2d %2d %2d %2d %5.2f %9.4f %9.4f %8.4f %5.2f %5d %s\n'%
                 (src.id,src.year,src.month,src.day,src.hour,src.minute,src.sec,src.lat,src.lon,src.dep,src.mag,Ndata,src.name))

# data 1
rec = recs[0]
phase = "P"
time  = cal_time(src,rec,constant_v)
noise = 1.0
doc_src_rec.write('%8d %8d %6s %9.4f %9.4f %9.1f %s %8.3f\n'%(src.id,rec.id,rec.name,rec.lat,rec.lon,rec.ele,phase,time+noise))

# data 2
rec = recs[2]
phase = "P"
time  = cal_time(src,rec,constant_v)
noise = 1.0
doc_src_rec.write('%8d %8d %6s %9.4f %9.4f %9.1f %s %8.3f\n'%(src.id,rec.id,rec.name,rec.lat,rec.lon,rec.ele,phase,time+noise))

# data 3
rec1 = recs[1]
rec2 = recs[2]
time1 = cal_time(src,rec1,constant_v)
time2 = cal_time(src,rec2,constant_v)
phase = "P,cs"
noise = 1.0
doc_src_rec.write('%8d %8d %6s %9.4f %9.4f %9.1f %8d %6s %9.4f %9.4f %9.1f %s %8.3f\n'%(src.id,rec1.id,rec1.name,rec1.lat,rec1.lon,rec1.ele,rec2.id,rec2.name,rec2.lat,rec2.lon,rec2.ele,phase,time1-time2+noise))

# data 4
rec = recs[3]
src2 = srcs[0]
time1 = cal_time(src,rec,constant_v)
time2 = cal_time(src2,rec,constant_v)
phase = "P,cr"
noise = 1.0
doc_src_rec.write('%8d %8d %6s %9.4f %9.4f %9.1f %8d %6s %9.4f %9.4f %9.4f %s %8.3f\n'%(src.id,rec.id,rec.name,rec.lat,rec.lon,rec.ele,src2.id,src2.name,src2.lat,src2.lon,src2.dep,phase,time1-time2+noise))


doc_src_rec.close()

        


