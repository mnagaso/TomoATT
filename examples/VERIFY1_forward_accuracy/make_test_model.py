
import numpy as np
import math
import h5py
import os

# grid
R_earth = 6371.0

rr1=5900 
rr2=6400
tt1=(30.0)/180*math.pi
tt2=(50.0)/180*math.pi
pp1=(15.0)/180*math.pi
pp2=(40.0)/180*math.pi

# source
r0=(rr1+rr2)/2.0+0.01
t0=(tt1+tt2)/2.0+0.0001
p0=(pp1+pp2)/2.0+0.0001

z0 = r0 * math.sin(t0)
x0 = r0 * math.cos(t0) * math.cos(p0)
y0 = r0 * math.cos(t0) * math.sin(p0)


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



try:
    os.mkdir("models")
except:
    print("dir models exists")

for Ngrid in [41,61,81,121,161]:
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

                [vel,time] = linear_velocity_time(x,y,z)
                vel_init[ir,it,ip] = vel

    # write out in hdf5 format

    fout_init = h5py.File('models/model_N%d_%d_%d.h5'%(n_rtp[0],n_rtp[1],n_rtp[2]), 'w')

    # write out the arrays eta_init, xi_init, zeta_init, fun_init, a_init, b_init, c_init, f_init
    fout_init.create_dataset('eta', data=eta_init)
    fout_init.create_dataset('xi', data=xi_init)
    fout_init.create_dataset('zeta', data=zeta_init)
    fout_init.create_dataset('vel', data=vel_init)

    fout_init.close()


# dummys
src_id      = 0
src_year    = 1998
src_month   = 1
src_day     = 1
src_hour    = 0
src_minute  = 0
src_sec     = 0
src_lat     = t0 / math.pi * 180
src_lon     = p0 / math.pi * 180
src_dep     = R_earth - r0
src_mag     = 3.0
src_name    = "s0"


Nr = 10; Nt = 10; Np = 10
N_rec       = Nr * Nt * Np



fname = 'src_rec_true.dat'
doc_src_rec = open(fname,'w')

doc_src_rec.write('%8d %8d %2d %2d %2d %2d %5.2f %9.4f %9.4f %8.4f %5.2f %5d %s\n'%
                 (src_id,src_year,src_month,src_day,src_hour,src_minute,src_sec,src_lat,src_lon,src_dep,src_mag,N_rec,src_name))

count = 0
for i in range(Nr):
    for j in range(Nt):
        for k in range(Np):
            rec_id      = count
            rec_name    = "r%d"%(count)
            rec_lat     = 32 + (48-32)/(Nt-1) * j   
            rec_lon     = 17 + (38-17)/(Np-1) * k   
            rec_dep     = (0 + 400/(Nr-1) * i)
            rec_ele     = -rec_dep*1000
            phase       = "P"

            r = R_earth - rec_dep
            t = rec_lat/180.0*math.pi
            p = rec_lon/180.0*math.pi
            z = r * math.sin(t)
            x = r * math.cos(t) * math.cos(p)
            y = r * math.cos(t) * math.sin(p)
            vel,time = linear_velocity_time(x,y,z)

            doc_src_rec.write('%8d %8d %6s %9.4f %9.4f %9.1f %s %8.3f\n'%(src_id,rec_id,rec_name,rec_lat,rec_lon,rec_ele,phase,time))

            count+=1
doc_src_rec.close()



        


