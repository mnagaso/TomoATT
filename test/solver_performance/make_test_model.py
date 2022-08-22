# %% [markdown]
# # notebook for create init and true test model

# %%
import numpy as np
import math

R_earth = 6378.1370
# grid
dep1 = 5740.6370
dep2 = 5790.6370
rr1=R_earth - dep1
rr2=R_earth - dep2
tt1=45.0/180*math.pi
tt2=55.0/180*math.pi
pp1=35.0/180*math.pi
pp2=45.0/180*math.pi

n_rtp = [40,40,40]
n_rtp.reverse()
dr = (rr2-rr1)/n_rtp[2]
dt = (tt2-tt1)/n_rtp[1]
dp = (pp2-pp1)/n_rtp[0]
rr = np.array([rr1 + x*dr for x in range(n_rtp[2])])
tt = np.array([tt1 + x*dt for x in range(n_rtp[1])])
pp = np.array([pp1 + x*dp for x in range(n_rtp[0])])

# initial model
#gamma = 0.0
#s0 = 1.0/6.0
#slow_p=0.06
#ani_p=0.0

def WGS84ToCartesian(R, lat, lon):
    # equatorial radius WGS84 major axis
    equRadius = R_earth
    flattening = 1.0 / 298.257222101

    sqrEccentricity = flattening * (2.0 - flattening);

    lat_rad = lat; # input should be already in radian
    lon_rad = lon;
    alt = R - equRadius; #// altitude (height above sea level)

    sinLat = math.sin(lat_rad);
    cosLat = math.cos(lat_rad);
    sinLon = math.sin(lon_rad);
    cosLon = math.cos(lon_rad);

    # Normalized radius
    normRadius = equRadius / math.sqrt(1.0 - sqrEccentricity * sinLat * sinLat);

    x = (normRadius + alt) * cosLat * cosLon;
    y = (normRadius + alt) * cosLat * sinLon;
    z = (normRadius * (1.0 - sqrEccentricity) + alt) * sinLat;
    return [x,y,z]

# source position
rs = R_earth - 5765.6370
ts = 50/180*math.pi
ps = 40/180*math.pi
x0, y0, z0 = WGS84ToCartesian(rs, ts, ps)
model_base = 5.0
gx = 0.01
gy = -0.006
gz = -0.004

eta_init = np.zeros(n_rtp)
xi_init  = np.zeros(n_rtp)
zeta_init = np.zeros(n_rtp)
fun_init = np.zeros(n_rtp)
a_init = np.zeros(n_rtp)
b_init = np.zeros(n_rtp)
c_init = np.zeros(n_rtp)
f_init = np.zeros(n_rtp)
vel_init = np.zeros(n_rtp)
u_true = np.zeros(n_rtp)

# true model
#eta_true = np.zeros(n_rtp)
#xi_true  = np.zeros(n_rtp)
#zeta_true = np.zeros(n_rtp)
#fun_true = np.zeros(n_rtp)
#a_true = np.zeros(n_rtp)
#b_true = np.zeros(n_rtp)
#c_true = np.zeros(n_rtp)
#f_true = np.zeros(n_rtp)

for ir in range(n_rtp[2]):
    for it in range(n_rtp[1]):
        for ip in range(n_rtp[0]):
            x, y, z = WGS84ToCartesian(rr[ir], tt[it], pp[ip])

            #eta_init[ip,it,ir] = 0.0
            #xi_init[ip,it,ir]  = 0.0
            #zeta_init[ip,it,ir] = 0.0
            a_init[ip,it,ir] = 1.0
            b_init[ip,it,ir] = 1.0 - 2.0*xi_init[ip,it,ir]
            c_init[ip,it,ir] = 1.0 + 2.0*xi_init[ip,it,ir]
            f_init[ip,it,ir] = -2.0 * eta_init[ip,it,ir]

            vel = model_base + gx*(x - x0) + gy*(y - y0) + gz*(z - z0); #/ m/s
            vel_init[ip,it,ir] = vel
            fun_init[ip,it,ir] = 1.0/vel
            normg            = math.sqrt(gx*gx + gy*gy + gz*gz);
            dist             = math.sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0));
            arcz             = 1.0+0.5*fun_init[ip,it,ir]*(1.0/model_base)*normg*normg*dist*dist;
            arccosh          = math.log(arcz + math.sqrt(arcz*arcz - 1.0));
            u_true[ip,it,ir] = arccosh/normg;




# %%
r_earth = 6378.1370
print("depminmax {} {}".format(rr1,rr2))
print("src depth {}".format(rs))

# %%
# write out
import h5py

#dtype_out = np.single
dtype_out = np.double

fout_init = h5py.File('test_model_init.h5', 'w')

# write out the arrays eta_init, xi_init, zeta_init, fun_init, a_init, b_init, c_init, f_init
fout_init.create_dataset('eta',   data=eta_init.T.astype(dtype_out), dtype=dtype_out)
fout_init.create_dataset('xi',    data=xi_init.T.astype(dtype_out), dtype=dtype_out)
fout_init.create_dataset('zeta',  data=zeta_init.T.astype(dtype_out), dtype=dtype_out)
fout_init.create_dataset('fun',   data=fun_init.T.astype(dtype_out), dtype=dtype_out)
fout_init.create_dataset('fac_a', data=a_init.T.astype(dtype_out), dtype=dtype_out)
fout_init.create_dataset('fac_b', data=b_init.T.astype(dtype_out), dtype=dtype_out)
fout_init.create_dataset('fac_c', data=c_init.T.astype(dtype_out), dtype=dtype_out)
fout_init.create_dataset('fac_f', data=f_init.T.astype(dtype_out), dtype=dtype_out)
fout_init.create_dataset('vel',   data=vel_init.T.astype(dtype_out), dtype=dtype_out)
fout_init.create_dataset('u',     data=u_true.T.astype(dtype_out), dtype=dtype_out)

fout_init.close()


# %%
import random
# write src_rec_file

# dummys
year_dummy = 1998
month_dummy = 1
day_dummy = 1
hour_dummy = 0
minute_dummy = 0
second_dummy = 0
mag_dummy = 3.0
id_dummy = 1000
st_name_dummy = 'AAAA'
phase_dummy = 'P'
dist_dummy = 100.0
arriv_t_dummy = 0.0

tt1deg = tt1 * 180.0/math.pi
tt2deg = tt2 * 180.0/math.pi
pp1deg = pp1 * 180.0/math.pi
pp2deg = pp2 * 180.0/math.pi


n_src = 1
n_rec = [2000]


lines = []

# create dummy src
for i_src in range(n_src):
    # define one point in the domain
    dep = 5765.6370
    lon = 40.0
    lat = 50.0

    src = [i_src, year_dummy, month_dummy, day_dummy, hour_dummy, minute_dummy, second_dummy, lat, lon, dep, mag_dummy, n_rec[i_src], id_dummy]
    lines.append(src)

    # create dummy station
    for i_rec in range(n_rec[i_src]):
        #elev_rec = random.uniform(dep1+0.01,dep1)*-1000.0 # elevation in m
        elev_rec = dep1*-1000.0
        lon_rec  = random.uniform(pp1deg,pp2deg)
        lat_rec  = random.uniform(tt1deg,tt2deg)

        rec = [i_src, i_rec, st_name_dummy, lat_rec, lon_rec, elev_rec, phase_dummy, dist_dummy, arriv_t_dummy]
        lines.append(rec)



# write out ev_arrivals file
fname = 'src_rec_test.dat'

with open(fname, 'w') as f:
    for line in lines:
        for elem in line:
            f.write('{} '.format(elem))
        f.write('\n')

# %%



