# %% [markdown]
# # notebook for create init and true test model

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
import math
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
# generate src_rec_file
import os

try:
    os.mkdir("src_rec_files")
except:
    print("dir models exists")

fname = 'src_rec_files/src_rec_config.dat'
doc_src_rec = open(fname,'w')

accept_thd_abs = 0.3
accept_thd_cs = 0.2
accept_thd_cr = 0.9

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
            if (accept > accept_thd_cs and azi_gap < azimuth_thd):
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
            if (accept > accept_thd_cr and azi_gap < azimuth_thd):
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


