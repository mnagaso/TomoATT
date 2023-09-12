# %% [markdown]
# # notebook for create syn and obs src_rec data

# %%
doc = open("OUTPUT_FILES/OUTPUT_FILES_signal/src_rec_file_forward.dat",'r')
info_srcrec = doc.readlines()
doc.close()


doc_obs = open("src_rec_obs.dat",'w')
doc_obs_info = []

# %%
import os
import math
from obspy import UTCDateTime
import numpy as np
import copy

class Event():
    def __init__(self):
        self.name = "nan"
        self.id = -1
        self.lat = 0.0
        self.lon = 0.0
        self.dep = 0.0
        self.mag = 0.0
        self.ortime = UTCDateTime(1999,1,1,0,0,0)
        self.ndata = 0     # 每个地震的有效 pick 数     
        self.data_info = []    



# %%
cc = 0
ev_info = []
for line in info_srcrec:
    tmp = line.split()
    if (cc == 0):   # event line
        ev = Event()
        #  1 2000  1  2 20 28  37.270   38.2843     39.0241  11.00  3.60    8   1725385
        # id_ev   = int(tmp[0])
        ev.id   = int(tmp[0])
        year    = int(tmp[1])
        month   = int(tmp[2])
        day     = int(tmp[3])
        hour    = int(tmp[4])
        minute  = int(tmp[5])
        second  = float(tmp[6])
        ev.ortime = UTCDateTime(year,month,day,hour,minute,0) + second
        ev.lat  = float(tmp[7])
        ev.lon  = float(tmp[8])
        ev.dep  = float(tmp[9])
        ev.mag  = float(tmp[10])
        ev.ndata   = int(tmp[11])            
        ev.name = tmp[12]
        cc += 1
    else:   # data line
        ev.data_info.append(line)
        if (cc == ev.ndata):
            cc = 0
            ev_info.append(ev)
        else:
            cc += 1

# %%
# add ortime into earthquake, adding deviations to the data
import random
seed_value = 42
random.seed(seed_value)

tau_shift = []
lat_shift = []
lon_shift = []
dep_shift = []

Nev = len(ev_info)

for i in range(Nev):
    tau_shift.append(random.uniform(-1.5,1.5))
    lat_shift.append(random.uniform(-0.2,0.2))
    lon_shift.append(random.uniform(-0.2,0.2))
    dep_shift.append(random.uniform(-5.0,5.0))


for ev in ev_info:
    ortime_obs  = ev.ortime + tau_shift[ev.id]
    lat_obs     = ev.lat    + lat_shift[ev.id]
    lon_obs     = ev.lon    + lon_shift[ev.id]
    dep_obs     = ev.dep    + dep_shift[ev.id]
    doc_obs.write('%8d %8d %2d %2d %2d %2d %5.2f %9.4f %9.4f %8.4f %5.2f %5d %s\n'%
                 (ev.id,ortime_obs.year,ortime_obs.month,ortime_obs.day,ortime_obs.hour,ortime_obs.minute,ortime_obs.second+ortime_obs.microsecond/1000000,
                  lat_obs,lon_obs,dep_obs,ev.mag,ev.ndata,ev.name))

    for data in ev.data_info:
        tmp = data.split()
        if(len(tmp) < 10):  # abs
            # src.id,rec.id,rec.name,rec.lat,rec.lon,rec.ele,phase,0.0
            evid    = tmp[0]; recid  = tmp[1]; recname = tmp[2]; 
            reclat  = tmp[3]; reclon = tmp[4]; recele  = tmp[5]; 
            phase   = tmp[6]; 
            time    = float(tmp[7]) - tau_shift[ev.id]

            doc_obs.write('%8s %8s %6s %9s %9s %9s %s %8.3f\n'%(evid,recid,recname,reclat,reclon,recele,phase,time))

        else:       # dif
            if(tmp[11].__contains__("cr")):       # common receiver
                # src.id,rec.id,rec.name,rec.lat,rec.lon,rec.ele,src2.id,src2.name,src2.lat,src2.lon,src2.dep,phase,0.0
                evid    = tmp[0]; rec1id  = tmp[1]; rec1name = tmp[2]; 
                rec1lat = tmp[3]; rec1lon = tmp[4]; rec1ele  = tmp[5]; 
                ev2id   = tmp[6]; ev2name = tmp[7]; 
                ev2lat  = float(tmp[8])     + lat_shift[int(ev2id)]; 
                ev2lon  = float(tmp[9])     + lon_shift[int(ev2id)]; 
                ev2dep  = float(tmp[10])    + dep_shift[int(ev2id)]; 
                phase   = tmp[11]; 
                time    = float(tmp[12]) - tau_shift[int(evid)] + tau_shift[int(ev2id)]
                doc_obs.write('%8s %8s %6s %9s %9s %9s %8s %6s %9.4f %9.4f %9.4f %s %8.3f\n'%
                                 (evid,rec1id,rec1name,rec1lat,rec1lon,rec1ele,ev2id,ev2name,ev2lat,ev2lon,ev2dep,phase,time))

            elif(tmp[11].__contains__("cs")):    # common source
                # src.id,rec.id,rec.name,rec.lat,rec.lon,rec.ele,src2.id,src2.name,src2.lat,src2.lon,src2.dep,phase,0.0
                evid    = tmp[0]; rec1id   = tmp[1]; rec1name = tmp[2]; 
                rec1lat = tmp[3]; rec1lon  = tmp[4]; rec1ele  = tmp[5]; 
                rec2id  = tmp[6]; rec2name = tmp[7]; 
                rec2lat = tmp[8]; rec2lon  = tmp[9]; rec2ele  = tmp[10]; 
                phase   = tmp[11]; 
                time    = float(tmp[12])
                doc_obs.write('%8s %8s %6s %9s %9s %9s %8s %6s %9s %9s %9s %s %8.3f\n'%
                                 (evid,rec1id,rec1name,rec1lat,rec1lon,rec1ele,rec2id,rec2name,rec2lat,rec2lon,rec2ele,phase,time))

doc_obs.close()


