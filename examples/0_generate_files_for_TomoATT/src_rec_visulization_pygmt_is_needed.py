
# %%
# -------------------- plot stations and receivers -------------------------
import os

try:
    os.mkdir("img")
except:
    pass

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
import math

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





