# %%
def relative_id_nosort(line_target,line_ref):
    N_target = len(line_target)
    id_target = []
    N_ref = len(line_ref)
    for ip in range(0,N_target):
        val_target = line_target[ip]
        for i in range(0,N_ref-1):
            if (line_ref[i]<=val_target and line_ref[i+1]>val_target):
                id_target.append(i)
                break
        if (len(id_target)<ip+1): 
            id_target.append(-1)    
        if (len(id_target)!=ip+1): 
            print('error')
            print(len(id_target),ip)
    return id_target 


# %%
import h5py
import numpy as np

def read_model_file_horizontal(grid_fname,model_fname,Nstep,stlat,edlat,stlon,edlon,Nlat,Nlon,all_dep):
    R_earth = 6371.0

    fmodel  = h5py.File(model_fname, 'r')
    fgrid   = h5py.File(grid_fname, 'r')

    print(fmodel['model'].keys())

    vel_inv  = np.array(fmodel['model']["vel_inv_%s"%(Nstep)])
    xi_inv  = np.array(fmodel['model']["xi_inv_%s"%(Nstep)])
    eta_inv  = np.array(fmodel['model']["eta_inv_%s"%(Nstep)])

    grid_r = np.array(fgrid["Mesh"]["node_coords_r"])
    grid_t = np.array(fgrid["Mesh"]["node_coords_t"])
    grid_p = np.array(fgrid["Mesh"]["node_coords_p"])
    
    fmodel.close()
    fgrid.close()

    model={}
    dir_lat={};dir_lon={};dir_dep={}
    np_lat=[];np_lon=[];np_dep=[]
    for i in range(len(vel_inv)):
        lat  = '%.8g'%(grid_t[i]); 
        lon  = '%.8g'%(grid_p[i]); 
        dep  = '%.8g'%(R_earth - grid_r[i]); 
        vel = vel_inv[i]
        xi   = xi_inv[i]
        eta  = eta_inv[i]
        sigma = vel_inv[i] 

        model['%s_%s_%s'%(lat,lon,dep)]=(vel,xi,eta,sigma)
        dir_lat[lat] = float(lat)
        dir_lon[lon] = float(lon)
        dir_dep[dep] = float(dep)



    for i in dir_lat:
        np_lat.append(dir_lat[i])
    np_lat.sort()
    for i in dir_lon:
        np_lon.append(dir_lon[i])
    np_lon.sort()
    for i in dir_dep:
        np_dep.append(dir_dep[i])
    np_dep.sort()



    all_lat=np.linspace(stlat,edlat,Nlat,dtype=float)
    all_lon=np.linspace(stlon,edlon,Nlon,dtype=float)



    np_lat_id=[];np_lon_id=[];np_dep_id=[]

    np_lat_id = relative_id_nosort(all_lat,np_lat) 
    np_lon_id = relative_id_nosort(all_lon,np_lon) 
    np_dep_id = relative_id_nosort(all_dep,np_dep) 

    target_dep_value = np.full((len(all_dep),Nlat*Nlon,7),np.nan)    # dep,lat,lon,vel,xi,eta,simga
    for idep in range(0,len(all_dep)):
        dep=all_dep[idep]

        for i in range(0,Nlat):
            lat=all_lat[i]
            
            for j in range(0,Nlon):
                lon=all_lon[j]

                (phy_lat, phy_lon) = (lat,lon)


                #try:

                # 确定了目标点，开始寻找周围八个点
                lat1=np_lat[np_lat_id[i]];lat2=np_lat[np_lat_id[i]+1]
                lon1=np_lon[np_lon_id[j]];lon2=np_lon[np_lon_id[j]+1]
                dep1=np_dep[np_dep_id[idep]];dep2=np_dep[np_dep_id[idep]+1]

                if (not(lat>=lat1 and lat<lat2 and lon>=lon1 and lon<lon2 and dep>=dep1 and dep<dep2)):
                    target_dep_value[idep,i*Nlon+j,:] = (dep,lat,lon,np.nan,np.nan,np.nan,np.nan)
                    # print(lat,lon,dep,"out of region")
                else:
                    # 三个方向比例
                    rdep=(dep-dep1)/(dep2-dep1); rdep=min(rdep,1);rdep=max(0,rdep)
                    rlat=(lat-lat1)/(lat2-lat1); rlat=min(rlat,1);rlat=max(0,rlat)
                    rlon=(lon-lon1)/(lon2-lon1); rlon=min(rlon,1);rlon=max(0,rlon)

                    vel=0
                    vel+=rlat*rlon*rdep*model['%.8g_%.8g_%.8g'%(lat2,lon2,dep2)][0]
                    vel+=rlat*rlon*(1-rdep)*model['%.8g_%.8g_%.8g'%(lat2,lon2,dep1)][0]
                    vel+=rlat*(1-rlon)*rdep*model['%.8g_%.8g_%.8g'%(lat2,lon1,dep2)][0]
                    vel+=rlat*(1-rlon)*(1-rdep)*model['%.8g_%.8g_%.8g'%(lat2,lon1,dep1)][0]
                    vel+=(1-rlat)*rlon*rdep*model['%.8g_%.8g_%.8g'%(lat1,lon2,dep2)][0]
                    vel+=(1-rlat)*rlon*(1-rdep)*model['%.8g_%.8g_%.8g'%(lat1,lon2,dep1)][0]
                    vel+=(1-rlat)*(1-rlon)*rdep*model['%.8g_%.8g_%.8g'%(lat1,lon1,dep2)][0]
                    vel+=(1-rlat)*(1-rlon)*(1-rdep)*model['%.8g_%.8g_%.8g'%(lat1,lon1,dep1)][0]


                    xi=0
                    xi+=rlat*rlon*rdep*model['%.8g_%.8g_%.8g'%(lat2,lon2,dep2)][1]
                    xi+=rlat*rlon*(1-rdep)*model['%.8g_%.8g_%.8g'%(lat2,lon2,dep1)][1]
                    xi+=rlat*(1-rlon)*rdep*model['%.8g_%.8g_%.8g'%(lat2,lon1,dep2)][1]
                    xi+=rlat*(1-rlon)*(1-rdep)*model['%.8g_%.8g_%.8g'%(lat2,lon1,dep1)][1]
                    xi+=(1-rlat)*rlon*rdep*model['%.8g_%.8g_%.8g'%(lat1,lon2,dep2)][1]
                    xi+=(1-rlat)*rlon*(1-rdep)*model['%.8g_%.8g_%.8g'%(lat1,lon2,dep1)][1]
                    xi+=(1-rlat)*(1-rlon)*rdep*model['%.8g_%.8g_%.8g'%(lat1,lon1,dep2)][1]
                    xi+=(1-rlat)*(1-rlon)*(1-rdep)*model['%.8g_%.8g_%.8g'%(lat1,lon1,dep1)][1]


                    eta=0
                    eta+=rlat*rlon*rdep*model['%.8g_%.8g_%.8g'%(lat2,lon2,dep2)][2]
                    eta+=rlat*rlon*(1-rdep)*model['%.8g_%.8g_%.8g'%(lat2,lon2,dep1)][2]
                    eta+=rlat*(1-rlon)*rdep*model['%.8g_%.8g_%.8g'%(lat2,lon1,dep2)][2]
                    eta+=rlat*(1-rlon)*(1-rdep)*model['%.8g_%.8g_%.8g'%(lat2,lon1,dep1)][2]
                    eta+=(1-rlat)*rlon*rdep*model['%.8g_%.8g_%.8g'%(lat1,lon2,dep2)][2]
                    eta+=(1-rlat)*rlon*(1-rdep)*model['%.8g_%.8g_%.8g'%(lat1,lon2,dep1)][2]
                    eta+=(1-rlat)*(1-rlon)*rdep*model['%.8g_%.8g_%.8g'%(lat1,lon1,dep2)][2]
                    eta+=(1-rlat)*(1-rlon)*(1-rdep)*model['%.8g_%.8g_%.8g'%(lat1,lon1,dep1)][2]

                    sigma = 0
                    sigma+=rlat*rlon*rdep*model['%.8g_%.8g_%.8g'%(lat2,lon2,dep2)][3]
                    sigma+=rlat*rlon*(1-rdep)*model['%.8g_%.8g_%.8g'%(lat2,lon2,dep1)][3]
                    sigma+=rlat*(1-rlon)*rdep*model['%.8g_%.8g_%.8g'%(lat2,lon1,dep2)][3]
                    sigma+=rlat*(1-rlon)*(1-rdep)*model['%.8g_%.8g_%.8g'%(lat2,lon1,dep1)][3]
                    sigma+=(1-rlat)*rlon*rdep*model['%.8g_%.8g_%.8g'%(lat1,lon2,dep2)][3]
                    sigma+=(1-rlat)*rlon*(1-rdep)*model['%.8g_%.8g_%.8g'%(lat1,lon2,dep1)][3]
                    sigma+=(1-rlat)*(1-rlon)*rdep*model['%.8g_%.8g_%.8g'%(lat1,lon1,dep2)][3]
                    sigma+=(1-rlat)*(1-rlon)*(1-rdep)*model['%.8g_%.8g_%.8g'%(lat1,lon1,dep1)][3]

                    target_dep_value[idep,i*Nlon+j,:] = (dep,phy_lat,phy_lon,vel,xi,eta,sigma)
    return target_dep_value


# %%
def ave_model(target_model):
    Ndep  = np.size(target_model,0)
    Ngrid = np.size(target_model,1)
    average_model = np.copy(target_model)    # dep,lat,lon,slowness,xi,eta,vel
    
    for idep in range(Ndep):
        ave_value = 0
        count = 0
        for igrid in range(Ngrid):
            if(np.isnan(target_model[idep,igrid,3])):
                # nan, 
                pass
            else:
                print(target_model[idep,igrid,:])
                print(ave_value,count)
                ave_value = ave_value + target_model[idep,igrid,3]
                count = count + 1
        
        ave_value = ave_value/count

        for igrid in range(Ngrid):
            if(np.isnan(target_model[idep,igrid,3])):
                # nan, 
                pass
            else:
                average_model[idep,igrid,3] = ave_value

    return average_model

# %%
# 初始化，类的定义和声明

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
        self.Nt = 0     # 每个地震的有效 pick 数         
        self.Ndt = 0       # 每个地震的 双差数量
        self.t = {}     # stname -> (stname, phase, time, stwt)
        self.dt = {}    # stname1+stname2 -> (stname1, stname2, phase, dif_time, stwt)
        self.azi_gap = 360.0
        self.misfit = {}    # 走时残差，真实数据和合成数据只差，用于评估。  stname or stname1+stname2 -> residual
        self.tag    = {}    # 用于描述地震本身的属性，例如 是否HQ: is_HQ; 残差均值: res_mean; 残差标准差 rec_sd; weight

class Station():
    def __init__(self):
        self.name = "nan"
        # self.network = "nan"
        self.id = -1
        self.lat = 0.0
        self.lon = 0.0
        self.ele = 0.0
        # self.Nev = 0       # 每个台站涉及到的地震 picks 数量


        




# %%
# 函数：read_src_rec_file(fname) 读取到时数据
def read_src_rec_file(fname):
    ev_info = {}
    st_info = {}

    doc = open(fname,'r')
    doc_input = doc.readlines()
    doc.close()

    cc = 0
    for info in doc_input:
        tmp=info.split()
        if (cc == 0):   # event line
            ev = Event()
            #  1 2000  1  2 20 28  37.270   38.2843     39.0241  11.00  3.60    8   1725385
            # id_ev   = int(tmp[0])
            ev.id   = len(ev_info)
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
            ev.Nt   = 0
            ev.t    = {}
            ev.Ndt  = 0
            ev.dt   = {}
            ndata   = int(tmp[11])
            name_ev = tmp[12]
            ev.name = name_ev
            cc += 1


        else:   # data line
            #  1      1      MYA       38.3261     38.4253   1050.0000  P   52.46   6.630 weight
            if(len(tmp) < 10):
                name_st = tmp[2]
                phase   = tmp[6]
                if (phase == "PG"):
                    phase = "Pg"
                if (phase == "PB"):
                    phase = "Pb"
                if (phase == "PN"):
                    phase = "Pn"    

                time    = float(tmp[7])
                if (not name_st in st_info):
                    st = Station()
                    st.name = name_st
                    st.id   = len(st_info)
                    st.lat  = float(tmp[3])
                    st.lon  = float(tmp[4])
                    st.ele  = float(tmp[5])
                    st_info[name_st] = st
                ev.t["%s+%s"%(name_st,phase)] = [name_st,phase,time]
                ev.Nt += 1
            else:   
            #  evid     stid1      stname1      lat1 lon1 eve1 stid2 stname2 lat2 lon2 ele2 phase diftime weight 
                name_st1 = tmp[2]
                if (not name_st1 in st_info):
                    st = Station()
                    st.name = name_st1
                    st.id   = len(st_info)
                    st.lat  = float(tmp[3])
                    st.lon  = float(tmp[4])
                    st.ele  = float(tmp[5])
                    st_info[name_st1] = st
 
                name_st2 = tmp[7]
                if (not name_st2 in st_info):
                    st = Station()
                    st.name = name_st2
                    st.id   = len(st_info)
                    st.lat  = float(tmp[8])
                    st.lon  = float(tmp[9])
                    st.ele  = float(tmp[10])
                    st_info[name_st2] = st
 
                phase = tmp[11]
                dif_time = float(tmp[12])
                ev.dt["%s+%s+%s"%(name_st1,name_st2,phase)] = [name_st1,name_st2,phase,dif_time]
                ev.Ndt += 1

            if (cc == ndata):
                cc = 0
                ev_info[name_ev] = ev
            else:
                cc += 1
    return [ev_info,st_info]

# 函数： data_lon_lat_dep_wt_ev(ev_info) 输出地震的 [lon,lat,dep,weight]
def data_lon_lat_dep_wt_ev(ev_info):
    lat = []
    lon = []
    dep = []
    wt = []
    for key in ev_info:
        lat.append(ev_info[key].lat)
        lon.append(ev_info[key].lon)
        dep.append(ev_info[key].dep)
        try:
            wt.append(ev_info[key].tag["weight"])
        except:
            wt.append(1.0)
    return [np.array(lon),np.array(lat),np.array(dep),np.array(wt)]

# 函数： data_lon_lat_dep_st(ev_info,st_info) 输出台站的 [lon,lat,dep]
def data_lon_lat_ele_st(ev_info,st_info):
    names = {}
    lat = [] 
    lon = [] 
    ele = []
    for key_ev in ev_info:
        for key_t in ev_info[key_ev].t:
            name_st = ev_info[key_ev].t[key_t][0]
            names[name_st] = name_st

    for name in names:
        lat.append(st_info[name].lat)
        lon.append(st_info[name].lon)
        ele.append(st_info[name].ele)

    return [np.array(lon),np.array(lat),np.array(ele)]


# %%
import sys

args = sys.argv


if (len(args) == 1 and args[0] == "ckb"):
    tag = "ckb"
elif (len(args) == 1 and args[0] == "inv"):
    tag = "inv"
else:
    exit('Error input')

stlat=30.0
edlat=32.0
stlon=30.0
edlon=32.0

# 数据区域范围
data_stlat = stlat
data_edlat = edlat
data_stlon = stlon
data_edlon = edlon

# 投影方式
projection = "M11c"

# 离散点数
Nlat=101
Nlon=101

# 画图指定深度
all_dep = [10,30]


print('loading file ...')
# 读取模型文件

# test_path = "ckb_N3_10_6"
#test_path = "t5_ckb_N3_8_4_增加各向异性+t4"



# ckb model
if tag == "ckb":
    test_path = "ckb_model"
    target_fname = "OUTPUT_FILES/OUTPUT_FILES_swap_signal/out_data_sim_group_0.h5"   # 检测板反演的最终模型
    target_Nstep = "0001"
    ref_fname = "OUTPUT_FILES/OUTPUT_FILES_swap_inv_abs/out_data_sim_group_0.h5"  # 检测板反演的初始模型
    ref_Nstep = "0001"
    grid_fname = "OUTPUT_FILES/OUTPUT_FILES_swap_signal/out_data_grid.h5"                # 网格文件
    fname = "OUTPUT_FILES/OUTPUT_FILES_swap_signal/src_rec_file_forward.dat"
    output = "%s"%(test_path)


# ckb inv model
if tag == "inv":
    test_path = "OUTPUT_FILES_swap_inv_abs"
    target_fname = "OUTPUT_FILES/%s/out_data_sim_group_0.h5"%(test_path)   # 检测板反演的最终模型
    target_Nstep = "0002"
    ref_fname = "OUTPUT_FILES/%s/out_data_sim_group_0.h5"%(test_path)  # 检测板反演的初始模型
    ref_Nstep = "0001"
    grid_fname = "OUTPUT_FILES/%s/out_data_grid.h5"%(test_path)                # 网格文件
    output = "%s_%s"%(test_path,target_Nstep)
    fname = "OUTPUT_FILES/OUTPUT_FILES_swap_signal/src_rec_file_forward.dat"


# 插值获得画图数据
target_model = read_model_file_horizontal(grid_fname,target_fname,target_Nstep,data_stlat,data_edlat,data_stlon,data_edlon,Nlat,Nlon,all_dep)
ref_model = read_model_file_horizontal(grid_fname,ref_fname,ref_Nstep,data_stlat,data_edlat,data_stlon,data_edlon,Nlat,Nlon,all_dep)     # initiai model
# ref_model = ave_model(target_model)

[ev_info, st_info] = read_src_rec_file(fname)

# %%
# -------------------- 深度剖面画图 -------------------------
import pygmt
from pygmt.clib import Session



dlat = (edlat-stlat)/(Nlat-1)
dlon = (edlon-stlon)/(Nlat-1)

# 图片位置
xshift = [ 3,   13,   13,  -26, 13, 13, -26, 13, 13 ]
yshift = [ 80,   0,    0,  -15,  0,  0, -15, 0,  0  ]

frame = [
    ["xa2","ya2","NsWe"],["xa2","ya2","Nswe"],["xa2","ya2","NswE"],
    ["xa2","ya2","nsWe"],["xa2","ya2","nswe"],["xa2","ya2","nswE"],
    ["xa2","ya2","nSWe"],["xa2","ya2","nSwe"],["xa2","ya2","nSwE"],
]

subfigure = ['(a)', '(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)']


# %%
# -------------------------- 速度结构 --------------------------------

vel_per     = 6     # 颜色棒范围

# 开始画图
with pygmt.clib.Session() as session:
    session.call_module('gmtset', 'FONT 20p')
fig = pygmt.Figure()
pygmt.config(IO_SEGMENT_MARKER="<<<")

for i in range(0,len(all_dep)):
    print('plotting: dep %3d'%(all_dep[i]))
    pert = np.full(Nlat*Nlon,-100.0)
    lon = np.full(Nlat*Nlon,-100.0)
    lat = np.full(Nlat*Nlon,-100.0)

    for igrid in range(0,Nlat*Nlon):
        pert[igrid] = (target_model[i,igrid,3]-ref_model[i,igrid,3])/ref_model[i,igrid,3]*100
        lon[igrid] = target_model[i,igrid,2]
        lat[igrid] = target_model[i,igrid,1]

    # 生成 grd 文件
    grid = pygmt.surface(x=lon, y=lat, z=pert, spacing=(dlon/2, dlat/2), region=[stlon, edlon, stlat, edlat])

    # ----------- 构造 colorbar ---------------- 
    pygmt.makecpt(cmap="seis", series=[-vel_per, vel_per], background = True)


    # ----------- 画布 ----------------
    fig.shift_origin(xshift=xshift[i],yshift = yshift[i])
    fig.basemap(
        frame=frame[i],                # 坐标轴 -B 的属性  "xa10f5+lLABELNAME" a 表示annotation 的间距, f 表示线的间距, l表示label名称
        projection=projection,  # 投影方式
        region=[stlon,edlon,stlat,edlat],       # 区域范围    
    )

    # ----------- 画 map ----------------
    fig.grdimage(
        frame=frame[i],
        grid = grid,
        projection=projection,  # 投影方式
        region=[stlon,edlon,stlat,edlat],       # 区域范围
    )



    

    # ---------- 画地震台站 ----------------
    [lon_ev,lat_ev,dep_ev,wt_ev] = data_lon_lat_dep_wt_ev(ev_info)
    fig.plot(x = lon_ev, y = lat_ev, size = [0.03]*len(lat_ev),style = "c", fill = "black")

    # ----------- 画区域台站 ----------------
    [lon_st,lat_st,ele_st] = data_lon_lat_ele_st(ev_info,st_info)
    fig.plot(x = lon_st, y = lat_st, size = [0.3]*len(lat_st),style = "t",fill = "black",pen = "0.5p,white")

    # ----------- 一些注记 ----------------
    fig.text(text="Depth = %3d km"%(all_dep[i]), x=38.5, y=35.5, font="18p,Helvetica-Bold,black", fill="white")
    fig.text(text="Vp = %4.2f km/s"%(ref_model[i,0,6]), x=38.5, y=36.0, font="18p,Helvetica-Bold,black", fill="white")
    # fig.text(text=subfigure[i], x=98, y=20.5, font="18p,Helvetica-Bold,black", fill="white")

fig.shift_origin(xshift = 2, yshift = -2)
fig.colorbar(frame = ["a2f1","x+lln(m%d/m00)"%(target_Nstep),"y+l(%)"], position="+w7c/0.5c+h", ) # +e,默认是双箭头，f表示forward，b表示background ，w表示长宽，h表示水平

fig.savefig('img/%s_vel.jpg'%(output))

# %%
# -------------------------- 各向异性结构 --------------------------------
import math

ani_per = 0.03  # 颜色棒范围

RAD2DEG = 180/math.pi

# 开始画图
with pygmt.clib.Session() as session:
    session.call_module('gmtset', 'FONT 20p')
fig = pygmt.Figure()
pygmt.config(IO_SEGMENT_MARKER="<<<")     

for i in range(0,len(all_dep)):
    
    # data arrangement
    print('plotting: dep %3d'%(all_dep[i]))
    perturbation = np.full(Nlat*Nlon,-100.0)
    lon = np.full(Nlat*Nlon,-100.0)
    lat = np.full(Nlat*Nlon,-100.0)
    eta = np.full(Nlat*Nlon,-100.0)
    xi = np.full(Nlat*Nlon,-100.0)
    epsilon = np.full(Nlat*Nlon,-100.0)
    for igrid in range(0,Nlat*Nlon):
        lon[igrid] = target_model[i,igrid,2]
        lat[igrid] = target_model[i,igrid,1]
        eta[igrid] = target_model[i,igrid,5]
        xi[igrid] = target_model[i,igrid,4]
        epsilon[igrid] = math.sqrt(xi[igrid]**2 + eta[igrid]**2)

    lat_gap = 5; lon_gap = 5
    N_lat_ani = int(Nlat/lat_gap)
    N_lon_ani = int(Nlon/lon_gap)
    
    lat_ani = []
    lon_ani = []
    xi_ani =  []
    eta_ani = []
    sigma = []
    epsilon_ani = []
    psi_ani =  []
    igrid_ani = 0
    for i_lat in range(N_lat_ani):
        for i_lon in range(N_lon_ani):
            igrid = i_lat*lat_gap*Nlon + i_lon*lon_gap

            tmp_xi = target_model[i,igrid,4]
            tmp_eta = target_model[i,igrid,5]
            tmp_epsilon = math.sqrt(tmp_xi**2 + tmp_eta**2)
            if (tmp_epsilon > 0.01):
                lat_ani.append(target_model[i,igrid,1]) 
                lon_ani.append(target_model[i,igrid,2])
                
                xi_ani.append(tmp_xi)
                eta_ani.append(tmp_eta)
                epsilon_ani.append(tmp_epsilon)
                if (tmp_epsilon == 0):
                    psi_ani.append(0)
                else:
                    tmp = math.acos(tmp_xi/tmp_epsilon)
                    if (tmp_eta>0):
                        psi_ani.append(tmp/2*RAD2DEG)
                        # psi_ani.append(tmp/2*RAD2DEG)
                    else:
                        psi_ani.append((2*math.pi-tmp)/2*RAD2DEG)
                        # psi_ani.append((2*math.pi-tmp)/2*RAD2DEG)

    # 生成 grd 文件
    grid = pygmt.surface(x=lon, y=lat, z=epsilon, spacing=(dlon/2, dlat/2), region=[stlon, edlon, stlat, edlat])

    # ----------- 构造 colorbar ---------------- 
    pygmt.makecpt(cmap="cool", series=[0, ani_per],background=True)


    # ----------- 画布 ----------------
    fig.shift_origin(xshift=xshift[i],yshift = yshift[i])
    fig.basemap(
        frame=frame[i],                # 坐标轴 -B 的属性  "xa10f5+lLABELNAME" a 表示annotation 的间距, f 表示线的间距, l表示label名称
        projection=projection,  # 投影方式
        region=[stlon,edlon,stlat,edlat],       # 区域范围
    )


    # ----------- 画 map ----------------
    fig.grdimage(
        frame=frame[i],
        grid = grid,
        projection=projection,  # 投影方式
        region=[stlon,edlon,stlat,edlat],       # 区域范围
    )



    try:
        # fig.plot( x = lon_ani, y = lat_ani, direction = [psi_ani,np.full(len(psi_ani),0.2)], pen = "1p,red", style ="v0c")
        # fig.plot( x = lon_ani, y = lat_ani, direction = [psi_ani,np.full(len(psi_ani),-0.2)], pen = "1p,red", style ="v0c")

        # 0.005 - 0.05  ->  0.2 - 0.6
        length = []
        for tmp_i in range(0,len(epsilon_ani)):
            length.append(min((epsilon_ani[tmp_i]-0.005)/0.045*0.4+0.2,0.6))
        length = np.array(length)

        fig.plot( x = lon_ani, y = lat_ani, direction = [psi_ani,np.full(len(psi_ani),length)], pen = "3p,black", style ="v0c")
        fig.plot( x = lon_ani, y = lat_ani, direction = [psi_ani,np.full(len(psi_ani),-length)], pen = "3p,black", style ="v0c")
        
        fig.plot( x = lon_ani, y = lat_ani, direction = [psi_ani,np.full(len(psi_ani),length-0.02)], pen = "2.6p,yellow1", style ="v0c")
        fig.plot( x = lon_ani, y = lat_ani, direction = [psi_ani,np.full(len(psi_ani),-length+0.02)], pen = "2.6p,yellow1", style ="v0c")
    
    except:
        pass


    # ----------- 一些注记 ----------------
    fig.text(text="Depth = %3d km"%(all_dep[i]), x=38.5, y=35.5, font="18p,Helvetica-Bold,black", fill="white")
    # fig.text(text="Vp = %4.2f km/s"%(ave_target_model[i,0,6]), x=39.5, y=36.5, font="18p,Helvetica-Bold,black", fill="white")
    # fig.text(text=subfigure[i], x=98, y=20.5, font="18p,Helvetica-Bold,black", fill="white")

fig.shift_origin(xshift=2,yshift = -2)
fig.colorbar(frame = ["a0.02","x+lm%d"%(target_Nstep)], position="+w7c/0.5c+h") # +e,默认是双箭头，f表示forward，b表示background ，w表示长宽，h表示水平

fig.savefig('img/%s_ani.jpg'%(output))


