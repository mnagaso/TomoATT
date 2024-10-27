# %%
# initilization 初始化，类的定义和声明

import os
import math
from obspy import UTCDateTime
import numpy as np
import copy

class Event():      # class of earthquake 地震事件类
    def __init__(self):
        self.name = "nan"       # evname1 地震名称，推荐为地震
        self.id = -1
        self.lat = 0.0
        self.lon = 0.0
        self.dep = 0.0
        self.mag = 0.0
        self.ortime = UTCDateTime(1999,1,1,0,0,0)
        self.Nt = 0             # Number of the abs traveltime of earthquake 每个地震的 绝对到时 数量
        self.Ncs_dt = 0         # Number of the commmon source differential traveltime of earthquake   每个地震的 共源双差到时 数量
        self.Ncr_dt = 0         # Number of the commmon receiver differential traveltime of earthquake 每个地震的 共台站双差到时 数量
        self.t = {}             # stname1+phase -> (stname1, phase, time, data_weight)
        self.cs_dt = {}         # stname1 + stname2 + phase -> (stname1, stname2, phase, dif_time, data_weight)
        self.cr_dt = {}         # stname1 + evname2 + phase -> (stname1, evname2, phase, dif_time, data_weight)
        self.azi_gap = 360.0    # the max azimuthal gap of this earthquake 该地震的最大方位角间隔
        self.misfit = {}        # traveltime residual of the data 走时残差，真实数据和合成数据只差，用于评估。  stname or stname1+stname2 or stname1+evname2 -> residual
        self.tag    = {}        # additional tags for the earthquake, e.g., azi_gap, weight. (azimuthal gap, weight of the earthquake)
     
class Station():
    def __init__(self):
        self.name = "nan"       # stname1, recommend: network.stname 台站名称, 推荐为 台网.台站名
        self.id = -1
        self.lat = 0.0
        self.lon = 0.0
        self.ele = 0.0
        self.tag = {}           # additional tags for the station, e.g., wright, 用于描述台站本身的属性，例如 台站权重：weight




# %% [markdown]
# 函数类：一些数学上的基本函数，用于辅助处理数据
# 
# Functions: some basic auxiliary functions for processing data

# %%
# 函数： cal_dis(lat1, lon1,lat2, lon2) (千米)， cal_azimuth(lat1, lon1, lat2, lon2) (°) 计算震中距和方位角 compute epicentral distance (km) and azimuth (degree)

def cal_dis(lat1, lon1,lat2, lon2, R = 6371):
    latitude1 = (math.pi/180)*lat1
    latitude2 = (math.pi/180)*lat2
    longitude1 = (math.pi/180)*lon1
    longitude2= (math.pi/180)*lon2
    #因此AB两点的球面距离为:{arccos[sinb*siny+cosb*cosy*cos(a-x)]}*R
    #地球半径
    if((lat1-lat2)**2+(lon1-lon2)**2<0.000001):
        return 0

    d = math.acos(math.sin(latitude1)*math.sin(latitude2)+ math.cos(latitude1)*math.cos(latitude2)*math.cos(longitude2-longitude1))/math.pi*180
    return d * 2 * math.pi * R / 360

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
# 函数: 坐标旋转 rotate_src_rec(ev_info,st_info,theta0,phi0,psi): rotate to the new coordinate, satisfying the center r0,t0,p0 -> r0,0,0 and a anticlockwise angle psi
# 满足中心点 r0,t0,p0 -> r0,0,0 以及 逆时针旋转角度 psi

import numpy as np

RAD2DEG = 180/np.pi
DEG2RAD = np.pi/180
R_earth = 6371.0

# Spherical coordinates to Cartesian coordinate
def rtp2xyz(r,theta,phi):
    x = r * np.cos(theta*DEG2RAD) * np.cos(phi*DEG2RAD)
    y = r * np.cos(theta*DEG2RAD) * np.sin(phi*DEG2RAD)
    z = r * np.sin(theta*DEG2RAD)
    return (x,y,z)

# Cartesian coordinates to Spherical coordinate
def xyz2rtp(x,y,z):
    # theta: -90~90;  phi: -180~180
    r       = np.sqrt(x**2+y**2+z**2)
    theta   = np.arcsin(z/r)
    phi     = np.arcsin(y/r/np.cos(theta))


    idx = np.where((phi > 0) & (x*y < 0))
    phi[idx] = np.pi - phi[idx]
    idx = np.where((phi < 0) & (x*y > 0))
    phi[idx] = -np.pi - phi[idx]
        

    # for i in range(phi.size):
    #     if(phi[i] > 0 and x[i]*y[i] < 0):
    #         phi[i] = np.pi - phi[i]
    #     if(phi[i] < 0 and x[i]*y[i] > 0):
    #         phi[i] = -np.pi - phi[i]
    
    return (r,theta*RAD2DEG,phi*RAD2DEG)

# anti-clockwise rotation along x-axis
def rotate_x(x,y,z,theta):
    new_x = x 
    new_y = y *  np.cos(theta*DEG2RAD) + z * -np.sin(theta*DEG2RAD)
    new_z = y *  np.sin(theta*DEG2RAD) + z *  np.cos(theta*DEG2RAD)
    return (new_x,new_y,new_z)
    
# anti-clockwise rotation along y-axis
def rotate_y(x,y,z,theta):
    new_x = x *  np.cos(theta*DEG2RAD) + z *  np.sin(theta*DEG2RAD)
    new_y = y
    new_z = x * -np.sin(theta*DEG2RAD) + z *  np.cos(theta*DEG2RAD)
    return (new_x,new_y,new_z)

# anti-clockwise rotation along z-axis
def rotate_z(x,y,z,theta):
    new_x = x *  np.cos(theta*DEG2RAD) + y * -np.sin(theta*DEG2RAD)
    new_y = x *  np.sin(theta*DEG2RAD) + y *  np.cos(theta*DEG2RAD)
    new_z = z 
    return (new_x,new_y,new_z)

# spherical Rotation

# rotate to the new coordinate, satisfying the center r0,t0,p0 -> r0,0,0 and a anticlockwise angle psi
def rtp_rotation(t,p,theta0,phi0,psi):
    # step 1: r,t,p -> x,y,z
    (x,y,z) = rtp2xyz(1.0,t,p)

    # step 2: anti-clockwise rotation with -phi0 along z-axis:   r0,t0,p0 -> r0,t0,0
    (x,y,z) = rotate_z(x,y,z,-phi0)

    # step 3: anti-clockwise rotation with theta0 along y-axis:  r0,t0,0 -> r0,0,0
    (x,y,z) = rotate_y(x,y,z,theta0)

    # # step 4: anti-clockwise rotation with psi along x-axis
    (x,y,z) = rotate_x(x,y,z,psi)

    # step 5: x,y,z -> r,t,p
    (new_r,new_t,new_p) = xyz2rtp(x,y,z)
    
    return (new_t,new_p)


def rtp_rotation_reverse(new_t,new_p,theta0,phi0,psi):
    # step 1: r,t,p -> x,y,z
    (x,y,z) = rtp2xyz(1.0,new_t,new_p)

    # step 2: anti-clockwise rotation with -psi along x-axis
    (x,y,z) = rotate_x(x,y,z,-psi)

    # step 3: anti-clockwise rotation with -theta0 along y-axis:  r0,0,0 -> r0,t0,0 
    (x,y,z) = rotate_y(x,y,z,-theta0)

    # step 4: anti-clockwise rotation with phi0 along z-axis:   r0,t0,0 -> r0,t0,p0 
    (x,y,z) = rotate_z(x,y,z,phi0)

    # step 5: x,y,z -> r,t,p
    (r,t,p) = xyz2rtp(x,y,z)
    
    return (t,p)

def rotate_src_rec(ev_info,st_info,theta0,phi0,psi):
    ev_info_rotate = {}
    st_info_rotate = {}

    # rotate earthquakes
    for key_ev in ev_info:
        ev = ev_info[key_ev]
        ev_lat = np.array([ev.lat]); ev_lon = np.array([ev.lon])
        (ev_lat,ev_lon,) = rtp_rotation(ev_lat,ev_lon,theta0,phi0,psi)
        ev.lat = ev_lat[0]; ev.lon = ev_lon[0]
        ev_info_rotate[key_ev] = ev

    # rotate stations
    for key_st in st_info:
        st = st_info[key_st]
        st_lat = np.array([st.lat]); st_lon = np.array([st.lon])
        (st_lat,st_lon) = rtp_rotation(st_lat,st_lon,theta0,phi0,psi)
        st.lat = st_lat[0]; st.lon = st_lon[0]
        st_info_rotate[key_st] = st
    
    return (ev_info_rotate,st_info_rotate)

def rotate_src_rec_reverse(ev_info_rotate,st_info_rotate,theta0,phi0,psi):
    ev_info = {}
    st_info = {}

    # rotate earthquakes
    for key_ev in ev_info_rotate:
        ev = ev_info_rotate[key_ev]
        ev_lat = np.array([ev.lat]); ev_lon = np.array([ev.lon])
        (ev_lat,ev_lon,) = rtp_rotation_reverse(ev_lat,ev_lon,theta0,phi0,psi)
        ev.lat = ev_lat[0]; ev.lon = ev_lon[0]
        ev_info[key_ev] = ev

    # rotate stations
    for key_st in st_info_rotate:
        st = st_info_rotate[key_st]
        st_lat = np.array([st.lat]); st_lon = np.array([st.lon])
        (st_lat,st_lon) = rtp_rotation_reverse(st_lat,st_lon,theta0,phi0,psi)
        st.lat = st_lat[0]; st.lon = st_lon[0]
        st_info[key_st] = st
    
    return (ev_info,st_info)

# %%
# # 函数: 坐标旋转 rotate_src_rec(ev_info,st_info,theta0,phi0,psi): rotate to the new coordinate, satisfying the center r0,t0,p0 -> r0,0,0 and a anticlockwise angle psi
# # 满足中心点 r0,t0,p0 -> r0,0,0 以及 逆时针旋转角度 psi

# import numpy as np

# RAD2DEG = 180/np.pi
# DEG2RAD = np.pi/180
# R_earth = 6371.0

# # Spherical coordinates to Cartesian coordinate
# def rtp2xyz(r,theta,phi):
#     x = r * np.cos(theta*DEG2RAD) * np.cos(phi*DEG2RAD)
#     y = r * np.cos(theta*DEG2RAD) * np.sin(phi*DEG2RAD)
#     z = r * np.sin(theta*DEG2RAD)
#     return (x,y,z)

# # Cartesian coordinates to Spherical coordinate
# def xyz2rtp(x,y,z):
#     # theta: -90~90;  phi: -180~180
#     r       = np.sqrt(x**2+y**2+z**2)
#     theta   = np.arcsin(z/r)
#     phi     = np.arcsin(y/r/np.cos(theta))


#     if(phi > 0 and x*y < 0):
#         phi = np.pi - phi
#     if(phi < 0 and x*y > 0):
#         phi = -np.pi - phi
    
#     return (r,theta*RAD2DEG,phi*RAD2DEG)

# # anti-clockwise rotation along x-axis
# def rotate_x(x,y,z,theta):
#     new_x = x 
#     new_y = y *  np.cos(theta*DEG2RAD) + z * -np.sin(theta*DEG2RAD)
#     new_z = y *  np.sin(theta*DEG2RAD) + z *  np.cos(theta*DEG2RAD)
#     return (new_x,new_y,new_z)
    
# # anti-clockwise rotation along y-axis
# def rotate_y(x,y,z,theta):
#     new_x = x *  np.cos(theta*DEG2RAD) + z *  np.sin(theta*DEG2RAD)
#     new_y = y
#     new_z = x * -np.sin(theta*DEG2RAD) + z *  np.cos(theta*DEG2RAD)
#     return (new_x,new_y,new_z)

# # anti-clockwise rotation along z-axis
# def rotate_z(x,y,z,theta):
#     new_x = x *  np.cos(theta*DEG2RAD) + y * -np.sin(theta*DEG2RAD)
#     new_y = x *  np.sin(theta*DEG2RAD) + y *  np.cos(theta*DEG2RAD)
#     new_z = z 
#     return (new_x,new_y,new_z)

# # spherical Rotation

# # rotate to the new coordinate, satisfying the center r0,t0,p0 -> r0,0,0 and a anticlockwise angle psi
# def rtp_rotation(t,p,theta0,phi0,psi):
#     # step 1: r,t,p -> x,y,z
#     (x,y,z) = rtp2xyz(1.0,t,p)

#     # step 2: anti-clockwise rotation with -phi0 along z-axis:   r0,t0,p0 -> r0,t0,0
#     (x,y,z) = rotate_z(x,y,z,-phi0)

#     # step 3: anti-clockwise rotation with theta0 along y-axis:  r0,t0,0 -> r0,0,0
#     (x,y,z) = rotate_y(x,y,z,theta0)

#     # # step 4: anti-clockwise rotation with psi along x-axis
#     (x,y,z) = rotate_x(x,y,z,psi)

#     # step 5: x,y,z -> r,t,p
#     (new_r,new_t,new_p) = xyz2rtp(x,y,z)
    
#     return (new_t,new_p)


# def rtp_rotation_reverse(new_t,new_p,theta0,phi0,psi):
#     # step 1: r,t,p -> x,y,z
#     (x,y,z) = rtp2xyz(1.0,new_t,new_p)

#     # step 2: anti-clockwise rotation with -psi along x-axis
#     (x,y,z) = rotate_x(x,y,z,-psi)

#     # step 3: anti-clockwise rotation with -theta0 along y-axis:  r0,0,0 -> r0,t0,0 
#     (x,y,z) = rotate_y(x,y,z,-theta0)

#     # step 4: anti-clockwise rotation with phi0 along z-axis:   r0,t0,0 -> r0,t0,p0 
#     (x,y,z) = rotate_z(x,y,z,phi0)

#     # step 5: x,y,z -> r,t,p
#     (r,t,p) = xyz2rtp(x,y,z)
    
#     return (t,p)

# def rotate_src_rec(ev_info,st_info,theta0,phi0,psi):
#     ev_info_rotate = {}
#     st_info_rotate = {}

#     # rotate earthquakes
#     for key_ev in ev_info:
#         ev = ev_info[key_ev]
#         (ev.lat,ev.lon,) = rtp_rotation(ev.lat,ev.lon,theta0,phi0,psi)
#         ev_info_rotate[key_ev] = ev

#     # rotate stations
#     for key_st in st_info:
#         st = st_info[key_st]
#         (st.lat,st.lon) = rtp_rotation(st.lat,st.lon,theta0,phi0,psi)
#         st_info_rotate[key_st] = st
    
#     return (ev_info_rotate,st_info_rotate)

# def rotate_src_rec_reverse(ev_info_rotate,st_info_rotate,theta0,phi0,psi):
#     ev_info = {}
#     st_info = {}

#     # rotate earthquakes
#     for key_ev in ev_info_rotate:
#         ev = ev_info_rotate[key_ev]
#         (ev.lat,ev.lon) = rtp_rotation_reverse(ev.lat,ev.lon,theta0,phi0,psi)
#         ev_info[key_ev] = ev

#     # rotate stations
#     for key_st in st_info_rotate:
#         st = st_info_rotate[key_st]
#         (st.lat,st.lon) = rtp_rotation_reverse(st.lat,st.lon,theta0,phi0,psi)
#         st_info[key_st] = st
    
#     return (ev_info,st_info)

# %%
# 函数：最小二乘线性回归 linear_regression(X,Y)
def linear_regression(X,Y):
    slope,intercept = np.polyfit(X,Y,deg=1)
    fitted_values = slope * X + intercept
    residual = Y - fitted_values
    SEE = np.std(residual)
    return (slope,intercept,SEE)

# %% [markdown]
# 函数类：从地震数据 ev_info 和台站数据 st_info 中获得指定信息
# 
# Functions: obtain target information from ev_info and st_info

# %%
# 函数：data_lon_lat_dep_wt_ev(ev_info) 输出地震的 [lon,lat,dep,weight]. function: output the [lon,lat,dep,weight] of the earthquake
def data_lon_lat_dep_wt_ev(ev_info):
    lat = []
    lon = []
    dep = []
    weight = []
    for key in ev_info:
        lat.append(ev_info[key].lat)
        lon.append(ev_info[key].lon)
        dep.append(ev_info[key].dep)
        try:
            weight.append(ev_info[key].tag["weight"])
        except:
            weight.append(1.0)
    return [np.array(lon),np.array(lat),np.array(dep),np.array(weight)]

# %%
# 函数：data_ev_loc(ev_info) 输出地震的 [lon, lat, dep, ortime]. function: output the [lon, lat, dep, ortime] of the earthquake
def data_ev_loc(ev_info):
    lat = []
    lon = []
    dep = []
    ortime = []
    for key in ev_info:
        lat.append(ev_info[key].lat)
        lon.append(ev_info[key].lon)
        dep.append(ev_info[key].dep)
        ortime.append(ev_info[key].ortime.timestamp)
    return [np.array(lon),np.array(lat),np.array(dep),np.array(ortime)]

# %%
# 函数: data_lon_lat_ele_wt_st(ev_info,st_info) 输出台站的 [lon,lat,dep,weight]. function: output the [lon,lat,dep,weight] of the station
def data_lon_lat_ele_wt_st(ev_info,st_info):
    names = {}
    lat = [] 
    lon = [] 
    ele = []
    weight  = []
    for key_ev in ev_info:
        for key_t in ev_info[key_ev].t:     # absolute traveltime data 绝对走时数据
            name_st = ev_info[key_ev].t[key_t][0]
            names[name_st] = name_st

        for key_t in ev_info[key_ev].cs_dt: # common source differential traveltime data 共源双差走时数据
            name_st = ev_info[key_ev].cs_dt[key_t][0]
            names[name_st] = name_st
            name_st = ev_info[key_ev].cs_dt[key_t][1]
            names[name_st] = name_st

        for key_t in ev_info[key_ev].cr_dt: # common receiver differential traveltime data 共台站双差走时数据
            name_st = ev_info[key_ev].cr_dt[key_t][0]
            names[name_st] = name_st
            
    for name in names:  # only output the station which has data 只输出有数据的台站
        lat.append(st_info[name].lat)
        lon.append(st_info[name].lon)
        ele.append(st_info[name].ele)
        try:
            weight.append(st_info[name].tag["weight"])
        except:
            weight.append(1.0)
    return [np.array(lon),np.array(lat),np.array(ele),np.array(weight)]

# %%
# 函数：data_dis_time(ev_info) 输出 [震中距,到时]. function: output the [dis,time] of all data
def data_dis_time(ev_info,st_info):
    all_dis = []
    all_time = []
    for key_ev in ev_info:
        lat_ev = ev_info[key_ev].lat
        lon_ev = ev_info[key_ev].lon
        dep_ev = ev_info[key_ev].dep
        for key_t in ev_info[key_ev].t:
            all_time.append(ev_info[key_ev].t[key_t][2])
            lat_st = st_info[ev_info[key_ev].t[key_t][0]].lat
            lon_st = st_info[ev_info[key_ev].t[key_t][0]].lon
            ele_st = st_info[ev_info[key_ev].t[key_t][0]].ele
            dis = math.sqrt(cal_dis(lat_ev,lon_ev,lat_st,lon_st)**2 + (dep_ev+ele_st/1000)**2)
            all_dis.append(dis)

    return [np.array(all_dis),np.array(all_time)]

# %%
# 函数：data_dis_time(ev_info) 输出 [双差到时]. function: output the [cs_dt] of all data
def data_cs_dt(ev_info):
    all_time = []
    for key_ev in ev_info:
        for key_dt in ev_info[key_ev].cs_dt:
            all_time.append(ev_info[key_ev].cs_dt[key_dt][3])

    return np.array(all_time)

# %%
# 函数: data_dis_time_phase(ev_info,st_info,phase_list) 给定震相列表，分别输出每个震相对应数据的 [震中距,到时] # output (dis,time) of each given phase
def data_dis_time_phase(ev_info,st_info,phase_list):
    all_dis  = {}
    all_time = {}
    for phase in phase_list:
        all_dis[phase] = []
        all_time[phase] = []
    
    for key_ev in ev_info:
        lat_ev = ev_info[key_ev].lat
        lon_ev = ev_info[key_ev].lon
        dep_ev = ev_info[key_ev].dep
        for key_t in ev_info[key_ev].t:
            phase = key_t.split("+")[1]
            if (not phase in phase_list):
                continue

            all_time[phase].append(ev_info[key_ev].t[key_t][2])
            lat_st = st_info[ev_info[key_ev].t[key_t][0]].lat
            lon_st = st_info[ev_info[key_ev].t[key_t][0]].lon
            ele_st = st_info[ev_info[key_ev].t[key_t][0]].ele

            dis = math.sqrt(cal_dis(lat_ev,lon_ev,lat_st,lon_st)**2 + (dep_ev+ele_st/1000)**2)
            all_dis[phase].append(dis)

    for phase in phase_list:
        all_dis[phase] = np.array(all_dis[phase])
        all_time[phase] = np.array(all_time[phase])
        
    return [all_dis,all_time]

# %%
# 函数：data_lon_lat_dep_wt_ev(ev_info) 输出数据连线的 [line_x,line_y] # output the lines connecting station and earthquake for traveltime data
def data_line(ev_info,st_info):
    line_x = []
    line_y = []

    for key_ev in ev_info:
        lat_ev = ev_info[key_ev].lat
        lon_ev = ev_info[key_ev].lon
        for key_t in ev_info[key_ev].t:
            lat_st = st_info[ev_info[key_ev].t[key_t][0]].lat
            lon_st = st_info[ev_info[key_ev].t[key_t][0]].lon

            line_x.append([lon_ev,lon_st])
            line_y.append([lat_ev,lat_st])

    return [line_x,line_y]

# %% [markdown]
# 函数类：使用一些筛选原则，删除 ev_info 和 st_info 里面的部分数据
# 
# Functions: discard some data in ev_info and st_info based on selection criteria

# %%
# 函数：limit_ev_region(ev_info,lat1,lat2,lon1,lon2,dep1,dep2) 地震限制在研究范围内，删除地震. function: delete the earthquake out of the region
def limit_ev_region(ev_info,lat_min,lat_max,lon_min,lon_max,dep_min,dep_max):
    count_delete = 0

    del_key_ev = []
    for key_ev in ev_info:
        ev = ev_info[key_ev]
        lat = ev.lat
        lon = ev.lon
        dep = ev.dep

        if (lat < min(lat_min,lat_max) or lat > max(lat_min,lat_max)    \
            or lon < min(lon_min,lon_max) or lon > max(lon_min,lon_max)  \
            or dep < min(dep_min,dep_max) or dep > max(dep_min,dep_max)):
            del_key_ev.append(key_ev)
            count_delete += 1

        del_key_t = []
        for key_t in ev_info[key_ev].cr_dt:
            name_ev2 = ev_info[key_ev].cr_dt[key_t][1]
            lat2 = ev_info[name_ev2].lat
            lon2 = ev_info[name_ev2].lon
            dep2 = ev_info[name_ev2].dep
            if (lat2 < min(lat_min,lat_max) or lat2 > max(lat_min,lat_max)    \
                or lon2 < min(lon_min,lon_max) or lon2 > max(lon_min,lon_max) \
                or dep2 < min(dep_min,dep_max) or dep2 > max(dep_min,dep_max)):

                del_key_t.append(key_t)

        for key_t in del_key_t:
            del ev_info[key_ev].cr_dt[key_t]

        ev_info[key_ev].Ncr_dt = len(ev_info[key_ev].cr_dt)

    for key_ev in del_key_ev:
        del ev_info[key_ev]

    print("delete %d events out of the region, now %d earthquakes are retained within the study region"%(count_delete,len(ev_info)))
    return ev_info

# %%
# 函数: limit_st_region(ev_info,st_info,lat1,lat2,lon1,lon2) 台站限制在研究范围内，删除台站. function: delete the station out of the region
def limit_st_region(ev_info,st_info,lat1,lat2,lon1,lon2):

    for key_ev in ev_info:
        # delete the station out of the region in the absolute traveltime data
        del_key_t = []
        for key_t in ev_info[key_ev].t:
            name_st = ev_info[key_ev].t[key_t][0]
            lat_st = st_info[name_st].lat
            lon_st = st_info[name_st].lon
            if(lat_st < min(lat1,lat2) or lat_st > max(lat1,lat2) or lon_st < min(lon1,lon2) or lon_st > max(lon1,lon2)):
                del_key_t.append(key_t)
                
        for key_t in del_key_t:
            del ev_info[key_ev].t[key_t]
        ev_info[key_ev].Nt = len(ev_info[key_ev].t)

        # delete the station out of the region in the common source differential traveltime data
        del_key_t = []
        for key_t in ev_info[key_ev].cs_dt:
            name_st1 = ev_info[key_ev].cs_dt[key_t][0]
            lat_st1 = st_info[name_st1].lat
            lon_st1 = st_info[name_st1].lon

            name_st2 = ev_info[key_ev].cs_dt[key_t][1]
            lat_st2 = st_info[name_st2].lat
            lon_st2 = st_info[name_st2].lon
            if(lat_st1 < min(lat1,lat2) or lat_st1 > max(lat1,lat2) or lon_st1 < min(lon1,lon2) or lon_st1 > max(lon1,lon2) \
                or lat_st2 < min(lat1,lat2) or lat_st2 > max(lat1,lat2) or lon_st2 < min(lon1,lon2) or lon_st2 > max(lon1,lon2)):
                del_key_t.append(key_t)

        for key_t in del_key_t:
            del ev_info[key_ev].cs_dt[key_t]
        ev_info[key_ev].Ncs_dt = len(ev_info[key_ev].cs_dt)

        # delete the station out of the region in the common receiver differential traveltime data
        del_key_st = []
        for key_t in ev_info[key_ev].cr_dt:
            name_st = ev_info[key_ev].cr_dt[key_t][0]
            lat_st = st_info[name_st].lat
            lon_st = st_info[name_st].lon
            if(lat_st < min(lat1,lat2) or lat_st > max(lat1,lat2) or lon_st < min(lon1,lon2) or lon_st > max(lon1,lon2)):
                del_key_st.append(key_t)

        for key_t in del_key_st:
            del ev_info[key_ev].cr_dt[key_t]
        ev_info[key_ev].Ncr_dt = len(ev_info[key_ev].cr_dt)

    return ev_info



# %%
# 函数：limit_epi_dis(ev_info,st_info,epi_dis1,epi_dis2) 删除震中距范围在 dis1 - dis2 之间的数据. function: delete the station with epicentral distance from epi_dis1 to epi_dis2
def limit_epi_dis(ev_info,st_info,epi_dis1,epi_dis2):

    for key_ev in ev_info:
        ev = ev_info[key_ev]
        
        lat_ev = ev.lat
        lon_ev = ev.lon

        # delete the absolute traveltime data
        del_key_t = []
        for key_t in ev.t:
            stname = ev.t[key_t][0]
            lat_st = st_info[stname].lat
            lon_st = st_info[stname].lon
            dis = cal_dis(lat_ev, lon_ev, lat_st, lon_st)
            if (dis > epi_dis1 and dis < epi_dis2):
                del_key_t.append(key_t)
        for key_t in del_key_t:
            del ev.t[key_t]
        ev.Nt = len(ev.t)

        # delete the common source differential traveltime data
        del_key_t = []
        for key_t in ev.cs_dt:
            for i in range(2):
                stname = ev.t[key_t][i]
                lat_st = st_info[stname].lat
                lon_st = st_info[stname].lon
                dis = cal_dis(lat_ev, lon_ev, lat_st, lon_st)
            
                if (dis > epi_dis1 and dis < epi_dis2):
                    del_key_t.append(key_t)
                    break
        for key_t in del_key_t:
            del ev.cs_dt[key_t]
        ev.Ncs_dt = len(ev.cs_dt)

        # delete the common receiver differential traveltime data
        del_key_t = []
        for key_t in ev.cr_dt:
            stname = ev.cr_dt[key_t][0]
            lat_st = st_info[stname].lat
            lon_st = st_info[stname].lon
            dis = cal_dis(lat_ev, lon_ev, lat_st, lon_st)
            if (dis > epi_dis1 and dis < epi_dis2):
                del_key_t.append(key_t)

            lat_ev2 = ev_info[ev.cr_dt[key_t][1]].lat
            lon_ev2 = ev_info[ev.cr_dt[key_t][1]].lon
            dis = cal_dis(lat_ev2, lon_ev2, lat_st, lon_st)
            if (dis > epi_dis1 and dis < epi_dis2):
                del_key_t.append(key_t)
        
        for key_t in del_key_t:
            del ev.cr_dt[key_t]
        ev.Ncr_dt = len(ev.cr_dt)
        
        
        ev_info[key_ev] = ev

    return ev_info

# %%
# 函数: limit_data_residual(ev_info,st_info,slope,intercept,up,down)  将数据限制在 直线 time = dis * slope + intercept 的两侧， up and down 的范围里面
# remove outliers, only retain data satisfying:     slope * dis + intercept + down < time < slope * dis + intercept + up
def limit_data_residual(ev_info,st_info,slope,intercept,up,down):
    for key_ev in ev_info:
        lat_ev = ev_info[key_ev].lat
        lon_ev = ev_info[key_ev].lon
        dep_ev = ev_info[key_ev].dep
        del_key_t = []
        for key_t in ev_info[key_ev].t:
            name_st = ev_info[key_ev].t[key_t][0]
            lat_st = st_info[name_st].lat
            lon_st = st_info[name_st].lon
            ele_st = st_info[name_st].ele
            dis = math.sqrt(cal_dis(lat_ev,lon_ev,lat_st,lon_st)**2 + (dep_ev+ele_st/1000)**2)
            residual = ev_info[key_ev].t[key_t][2] - (slope*dis+intercept)

            if (residual < down or residual > up):
                del_key_t.append(key_t)

        for key_t in del_key_t:
            del ev_info[key_ev].t[key_t]

    for key_ev in ev_info:
        ev_info[key_ev].Nt = len(ev_info[key_ev].t)

    return ev_info



# %%
# 函数：limit_data_phase(ev_info,phase_list)   仅保留指定震相
def limit_data_phase(ev_info,phase_list):
    for key_ev in ev_info:
        # process the absolute traveltime data
        new_t = {}
        for key_t in ev_info[key_ev].t:
            phase = ev_info[key_ev].t[key_t][1]
            if phase in phase_list:
                new_t[key_t] = ev_info[key_ev].t[key_t]
        
        ev_info[key_ev].t = new_t
        ev_info[key_ev].Nt = len(ev_info[key_ev].t)

        # process the common source differential traveltime data
        new_t = {}
        for key_t in ev_info[key_ev].cs_dt:
            phase = ev_info[key_ev].cs_dt[key_t][2]
            phase = phase.split(",")[0]
            if phase in phase_list:
                new_t[key_t] = ev_info[key_ev].cs_dt[key_t]

        ev_info[key_ev].cs_dt = new_t
        ev_info[key_ev].Ncs_dt = len(ev_info[key_ev].cs_dt)

        # process the common receiver differential traveltime data
        new_t = {}
        for key_t in ev_info[key_ev].cr_dt:
            phase = ev_info[key_ev].cr_dt[key_t][2]
            phase = phase.split(",")[0]
            if phase in phase_list:
                new_t[key_t] = ev_info[key_ev].cr_dt[key_t]

        ev_info[key_ev].cr_dt = new_t
        ev_info[key_ev].Ncr_dt = len(ev_info[key_ev].cr_dt)

    return ev_info

# %%
# 函数：limit_min_Nt(min_Nt_thd, ev_info) 删除到时数量小于阈值的地震 function: delete the earthquake with the number of data less than min_Nt_thd
def limit_min_Nt(min_Nt_thd, ev_info):
    Nev = len(ev_info)

    del_key_ev = []
    for key_ev in ev_info:
        if(ev_info[key_ev].Nt < min_Nt_thd):
            del_key_ev.append(key_ev)
    
    for key_ev in del_key_ev:
        del ev_info[key_ev]

    print("Original data set has %d earthquakes, %d earthquakes are deleted, %d earthquakes are retained"%(Nev,len(del_key_ev),len(ev_info)))

    return ev_info

# %%
# 函数：limit_azi_gap(gap_thd) # 计算所有地震的方位角gap，并删除 gap > gap_thd 的地震. function: calculate azimuthal gap for all events and delete events with gap > gap_thd. 
def limit_azi_gap(gap_thd,ev_info,st_info):
    Nev = len(ev_info)

    del_key_ev = []
    for key_ev in ev_info:
        ev = ev_info[key_ev]
        gap = cal_azi_gap(ev,st_info)
        if (gap > gap_thd):
            del_key_ev.append(key_ev)
        else:
            ev_info[key_ev].tag["azi_gap"] = gap
    for key_ev in del_key_ev:
        del ev_info[key_ev]

    print("Original data set has %d earthquakes, %d earthquakes are deleted, %d earthquakes are retained"%(Nev,len(del_key_ev),len(ev_info)))

    return ev_info

# 函数：cal_azi_gap(ev,st_info) 计算单个地震的方位角空缺 calculate the azimuthal gap of the earthquake
def cal_azi_gap(ev,st_info):
    azi_all = []
    lat_ev = ev.lat
    lon_ev = ev.lon
    stlist = {}
    for key in ev.t:
        stname = ev.t[key][0]
        if (not stname in stlist):
            lat_st = st_info[stname].lat
            lon_st = st_info[stname].lon
            azi = cal_azimuth(lat_ev, lon_ev, lat_st, lon_st)
            azi_all.append(azi)
            stlist[stname] = 1

    azi_all.sort()
    if(len(azi_all) < 2):
        return 360.0
    else:
        gap = 0.0
        for i in range(len(azi_all)-1):
            gap = max(gap,azi_all[i+1] - azi_all[i])
        gap = max(gap,azi_all[0] + 360 - azi_all[-1])
        return gap



# %%
# 函数：limit_earthquake_decluster_Nt(ev_info,dlat,dlon,ddep,Top_N) 将区域划分为若干个subdomain，按照到时数量排序，每个box内仅保留到时数量最多的Top_N个地震，
# option 3, declustering. Divide the region into several subdomains, retain the Top N earthquakes in terms of the number of arrival times in each subdomain.
def limit_earthquake_decluster_Nt(ev_info,dlat,dlon,ddep,Top_N):
    # catagrate earthquakes into different subdomains 细分地震
    [ev_info,tag2name] = tag_event_cluster(dlat,dlon,ddep,ev_info)

    # sort earthquakes in the same subdomain 在每个tag中，根据接收数量，对地震质量进行排序
    tag2name = sort_cluster_Nt(ev_info, tag2name)

    # only retain Top_N earthquakes in each subdomain 在每个tag中，优先选择 Top_N 个地震
    [ev_info,tag2name] = limit_decluster(ev_info, tag2name,Top_N)

    return ev_info



# 函数： tag_event_cluster(size_lat,size_lon,size_dep,ev_info)细分研究区域，每个地震属于一个子区域，放在tag里面
def tag_event_cluster(size_lat,size_lon,size_dep,ev_info):
    tag2name = {}
    for key_ev in ev_info:
        name = ev_info[key_ev].name
        lat = ev_info[key_ev].lat
        lon = ev_info[key_ev].lon
        dep = ev_info[key_ev].dep
        tag = "%d_%d_%d"%(math.floor(lon/size_lon),math.floor(lat/size_lat),math.floor(dep/size_dep))
        ev_info[key_ev].tag["cluster"] = tag

        if (tag in tag2name):
            tag2name[tag].append(name)
        else:
            tag2name[tag] = []
            tag2name[tag].append(name)

    return [ev_info,tag2name]

# 函数： sort_cluster_Nt(ev_info, tag2name) 在每个tag中，根据接收数量，对地震质量进行排序
def sort_cluster_Nt(ev_info, tag2name):
    for key_tag in tag2name:
        names_ev = tag2name[key_tag]
        Nt = []
        for key_ev in names_ev:
            Nt.append(len(ev_info[key_ev].t))
        
        # 根据地震到时数量，给该tag下的地震排序
        sorted_Nt = sorted(enumerate(Nt), key=lambda x: x[1], reverse=True)
        tag2name[key_tag] = []
        for index, Nt in sorted_Nt:
            tag2name[key_tag].append(names_ev[index])
        
    return tag2name
        
# 函数： limit_cluster(ev_info, tag2name, Max) 在每个tag中，优先选择 threshold 个地震
def limit_decluster(ev_info, tag2name, Max):
    del_key_ev = []
    for key_tag in tag2name:
        names_ev = tag2name[key_tag]
        
        if(len(names_ev) > Max):
            tag2name[key_tag] = names_ev[0:Max]
            for i in range(Max,len(names_ev)):    # 删除排序超过 threshold 的地震
                del_key_ev.append(names_ev[i])
    
    for key_ev in del_key_ev:
        del ev_info[key_ev]
    
    return [ev_info,tag2name]



# %% [markdown]
# 函数类：给地震、台站、数据添加权重
# 
# Functions: assign weights to earthquakes, stations, and data

# %%
# 函数：box_weighting_ev (ev_info,dlat,dlon,ddep) 赋予地震权重. Assign box-weight to the earthquakes 
def box_weighting_ev(ev_info,dlon,dlat,ddep):

    # 归类
    distribute = {}
    all_tag_wt = {}

    for key_ev in ev_info:
        lat_id = math.floor((ev_info[key_ev].lat) / dlat)
        lon_id = math.floor((ev_info[key_ev].lon) / dlon)
        dep_id = math.floor((ev_info[key_ev].dep) / ddep)

        tag = '%d_%d_%d'%(lat_id,lon_id,dep_id)
        if (tag in distribute):
            distribute[tag] += 1
        else:
            distribute[tag] = 1

    max_weight = 0
    for tag in distribute:
        all_tag_wt[tag] = 1.0/math.sqrt(distribute[tag])
        max_weight = max(max_weight,all_tag_wt[tag])

    for key_ev in ev_info:
        lat_id = math.floor((ev_info[key_ev].lat) / dlat)
        lon_id = math.floor((ev_info[key_ev].lon) / dlon)
        dep_id = math.floor((ev_info[key_ev].dep) / ddep)
        
        tag = '%d_%d_%d'%(lat_id,lon_id,dep_id)
        
        ev_info[key_ev].tag["weight"] = all_tag_wt[tag]/max_weight

    return ev_info

# %%
# 函数：geographical_weighting_ev_rough(ev_info,dlat,dlon,ddep) 赋予地震权重 Assign geographical_weighting to the earthquakes roughly
def geographical_weighting_ev_rough(ev_info,dlat,dlon,ddep,coefficient = 0.5):
    # 归类
    distribute = {}
    all_tag_wt = {}

    for key_ev in ev_info:
        lat_id = int(ev_info[key_ev].lat/dlat)
        lon_id = int(ev_info[key_ev].lon/dlat)
        dep_id = int(ev_info[key_ev].dep/ddep)


        tag = '%d_%d_%d'%(lat_id,lon_id,dep_id)
        if (tag in distribute):
            distribute[tag] += 1
        else:
            distribute[tag] = 1
    
    # 计算 每一类的 weight
    delta0 = 0
    for tag1 in distribute:
        tmp1 = tag1.split('_')
        # evlat1 = float(tmp1[0])*dlat; evlon1 = float(tmp1[1])*dlon; evdep1 = float(tmp1[2])*ddep

        for tag2 in distribute:
            tmp2 = tag2.split('_')
            # evlat2 = float(tmp2[0])*dlat; evlon2 = float(tmp2[1])*dlon; evdep2 = float(tmp2[2])*ddep

            # distance of id
            delta_tp = math.sqrt((int(tmp1[0]) - int(tmp2[0]))**2 + (int(tmp1[1]) - int(tmp2[1]))**2 + (int(tmp1[2]) - int(tmp2[2]))**2)
            delta0 = delta0 + distribute[tag1] * distribute[tag2] * delta_tp
            
    delta0 = delta0/(len(ev_info)**2) * coefficient

    max_weight = 0.0
    for tag1 in distribute:
        tmp1 = tag1.split('_')
        # evlat1 = float(tmp1[0])*dlat; evlon1 = float(tmp1[1])*dlon; evdep1 = float(tmp1[2])*ddep

        weight = 0
        for tag2 in distribute:
            tmp2 = tag2.split('_')
            # evlat2 = float(tmp2[0])*dlat; evlon2 = float(tmp2[1])*dlon; evdep2 = float(tmp2[2])*ddep

            delta_tp = math.sqrt((int(tmp1[0]) - int(tmp2[0]))**2 + (int(tmp1[1]) - int(tmp2[1]))**2 + (int(tmp1[2]) - int(tmp2[2]))**2)
            weight = weight + math.exp(-(delta_tp/delta0)**2) * distribute[tag2]

        all_tag_wt[tag1] = (1.0/weight)
        max_weight = max(max_weight,1.0/weight)

    # 根据地震的tag，给每个地震添加权重
    for key_ev in ev_info:
        lat_id = int(ev_info[key_ev].lat/dlat)
        lon_id = int(ev_info[key_ev].lon/dlon)
        dep_id = int(ev_info[key_ev].dep/ddep)
        
        tag = '%d_%d_%d'%(lat_id,lon_id,dep_id)
        
        ev_info[key_ev].tag["weight"] = all_tag_wt[tag]/max_weight
    return ev_info


# %%
# 函数：box_weighting_st(ev_info,st_info,dlat,dlon) 赋予台站权重 Assign geographical_weighting to the stations roughly
def box_weighting_st(ev_info,st_info,dlon,dlat):

    [lon_ev,lat_ev,dep_ev,wt_ev] = data_lon_lat_dep_wt_ev(ev_info)

    # 整合所有涉及到的台站
    wt_st = {}
    name_st = {}
    for key_ev in ev_info:
        for key_t in ev_info[key_ev].t:
            name_rec = ev_info[key_ev].t[key_t][0]
            wt_st[name_rec] = -1.0
            name_st[name_rec] = 1

    # 归类
    distribute = {}
    all_tag_wt = {}

    # 统计每个小subdomain内的台站数量
    for key_st in name_st:
        lat_id = math.floor((st_info[key_st].lat) / dlat)
        lon_id = math.floor((st_info[key_st].lon) / dlon)

        tag = '%d_%d'%(lat_id,lon_id)
        if (tag in distribute):
            distribute[tag] += 1
        else:
            distribute[tag] = 1

    max_weight = 0
    for tag in distribute:
        all_tag_wt[tag] = 1.0/math.sqrt(distribute[tag])
        max_weight = max(max_weight,all_tag_wt[tag])

    # 赋予台站权重
    for key_st in name_st:
        lat_id = math.floor((st_info[key_st].lat) / dlat)
        lon_id = math.floor((st_info[key_st].lon) / dlon)
        tag = '%d_%d'%(lat_id,lon_id)       
        wt_st[key_st] = all_tag_wt[tag]/max_weight

    # modify weight tag in st_info
    for key_t in wt_st:
        st_info[key_t].tag["weight"] = wt_st[key_t]

    # modify weight of abs data ev_info
    for key_ev in ev_info:
        for key_t in ev_info[key_ev].t:
            name_rec = ev_info[key_ev].t[key_t][0]
            ev_info[key_ev].t[key_t][3] = wt_st[name_rec]    

    return [ev_info,st_info]


# %%
# 函数：geographical_weighting_st(ev_info,st_info,) 赋予台站权重 Assign geographical_weighting to the stations roughly
def geographical_weighting_st(ev_info,st_info,coefficient = 0.5):

    # 整合所有涉及到的台站
    wt_st = {}
    name_st = {}
    for key_ev in ev_info:
        for key_t in ev_info[key_ev].t:
            name_rec = ev_info[key_ev].t[key_t][0]
            wt_st[name_rec] = -1.0
            name_st[name_rec] = 1

    # 计算每个台站的weight
    delta0 = 0
    for key_st1 in name_st:
        stlat1 = st_info[key_st1].lat
        stlon1 = st_info[key_st1].lon

        for key_st2 in name_st:
            stlat2 = st_info[key_st2].lat
            stlon2 = st_info[key_st2].lon

            delta_tp = cal_dis(stlat1,stlon1,stlat2,stlon2)
            delta0 = delta0 + delta_tp
            
    delta0 = delta0/(len(wt_st)**2)*coefficient

    max_weight = 0.0
    for key_st1 in name_st:
        stlat1 = st_info[key_st1].lat
        stlon1 = st_info[key_st1].lon

        weight = 0
        for key_st2 in name_st:
            stlat2 = st_info[key_st2].lat
            stlon2 = st_info[key_st2].lon

            delta_tp = cal_dis(stlat1,stlon1,stlat2,stlon2)
            weight = weight + math.exp(-(delta_tp/delta0)**2)

        wt_st[key_st1] = (1.0/weight)
        max_weight = max(max_weight,1.0/weight)

    for key_st1 in wt_st:
        wt_st[key_st1] = wt_st[key_st1]/max_weight

    # 给每个地震中的数据添加weight
    for key_ev in ev_info:
        for key_t in ev_info[key_ev].t:
            name_rec = ev_info[key_ev].t[key_t][0]
            if (not name_rec in wt_st):
                ValueError("数据的台站不在计算列表中 station is excluded in the list")
            
            if (len(ev_info[key_ev].t[key_t])==3):
                ev_info[key_ev].t[key_t].append(wt_st[name_rec])
            elif (len(ev_info[key_ev].t[key_t])==4):
                ev_info[key_ev].t[key_t][3] = wt_st[name_rec]
            else:
                ValueError("数据的权重信息有误, error in the weight information of the absoulte traveltime")

    # modify weight tag in st_info
    for key_t in wt_st:
        st_info[key_t].tag["weight"] = wt_st[key_t]

    # modify weight of abs data ev_info
    for key_ev in ev_info:
        for key_t in ev_info[key_ev].t:
            name_rec = ev_info[key_ev].t[key_t][0]
            ev_info[key_ev].t[key_t][3] = wt_st[name_rec] 
            
    return [ev_info,st_info]


# %% [markdown]
# 函数类：数据添加噪声
# 
# Function: add noise into data

# %%
# 函数：assign_gaussian_noise():
def assign_gaussian_noise(ev_info,sigma):

    # 记录一下对应台站有哪些震相
    st2phase = {}   # 台站名字 -> [与这个台站有关的绝对到时数据的键值]

    
    for key_ev in ev_info:
        # 绝对到时噪声
        for key_t in ev_info[key_ev].t:
            stname = ev_info[key_ev].t[key_t][0]
            ev_info[key_ev].t[key_t][2] = ev_info[key_ev].t[key_t][2] + np.random.normal(0,sigma)
            if(stname in st2phase):
                st2phase[stname].append(key_t)
            else:
                st2phase[stname] = [key_t]
    
    for key_ev in ev_info:
        # 共源双差到时噪声
        for key_dt in ev_info[key_ev].cs_dt:
            stname1 = ev_info[key_ev].cs_dt[key_dt][0]
            stname2 = ev_info[key_ev].cs_dt[key_dt][1]
            t1 = -999
            t2 = -999
            # 寻找有没有该数据的到时
            if (stname1 in st2phase):
                for key_t in st2phase[stname1]:
                    if (key_t in ev_info[key_ev].t):
                        t1 = ev_info[key_ev].t[key_t][2]
                        break
            if (stname2 in st2phase):
                for key_t in st2phase[stname2]:
                    if (key_t in ev_info[key_ev].t):
                        t2 = ev_info[key_ev].t[key_t][2]
                        break
            
            if (t1 == -999 or t2 == -999):
                # 没有绝对到时数据，，共源双差数据残差增加 sqrt(2)倍的噪声
                ev_info[key_ev].cs_dt[key_dt][3] = ev_info[key_ev].cs_dt[key_dt][3] + np.random.normal(0,sigma*np.sqrt(2))
                print('no data: ', key_ev, key_dt)
            else:
                # 有绝对到时数据，共源双差数据由相减得的
                ev_info[key_ev].cs_dt[key_dt][3] = t1 - t2

        # 共台站双差到时
        for key_dt in ev_info[key_ev].cr_dt:
            stname  = ev_info[key_ev].cr_dt[key_dt][0]
            key_ev2 = ev_info[key_ev].cr_dt[key_dt][1]

            t1 = -999
            t2 = -999
            # 寻找有没有该数据的到时
            if (stname in st2phase):
                for key_t in st2phase[stname]:
                    if (key_t in ev_info[key_ev].t):
                        t1 = ev_info[key_ev].t[key_t][2]
                        break
                    else:
                        print('not found 1: ', key_ev, key_t)

                for key_t in st2phase[stname]:
                    if (key_t in ev_info[key_ev2].t):
                        t2 = ev_info[key_ev2].t[key_t][2]
                        break
                    else:
                        print('not found 2: ', key_ev, key_t)

            if (t1 == -999 or t2 == -999):
                # 没有绝对到时数据，，共源双差数据残差增加 sqrt(2)倍的噪声
                ev_info[key_ev].cr_dt[key_dt][3] = ev_info[key_ev].cr_dt[key_dt][3] + np.random.normal(0,sigma*np.sqrt(2))
                print('no data: ', key_ev, key_dt)
            else:
                # 有绝对到时数据，共源双差数据由相减得的
                ev_info[key_ev].cr_dt[key_dt][3] = t1 - t2

    return ev_info

# %%
# 函数：assign_uniform_noise_to_ev():
def assign_uniform_noise_to_ev(ev_info, range_lat, range_lon, range_dep, range_time):

    # 循环所有地震，给这些地震赋予噪声
    ev_noise = {}   # 地震名字 -> noise of [lat,lon,dep,ortime]
    # loop 地震列表
    for key_ev in ev_info:
        evname = key_ev
        if (evname in ev_noise):
            print("error: repeated earthquake name")
            exit()
        else:
            # 生成噪声
            ev_noise[evname] = np.random.uniform(-1,1,4) * np.array([range_lat,range_lon,range_dep,range_time])

    # 给每个 数据 添加噪声
    for key_ev in ev_info:

        # 绝对到时噪声
        for key_t in ev_info[key_ev].t:

            ev_info[key_ev].t[key_t][2] = ev_info[key_ev].t[key_t][2] - ev_noise[key_ev][3]

    
        # 共源双差到时噪声 (双差到时没有变化)

        # 共台站双差到时
        for key_dt in ev_info[key_ev].cr_dt:
            key_ev2 = ev_info[key_ev].cr_dt[key_dt][1]

            if (key_ev2 in ev_noise):
                ev_info[key_ev].cr_dt[key_dt][3] = ev_info[key_ev].cr_dt[key_dt][3] - ev_noise[key_ev][3] + ev_noise[key_ev2][3]
            else:
                print("earthquake %s is not included in ev_list"%(key_ev2))
                ev_noise[key_ev2] = np.random.uniform(-1,1,4) * np.array([range_lat,range_lon,range_dep,range_time])
                ev_info[key_ev].cr_dt[key_dt][3] = ev_info[key_ev].cr_dt[key_dt][3] - ev_noise[key_ev][3] + ev_noise[key_ev2][3]


    # 给每个地震添加噪声
    for key_ev in ev_noise:
        ev_info[key_ev].lat = ev_info[key_ev].lat + ev_noise[key_ev][0]
        ev_info[key_ev].lon = ev_info[key_ev].lon + ev_noise[key_ev][1]
        ev_info[key_ev].dep = abs(ev_info[key_ev].dep + ev_noise[key_ev][2])
        ev_info[key_ev].ortime = ev_info[key_ev].ortime + ev_noise[key_ev][3]


    return ev_info

# %% [markdown]
# 函数类：生成差分数据
# 
# Functions: generate differential traveltime

# %%
# 函数：generate_cs_dif(ev_info,st_info,dis_thd,azi_thd): 从绝对到时生成共源双差到时, 台间距小于 dis_thd, 方位差小于 azi_thd
# function: generate common source differential traveltime data from absolute traveltime data, the stations separation is less than dis_thd, the azimuth difference is less than azi_thd
def generate_cs_dif(ev_info,st_info,dis_thd,azi_thd):
    count_t = 0
    count_cs_dt = 0

    for key_ev in ev_info:
        ev = ev_info[key_ev]

        lat_ev = ev.lat
        lon_ev = ev.lon
        
        # 遍历所有的到时
        name_st_list = []   # names of stations
        t_list = []         # traveltime
        wt_list = []        # weight
        for key_t in ev.t:
            name_st_list.append(ev.t[key_t][0])
            t_list.append(ev.t[key_t][2])
            wt_list.append(ev.t[key_t][3])
            count_t += 1
            
        # 寻找可能的双差到时
        for id_st1 in range(len(name_st_list)-1):
            name_st1 = name_st_list[id_st1]
            lat_st1  = st_info[name_st1].lat
            lon_st1  = st_info[name_st1].lon
            t_st1    = t_list[id_st1]
            wt_st1   = wt_list[id_st1]

            for id_st2 in range(id_st1+1,len(name_st_list)):
                name_st2 = name_st_list[id_st2]
                lat_st2  = st_info[name_st2].lat
                lon_st2  = st_info[name_st2].lon
                t_st2    = t_list[id_st2]
                wt_st2   = wt_list[id_st2]

                dis = cal_dis(lat_st1,lon_st1,lat_st2,lon_st2)
                azi_st1 = cal_azimuth(lat_ev,lon_ev,lat_st1,lon_st1)
                azi_st2 = cal_azimuth(lat_ev,lon_ev,lat_st2,lon_st2)
                
                azi_dif = abs(azi_st1 - azi_st2)

                if(dis < dis_thd and (azi_dif < azi_thd or (360-azi_dif) < azi_thd )):
                    ev.cs_dt["%s+%s+%s"%(name_st1,name_st2,"P,cs")] = [name_st1,name_st2,"P,cs",t_st1-t_st2,(wt_st1+wt_st2)/2]
                    count_cs_dt += 1

        ev_info[key_ev].Ncs_dt = len(ev.cs_dt)

    print('we generate %d common source differential traveltimes from %s absolute traveltimes'%(count_cs_dt,count_t))
    return ev_info


# %%
# 函数：generate_cr_dif(ev_info,st_info,dis_thd): 从绝对到时生成共台双差到时, 震源间距小于 dis_thd
# function: generate common receiver differential traveltime data from absolute traveltime data, the earthquake separation is less than dis_thd
def generate_cr_dif(ev_info,st_info,dis_thd,azi_thd):
    
    # 构造映射：rec2src[name_ev] -> {name_st: [name_ev, name_st, t, wt]; name_st: [name_ev, name_st, t, wt]; ...}
    rec2src = build_rec_src_map(ev_info,dis_thd)
    print("rec to src map generation finished")

    # 构造双差数据关联映射：rec2src_pair[key_t]
    rec2src_pair = build_rec_src_pair_map(rec2src)
    print("rec to src_pair map generation finished")

    for key_t in rec2src_pair:
        name_st = key_t.split('+')[0]
        lat_st = st_info[name_st].lat
        lon_st = st_info[name_st].lon

        for ev_tag in rec2src_pair[key_t]:
            name_ev1 = rec2src_pair[key_t][ev_tag][0]
            lat_ev1  = ev_info[name_ev1].lat
            lon_ev1  = ev_info[name_ev1].lon
            dep_ev1  = ev_info[name_ev1].dep

            name_ev2 = rec2src_pair[key_t][ev_tag][1]
            lat_ev2  = ev_info[name_ev2].lat
            lon_ev2  = ev_info[name_ev2].lon
            dep_ev2  = ev_info[name_ev2].dep

            dis_xy = cal_dis(lat_ev1,lon_ev1,lat_ev2,lon_ev2)
            dis_z  = abs(dep_ev1 - dep_ev2)
            dis = math.sqrt(dis_xy**2 + dis_z**2)
            if(dis > dis_thd):      # limit of the distance between two earthquakes
                continue

            azi1 = cal_azimuth(lat_ev1,lon_ev1,lat_st,lon_st)
            azi2 = cal_azimuth(lat_ev2,lon_ev2,lat_st,lon_st)
            azi_dif = abs(azi1 - azi2)

            if(azi_dif > azi_thd and (360-azi_dif) > azi_thd):     # limit of the azimuth difference between two earthquakes
                continue

            t_ev1 = ev_info[name_ev1].t[key_t][2]
            t_ev2 = ev_info[name_ev2].t[key_t][2]
            wt_ev1 = ev_info[name_ev1].t[key_t][3] * ev_info[name_ev1].tag["weight"]
            wt_ev2 = ev_info[name_ev2].t[key_t][3] * ev_info[name_ev2].tag["weight"]
            wt = (wt_ev1 + wt_ev2)/2/ev_info[name_ev1].tag["weight"]        # 实际的数据权重是 wt_ev1 + wt_ev2，但是在TomoATT计算中，我们需要将其除以 ev_info[name_ev1].tag["weight"]
            
            ev_info[name_ev1].cr_dt["%s+%s+%s"%(name_st,name_ev2,"P,cr")] = [name_st,name_ev2,"P,cr",t_ev1-t_ev2,wt]

    # 统计双差数量
    count_cr_dt = 0
    count_t     = 0
    for key_ev in ev_info:
        ev_info[key_ev].Ncr_dt = len(ev_info[key_ev].cr_dt)
        count_cr_dt += ev_info[key_ev].Ncr_dt
        count_t     += ev_info[key_ev].Nt

    print('we generate %d common receiver differential traveltimes from %s absolute traveltimes'%(count_cr_dt,count_t))

    return ev_info


# 构造映射：rec2src      = {key_t: dict_tag; key_t: dict_tag; ...} 
#          dict_tag     = {tag: list_name_ev; tag:list_name_ev; ...}
#          list_name_ev = [name_ev1, name_ev2, ...]
# 根据地震位置，将其分配进入不同的子区域。子区域大小为 dlat*dlon*ddep. 在做共台站双差的时候，仅在同一个子区域或相邻子区域内的地震对才会被考虑
def build_rec_src_map(ev_info,dis_thd):
    rec2src = {}
    for key_ev in ev_info:
        name_ev = ev_info[key_ev].name
        lat = ev_info[key_ev].lat
        lon = ev_info[key_ev].lon
        dep = ev_info[key_ev].dep
        tag_dep = math.floor(dep/dis_thd)
        tag_lat = math.floor(lat/180*math.pi*R_earth/dis_thd)
        tag_lon = math.floor(lon/180*math.pi*R_earth*math.cos(lat)/dis_thd)
        tag = "%d_%d_%d"%(tag_lon,tag_lat,tag_dep)


        for key_t in ev_info[key_ev].t:

            # 创建 dictionary
            if (not key_t in rec2src):
                rec2src[key_t] = {tag:[]}
            elif (not tag in rec2src[key_t]):
                rec2src[key_t][tag] = []
            
            # 添加数据
            rec2src[key_t][tag].append(name_ev)

    return rec2src

# 函数：generate_adjacent_tag(tag) 生成 tag 周围的 tag
def generate_adjacent_tag(tag): # 不包含 tag 自己
    adjacent_tag_list = []
    tmp = tag.split('_')
    tag_lon = int(tmp[0])
    tag_lat = int(tmp[1])
    tag_dep = int(tmp[2])

    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                if(i == 0 and j == 0 and k == 0):
                    continue
                adjacent_tag_list.append("%d_%d_%d"%(tag_lon+i,tag_lat+j,tag_dep+k))

    return adjacent_tag_list


# 构造映射：rec2src_pair
def build_rec_src_pair_map(rec2src):
    rec2src_pair = {}
    
    for key_t in rec2src:
        rec2src_pair[key_t] = {}

        for tag in rec2src[key_t]:
            name_ev_list1 = rec2src[key_t][tag]

            name_ev_list2 = rec2src[key_t][tag]
            adjacent_tag_list = generate_adjacent_tag(tag)
            for adjacent_tag in adjacent_tag_list:
                if (adjacent_tag in rec2src[key_t]):      # 如果周围 tag 的区域 有地震，就加入地震列表
                    name_ev_list2 = name_ev_list2 + rec2src[key_t][adjacent_tag]

            # 寻找可能的地震对
            for id_ev1 in range(len(name_ev_list1)-1):
                name_ev1 = name_ev_list1[id_ev1]

                for id_ev2 in range(id_ev1+1,len(name_ev_list2)):   # 从 id_ev1+1 开始，就已经排除了 tag 内部地震的重复
                    name_ev2 = name_ev_list2[id_ev2]

                    ev_tag1 = "%s+%s"%(name_ev1,name_ev2)
                    ev_tag2 = "%s+%s"%(name_ev2,name_ev1)

                    if(ev_tag1 in rec2src_pair[key_t] or ev_tag2 in rec2src_pair[key_t]):
                        continue

                    rec2src_pair[key_t][ev_tag1] = [name_ev1,name_ev2]        


    return rec2src_pair

# %% [markdown]
# 函数类：读取和输出数据文件
# 
# Functions: read and write src_rec.dat file

# %%
# 函数：reorder_src(ev) 按照顺序给地震重新编号. 如果地震没有数据，则编号为-999. function: reorder the earthquake id. if the earthquake has no data, the id is -999
def reorder_src(ev_info):
    
    ev_id = 0
    for key_ev in ev_info:
        ev      = ev_info[key_ev]

        if(ev.Nt + ev.Ncs_dt + ev.Ncr_dt  == 0):
            ev.id = -999
        else:
            ev_info[key_ev].id = ev_id
            ev_id += 1

    return ev_info


# %%
# 函数：read_src_rec_file(fname), 读取到时数据. function:  read src_rec.dat file.  
# 
def read_src_rec_file(fname):
    ev_info = {}
    st_info = {}

    tmp_ev_info = {}

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
            ev.Nt   = 0
            ev.Ncs_dt  = 0
            ev.Ncr_dt  = 0
            ev.t    = {}
            ev.cs_dt   = {}
            ev.cr_dt   = {}
            ndata   = int(tmp[11])
            name_ev = tmp[12]
            ev.name = name_ev
            cc += 1
            try:
                ev.tag["weight"] = float(tmp[13])
            except:
                pass

            if (ndata == 0):
                cc = 0
                ev_info[name_ev] = ev

        else:   # data line
            #  1      1      MYA       38.3261     38.4253   1050.0000  P   52.46   6.630 weight
            if(len(tmp) < 10):  # absolue traveltime data  绝对走时数据
                name_st = tmp[2]
                phase   = tmp[6]
                if (phase == "PG"):
                    phase = "Pg"
                if (phase == "PB"):
                    phase = "Pb"
                if (phase == "PN"):
                    phase = "Pn"    

                if (not name_st in st_info):
                    st = Station()
                    st.name = name_st
                    st.id   = float(tmp[1])
                    st.lat  = float(tmp[3])
                    st.lon  = float(tmp[4])
                    st.ele  = float(tmp[5])
                    st_info[name_st] = st
                
                time    = float(tmp[7])
                if(len(tmp) == 9):
                    weight = float(tmp[8])
                else:
                    weight = 1.0
                ev.t["%s+%s"%(name_st,phase)] = [name_st,phase,time,weight]
                ev.Nt += 1

            else:   # differential traveltime data  双差走时数据
                phase = tmp[11]
                if (phase.__contains__("cr")):  # common receiver differential traveltime 共台站双差
                    #  evid     stid1      stname1      lat1 lon1 eve1 evid2 evname2 lat2 lon2 dep2 phase,cr diftime weight 

                    name_st1 = tmp[2]
                    if (not name_st1 in st_info):   # add station to the station list 把台站加入到台站列表中
                        st = Station()
                        st.name = name_st1
                        st.id   = float(tmp[1])
                        st.lat  = float(tmp[3])
                        st.lon  = float(tmp[4])
                        st.ele  = float(tmp[5])
                        st_info[name_st1] = st
                    
                    name_ev2 = tmp[7]
                                                    # add earthquake to the temp earthquake list 把地震加入到地震列表中
                    ev2 = Event()   
                    ev2.name = name_ev2
                    ev2.id   = float(tmp[6])
                    ev2.lat  = float(tmp[8])
                    ev2.lon  = float(tmp[9])
                    ev2.dep  = float(tmp[10])
                    tmp_ev_info[name_ev2] = ev2 


                    dif_time = float(tmp[12])
                    if(len(tmp) == 14):
                        weight = float(tmp[13])
                    else:
                        weight = 1.0
                    ev.cr_dt["%s+%s+%s"%(name_st1,name_ev2,phase)] = [name_st1,name_ev2,phase,dif_time,weight]
                    ev.Ncr_dt += 1

                else:                           # common source differential traveltime 共台站双差
                    #  evid     stid1      stname1      lat1 lon1 eve1 stid2 stname2 lat2 lon2 ele2 phase,cs diftime weight 

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
                        st.id   = float(tmp[6])
                        st.lat  = float(tmp[8])
                        st.lon  = float(tmp[9])
                        st.ele  = float(tmp[10])
                        st_info[name_st2] = st
 
                    dif_time = float(tmp[12])
                    if(len(tmp) == 14):
                        weight = float(tmp[13])
                    else:
                        weight = 1.0
                    ev.cs_dt["%s+%s+%s"%(name_st1,name_st2,phase)] = [name_st1,name_st2,phase,dif_time,weight]
                    ev.Ncs_dt += 1

            if (cc == ndata):   # end of the event data
                cc = 0
                ev_info[name_ev] = ev
            else:
                cc += 1

    # 将临时地震列表中的地震加入到地震列表中
    for key_ev in tmp_ev_info:
        if (not key_ev in ev_info):
            ev_info[key_ev] = tmp_ev_info[key_ev]

    return [ev_info,st_info]

# %%
# 函数: write_src_rec_file(fname,ev_info,st_info) 输出src_rec到时数据. function: output the src_rec.dat file
def write_src_rec_file(fname,ev_info,st_info):
    ev_info = reorder_src(ev_info)
    doc_src_rec = open(fname,'w')

    min_lat =  9999
    max_lat = -9999
    min_lon =  9999
    max_lon = -9999
    min_dep =  9999
    max_dep = -9999
    
    record_ev = {}
    record_st = {}
    Nt_total  = 0
    Ncs_dt_total = 0
    Ncr_dt_total = 0

    for key_ev in ev_info:
        ev      = ev_info[key_ev]
        evid    = ev.id
        year    = ev.ortime.year
        month   = ev.ortime.month
        day     = ev.ortime.day
        hour    = ev.ortime.hour
        minute  = ev.ortime.minute
        second  = ev.ortime.second
        msec    = ev.ortime.microsecond
        lat_ev  = ev.lat
        lon_ev  = ev.lon
        dep_ev  = ev.dep
        mag     = ev.mag
        ndata   = ev.Nt + ev.Ncs_dt + ev.Ncr_dt
        name_ev = ev.name
        try:
            weight_ev = ev.tag["weight"]
        except:
            weight_ev = 1.0

        if(ndata == 0):     # 如果这个地震没有数据，则不输出这个地震. if the earthquake has no data, do not output it
            continue

        doc_src_rec.write('%7d %6d %2d %2d %2d %2d %5.2f %9.4f %9.4f %9.4f %5.2f %7d %s %7.3f\n'%(\
            evid,year,month,day,hour,minute,second+msec/1000000,lat_ev,lon_ev,dep_ev,mag,ndata,name_ev,weight_ev))
        
        min_lat =  min(min_lat, lat_ev)
        max_lat =  max(max_lat, lat_ev)
        min_lon =  min(min_lon, lon_ev)
        max_lon =  max(max_lon, lon_ev)
        min_dep =  min(min_dep, dep_ev)
        max_dep =  max(max_dep, dep_ev)

        record_ev[name_ev] = 1  # record this earthquake 记录下这个地震
        Nt_total += ev.Nt
        Ncs_dt_total += ev.Ncs_dt
        Ncr_dt_total += ev.Ncr_dt

        for key_t in ev.t:
            data    = ev.t[key_t]
            st      = st_info[data[0]]
            stid    = st.id
            name_st = st.name
            lat_st  = st.lat
            lon_st  = st.lon
            ele_st  = st.ele
            phase   = data[1]
            time    = data[2]
            try:
                weight_data = data[3]
            except:
                weight_data = 1.0
            doc_src_rec.write('%7d %7d %6s %9.4f %9.4f %9.4f %s %8.4f %7.3f \n'%(evid,stid,name_st,lat_st,lon_st,ele_st,phase,time,weight_data))

            min_lat =  min(min_lat, lat_st)
            max_lat =  max(max_lat, lat_st)
            min_lon =  min(min_lon, lon_st)
            max_lon =  max(max_lon, lon_st)
            min_dep =  min(min_dep, -ele_st/1000)
            max_dep =  max(max_dep, -ele_st/1000)

            record_st[name_st] = 1  # record this station 记录下这个台站

        for key_t in ev.cs_dt:
            data    = ev.cs_dt[key_t]
            st1     = st_info[data[0]]
            stid1   = st1.id
            name_st1= st1.name
            lat_st1 = st1.lat
            lon_st1 = st1.lon
            ele_st1 = st1.ele
            st2     = st_info[data[1]]
            stid2   = st2.id
            name_st2= st2.name
            lat_st2 = st2.lat
            lon_st2 = st2.lon
            ele_st2 = st2.ele
            phase   = data[2]
            time    = data[3]
            try:
                weight_data = data[4]
            except:
                weight_data = 1.0
            doc_src_rec.write('%7d %7d %6s %9.4f %9.4f %9.4f %7d %6s %9.4f %9.4f %9.4f %s %8.4f %7.3f \n'%(\
                evid,stid1,name_st1,lat_st1,lon_st1,ele_st1,stid2,name_st2,lat_st2,lon_st2,ele_st2,phase,time,weight_data))
            
            min_lat =  min(min_lat, lat_st1)
            max_lat =  max(max_lat, lat_st1)
            min_lon =  min(min_lon, lon_st1)
            max_lon =  max(max_lon, lon_st1)
            min_dep =  min(min_dep, -ele_st1/1000)
            max_dep =  max(max_dep, -ele_st1/1000)

            min_lat =  min(min_lat, lat_st2)
            max_lat =  max(max_lat, lat_st2)
            min_lon =  min(min_lon, lon_st2)
            max_lon =  max(max_lon, lon_st2)
            min_dep =  min(min_dep, -ele_st2/1000)
            max_dep =  max(max_dep, -ele_st2/1000)

            record_st[name_st1] = 1  # record this station 记录下这个台站
            record_st[name_st2] = 1  # record this station 记录下这个台站

        for key_t in ev.cr_dt:
            data    = ev.cr_dt[key_t]
            st      = st_info[data[0]]
            stid    = st.id
            name_st = st.name
            lat_st  = st.lat
            lon_st  = st.lon
            ele_st  = st.ele
            ev2     = ev_info[data[1]]
            evid2   = ev2.id
            name_ev2= ev2.name
            lat_ev2 = ev2.lat
            lon_ev2 = ev2.lon
            dep_ev2 = ev2.dep
            phase   = data[2]
            time    = data[3]
            try:
                weight_data = data[4]
            except:
                weight_data = 1.0
            doc_src_rec.write('%7d %7d %6s %9.4f %9.4f %9.4f %7d %6s %9.4f %9.4f %9.4f %s %8.4f %7.3f \n'%(\
                evid,stid,name_st,lat_st,lon_st,ele_st,evid2,name_ev2,lat_ev2,lon_ev2,dep_ev2,phase,time,weight_data))
            
            min_lat =  min(min_lat, lat_st)
            max_lat =  max(max_lat, lat_st)
            min_lon =  min(min_lon, lon_st)
            max_lon =  max(max_lon, lon_st)
            min_dep =  min(min_dep, -ele_st/1000)
            max_dep =  max(max_dep, -ele_st/1000)

            min_lat =  min(min_lat, lat_ev2)
            max_lat =  max(max_lat, lat_ev2)
            min_lon =  min(min_lon, lon_ev2)
            max_lon =  max(max_lon, lon_ev2)
            min_dep =  min(min_dep, dep_ev2)
            max_dep =  max(max_dep, dep_ev2)
            
            record_ev[name_ev2] = 1  # record this station 记录下这个台站
            record_st[name_st] = 1   # record this station 记录下这个台站

    doc_src_rec.close() 

    print("src_rec.dat has been outputed: %d events, %d stations, %d abs traveltime, %d cs_dif traveltime, %d cr_dif traveltime. " \
          %(len(record_ev),len(record_st),Nt_total,Ncs_dt_total,Ncr_dt_total))
    print("earthquake and station region, lat: %6.1f - %6.1f, lon: %6.1f - %6.1f, dep: %6.1f - %6.1f"%(min_lat,max_lat,min_lon,max_lon,min_dep,max_dep)  )


# %%
# 函数: write_src_list_file(fname,ev_info) 输出地震列表. function: output the event list file
def write_src_list_file(fname,ev_info):
    doc_ev_list = open(fname,'w')

    for key_ev in ev_info:
        ev      = ev_info[key_ev]
        evid    = ev.id
        lat_ev  = ev.lat
        lon_ev  = ev.lon
        dep_ev  = ev.dep
        mag     = ev.mag
        name_ev = ev.name
        if (ev.id == -999):     # 如果这个地震没有数据，则不输出这个地震. if the earthquake has no data, do not output it
            continue
        doc_ev_list.write("%7d %s %s %9.4f %9.4f %9.4f %5.2f \n"%(evid,name_ev,ev.ortime,lat_ev,lon_ev,dep_ev,mag))
    doc_ev_list.close()

# %%
# 函数: write_rec_list_file(fname,ev_info,st_info) 输出台站数据. function: output the station list file
def write_rec_list_file(fname,ev_info,st_info):
    doc_st_list = open(fname,'w')

    st_list = {}
    for key_ev in ev_info:
        ev      = ev_info[key_ev]
        
        for key_t in ev.t:
            data    = ev.t[key_t]
            st      = st_info[data[0]]
            name_st = st.name
            lat_st  = st.lat
            lon_st  = st.lon
            ele_st  = st.ele
            if(not name_st in st_list):
                doc_st_list.write("%6s %9.4f %9.4f %10.4f \n"%(name_st,lat_st,lon_st,ele_st))
                st_list[name_st] = 1

        for key_t in ev.cs_dt:
            data    = ev.cs_dt[key_t]
            st1     = st_info[data[0]]
            name_st1= st1.name
            lat_st1 = st1.lat
            lon_st1 = st1.lon
            ele_st1 = st1.ele
            st2     = st_info[data[0]]
            name_st2= st2.name
            lat_st2 = st2.lat
            lon_st2 = st2.lon
            ele_st2 = st2.ele
            if(not name_st1 in st_list):
                doc_st_list.write("%6s %9.4f %9.4f %10.4f \n"%(name_st1,lat_st1,lon_st1,ele_st1))
                st_list[name_st1] = 1
            if(not name_st2 in st_list):
                doc_st_list.write("%6s %9.4f %9.4f %10.4f \n"%(name_st2,lat_st2,lon_st2,ele_st2))
                st_list[name_st2] = 1

        for key_t in ev.cr_dt:
            data    = ev.cr_dt[key_t]
            st      = st_info[data[0]]
            name_st = st.name
            lat_st  = st.lat
            lon_st  = st.lon
            ele_st  = st.ele
            if(not name_st in st_list):
                doc_st_list.write("%6s %9.4f %9.4f %10.4f \n"%(name_st,lat_st,lon_st,ele_st))
                st_list[name_st] = 1

    doc_st_list.close()

# %% [markdown]
# 函数类：读取反演网格文件
# 
# Functions: read inversion grid file

# %%
def read_inversion_grid_file(path):
    
    inv_grid_vel = []
    inv_grid_ani = []

    switch = False
    igrid = -1
    with open('%s/inversion_grid.txt'%(path)) as f:
        tmp_inv_grid = []
        for i,line in enumerate(f):

            # read the number of inversion grid in dep, lat, lon directions
            if(i==0):
                tmp = line.split()
                ndep = int(tmp[1])
                nlines = 3*ndep+1   # 每组反演网格的行数为 3*ndep+1
            
            iline = i % nlines

            if(iline == 0):    # info: number of inversion grid
                tmp = line.split()
                if (int(tmp[0]) > igrid):
                    igrid = int(tmp[0])
                else:   # change from vel to ani
                    switch = True
                    igrid = int(tmp[0])

            else:               # info location of inversion grid
                iline_sub = (iline-1) % 3
                if(iline_sub == 0): # dep
                    tmp = line.split()
                    dep = float(tmp[0])
                if(iline_sub == 1): # list of lat 
                    lat_list = line.split()
                if(iline_sub == 2): # list of lon
                    lon_list = line.split()
                    
                    # add inversion grid
                    for lat in lat_list:
                        for lon in lon_list:
                            tmp_inv_grid.append([float(lon), float(lat), dep])
                
                if(iline == nlines-1): # the last line of inversion grid
                    if(switch):
                        inv_grid_ani.append(tmp_inv_grid)
                    else:
                        inv_grid_vel.append(tmp_inv_grid)
                    tmp_inv_grid = []
    
    return [np.array(inv_grid_vel),np.array(inv_grid_ani)]

# %% [markdown]
# 函数类：包含画图函数
# 
# Functions: for plotting

# %%
# 画图：地震和台站分布，按照深度染色 fig_ev_st_distribution_dep(ev_info,st_info) plot: plot the distribution of the earthquake and stations, color-coded by ev depth
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def fig_ev_st_distribution_dep(ev_info,st_info):
    
    [lon_ev,lat_ev,dep_ev,wt_ev] = data_lon_lat_dep_wt_ev(ev_info)
    [lon_st,lat_st,ele_st,wt_st] = data_lon_lat_ele_wt_st(ev_info,st_info)

    min_lon = min(min(lon_ev),min(lon_st))
    max_lon = max(max(lon_ev),max(lon_st))
    
    min_lat = min(min(lat_ev),min(lat_st))
    max_lat = max(max(lat_ev),max(lat_st))
    
    max_dep = max(dep_ev)

    # 插入一个不影响画图的值来使得 colorbar 的范围更好看
    lon_ev = np.insert(lon_ev, 0, 9999); lat_ev = np.insert(lat_ev, 0, 9999); dep_ev = np.insert(dep_ev, 0, 0);

    fig = plt.figure(figsize=(12,12))
    gridspace = GridSpec(12,12,figure = fig)

    xrange = max_lon - min_lon + 1.0
    yrange = max_lat - min_lat + 1.0
    
    if (xrange > yrange):
        fig_x_size = 6
        fig_y_size = round(6*yrange/xrange)
    else:
        fig_x_size = round(6*xrange/yrange)
        fig_y_size = 6


    ax1 = fig.add_subplot(gridspace[0:fig_y_size,0:fig_x_size])

    bar_ev = ax1.scatter(lon_ev,lat_ev,c=dep_ev,cmap="jet",label = "src",s = 3)
    # ax1.plot(lon_st,lat_st,'rv',label = "rec",markersize = 6)
    bar_st = ax1.scatter(lon_st,lat_st,c="red",label = "rec",s = 100,marker='v',edgecolors='white')

    ax1.legend(fontsize = 14)
    ax1.tick_params(axis='x',labelsize=18)
    ax1.tick_params(axis='y',labelsize=18)
    ax1.set_xlabel('Lon',fontsize=18)
    ax1.set_ylabel('Lat',fontsize=18)
    ax1.set_xlim((min_lon - (max_lon - min_lon)*0.1,max_lon + (max_lon - min_lon)*0.1))
    ax1.set_ylim((min_lat - (max_lat - min_lat)*0.1,max_lat + (max_lat - min_lat)*0.1))


    ax2 = fig.add_subplot(gridspace[0:fig_y_size, fig_x_size+1 : fig_x_size+3])

    ax2.scatter(dep_ev,lat_ev,c=dep_ev,cmap="jet",label = "src",s = 3)
 
    ax2.tick_params(axis='x',labelsize=18)
    ax2.tick_params(axis='y',labelsize=18)
    ax2.set_xlabel('Dep',fontsize=18)
    ax2.set_ylabel('Lat',fontsize=18)
    ax2.set_xlim((-max_dep*0.05,max_dep*1.1))
    ax2.set_ylim((min_lat - (max_lat - min_lat)*0.1,max_lat + (max_lat - min_lat)*0.1))


    ax3 = fig.add_subplot(gridspace[fig_y_size+1:fig_y_size+3,0:fig_x_size])

    ax3.scatter(lon_ev,dep_ev,c=dep_ev,cmap="jet",label = "src",s = 3)
 
    ax3.tick_params(axis='x',labelsize=18)
    ax3.tick_params(axis='y',labelsize=18)
    ax3.set_xlabel('Lon',fontsize=18)
    ax3.set_ylabel('Dep',fontsize=18)
    ax3.set_xlim((min_lon - (max_lon - min_lon)*0.1,max_lon + (max_lon - min_lon)*0.1))
    ax3.set_ylim((-max_dep*0.05,max_dep*1.1))
    ax3.invert_yaxis()

    # 将 colorbar 放置在新的轴上
    ax4 = fig.add_subplot(gridspace[fig_y_size+2:fig_y_size+3,fig_x_size+1:fig_x_size+3])
    cbar1 = plt.colorbar(bar_ev, ax=ax4,orientation='horizontal')
    cbar1.set_label('Depth of earthquakes',fontsize=16)
    cbar1.ax.tick_params(axis='x', labelsize=16)    # colorbar的字体大小
    #隐藏坐标轴的边框
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.spines['bottom'].set_visible(False)
    ax4.spines['left'].set_visible(False)

    # 隐藏坐标轴刻度值
    ax4.set_xticks([])
    ax4.set_yticks([])


# %%
# 画图：地震和台站分布，按照权重染色 fig_ev_st_distribution_wt(ev_info,st_info) plot: plot the distribution of the earthquake and stations, color-coded by weight
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math

def fig_ev_st_distribution_wt(ev_info,st_info):
    
    
    [lon_ev,lat_ev,dep_ev,wt_ev] = data_lon_lat_dep_wt_ev(ev_info)
    [lon_st,lat_st,ele_st,wt_st] = data_lon_lat_ele_wt_st(ev_info,st_info)

    # 插入一个不影响画图的值来使得 colorbar 的范围更好看
    lon_ev = np.insert(lon_ev, 0, lon_ev[0]); lat_ev = np.insert(lat_ev, 0, lat_ev[0]); dep_ev = np.insert(dep_ev, 0, dep_ev[0]); wt_ev = np.insert(wt_ev, 0, 0.0)
    lon_ev = np.insert(lon_ev, 0, lon_ev[0]); lat_ev = np.insert(lat_ev, 0, lat_ev[0]); dep_ev = np.insert(dep_ev, 0, dep_ev[0]); wt_ev = np.insert(wt_ev, 0, 1.0)
    lon_st = np.insert(lon_st, 0, lon_st[0]); lat_st = np.insert(lat_st, 0, lat_st[0]); ele_st = np.insert(ele_st, 0, ele_st[0]); wt_st = np.insert(wt_st, 0, 0.0)
    lon_st = np.insert(lon_st, 0, lon_st[0]); lat_st = np.insert(lat_st, 0, lat_st[0]); ele_st = np.insert(ele_st, 0, ele_st[0]); wt_st = np.insert(wt_st, 0, 1.0)

    min_lon = min(min(lon_ev),min(lon_st))
    max_lon = max(max(lon_ev),max(lon_st))
    
    min_lat = min(min(lat_ev),min(lat_st))
    max_lat = max(max(lat_ev),max(lat_st))
    
    max_dep = max(dep_ev)

    fig = plt.figure(figsize=(12,12))
    gridspace = GridSpec(12,12,figure = fig)

    xrange = max_lon - min_lon + 1.0
    yrange = max_lat - min_lat + 1.0

    if (xrange > yrange):
        fig_x_size = 6
        fig_y_size = round(6*yrange/xrange)
    else:
        fig_x_size = round(6*xrange/yrange)
        fig_y_size = 6


    ax1 = fig.add_subplot(gridspace[0:fig_y_size,0:fig_x_size])

    bar_ev = ax1.scatter(lon_ev,lat_ev,c=wt_ev,cmap="jet",label = "src",s = 3)
    bar_st = ax1.scatter(lon_st,lat_st,c=wt_st,cmap="jet",label = "rec",s = 100,marker='^',edgecolors='white')
    
    ax1.legend(fontsize = 14)
    ax1.tick_params(axis='x',labelsize=18)
    ax1.tick_params(axis='y',labelsize=18)
    ax1.set_xlabel('Lon',fontsize=18)
    ax1.set_ylabel('Lat',fontsize=18)
    ax1.set_xlim((min_lon - (max_lon - min_lon)*0.1,max_lon + (max_lon - min_lon)*0.1))
    ax1.set_ylim((min_lat - (max_lat - min_lat)*0.1,max_lat + (max_lat - min_lat)*0.1))


    ax2 = fig.add_subplot(gridspace[0:fig_y_size, fig_x_size+1 : fig_x_size+3])

    ax2.scatter(dep_ev,lat_ev,c=wt_ev,cmap="jet",label = "src",s = 3)
 
    ax2.tick_params(axis='x',labelsize=18)
    ax2.tick_params(axis='y',labelsize=18)
    ax2.set_xlabel('Dep',fontsize=18)
    ax2.set_ylabel('Lat',fontsize=18)
    ax2.set_xlim((-max_dep*0.05,max_dep*1.1))
    ax2.set_ylim((min_lat - (max_lat - min_lat)*0.1,max_lat + (max_lat - min_lat)*0.1))


    ax3 = fig.add_subplot(gridspace[fig_y_size+1:fig_y_size+3,0:fig_x_size])

    ax3.scatter(lon_ev,dep_ev,c=wt_ev,cmap="jet",label = "src",s = 3)
 
    ax3.tick_params(axis='x',labelsize=18)
    ax3.tick_params(axis='y',labelsize=18)
    ax3.set_xlabel('Lon',fontsize=18)
    ax3.set_ylabel('Dep',fontsize=18)
    ax3.set_xlim((min_lon - (max_lon - min_lon)*0.1,max_lon + (max_lon - min_lon)*0.1))
    ax3.set_ylim((-max_dep*0.05,max_dep*1.1))
    ax3.invert_yaxis()

    # 将 colorbar 放置在新的轴上
    ax4 = fig.add_subplot(gridspace[fig_y_size+2:fig_y_size+3,fig_x_size+1:fig_x_size+3])
    cbar1 = plt.colorbar(bar_st, ax=ax4,orientation='horizontal')
    cbar1.set_label('Weight of stations',fontsize=16)
    cbar1.ax.tick_params(axis='x', labelsize=16)    # colorbar的字体大小
    #隐藏坐标轴的边框
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.spines['bottom'].set_visible(False)
    ax4.spines['left'].set_visible(False)

    # 隐藏坐标轴刻度值
    ax4.set_xticks([])
    ax4.set_yticks([])

    # 将 colorbar 放置在新的轴上
    ax5 = fig.add_subplot(gridspace[fig_y_size+1:fig_y_size+2,fig_x_size+1:fig_x_size+3])
    cbar1 = plt.colorbar(bar_ev, ax=ax5,orientation='horizontal')
    cbar1.set_label('Weight of earthquakes',fontsize=16)    
    cbar1.ax.tick_params(axis='x', labelsize=16)    # colorbar的字体大小
    #隐藏坐标轴的边框
    ax5.spines['top'].set_visible(False)
    ax5.spines['right'].set_visible(False)
    ax5.spines['bottom'].set_visible(False)
    ax5.spines['left'].set_visible(False)

    # 隐藏坐标轴刻度值
    ax5.set_xticks([])
    ax5.set_yticks([])


# %%
# 画图和函数：去除走时数据的异常值，画出走时距离散点图 plot and function: plot the distance-time scatter plot, remove the outliers.
# 将数据限制在 直线 time = dis * slope + intercept 的两侧， up and down 的范围里面
# remove outliers, only retain data satisfying:     slope * dis + intercept + down < time < slope * dis + intercept + up

def fig_data_plot_remove_outliers(ev_info,st_info,slope,intercept,up,down,dis_min,dis_max):   

    fig = plt.figure(figsize=(10,10))
    gridspace = GridSpec(6,6,figure = fig)
    ax2 = fig.add_subplot(gridspace[0:6, 0:6])

    # plot original data
    [dis_obs,time_obs] = data_dis_time(ev_info,st_info)
    ax2.plot(dis_obs,time_obs,'r.',markersize=1.5,label = "discarded")

    # 去除异常值，只保留满足条件的数据
    # remove outliers, only retain data satisfying:     slope * dis + intercept + down < time < slope * dis + intercept + up
    ev_info = limit_data_residual(ev_info,st_info,slope,intercept,up,down)

    [dis_obs,time_obs] = data_dis_time(ev_info,st_info)
    ax2.plot(dis_obs,time_obs,'b.',markersize=1.5,label = "retained")

    ax2.plot([dis_min,dis_max],[slope*dis_min+intercept+up,slope*dis_max+intercept+up],'b-',linewidth=2)
    ax2.plot([dis_min,dis_max],[slope*dis_min+intercept+down,slope*dis_max+intercept+down],'b-',linewidth=2)
    ax2.plot([dis_min,dis_max],[slope*dis_min+intercept,slope*dis_max+intercept],'k-',linewidth=2)

    ax2.legend(fontsize = 14)
    ax2.tick_params(axis='x',labelsize=18)
    ax2.tick_params(axis='y',labelsize=18)
    ax2.set_xlabel('Distance (km)',fontsize=18)
    ax2.set_ylabel('Traveltime',fontsize=18)
    ax2.set_xlim((dis_min,dis_max))
    ax2.set_ylim((intercept+down-5,slope*dis_max+intercept+up+5))

    return ev_info

# %%
# 画图: 画出不同震相的走时距离散点图 plot: distance-time scatter plot of given phases

def fig_data_plot_phase(ev_info,st_info,phase_list,color_list,dis_min,dis_max):   

    [dis_obs_phase,time_obs_phase] = data_dis_time_phase(ev_info,st_info,phase_list)

    regression = {}
    # 计算最小二乘 y = ax+b
    for key_phase in phase_list:
        X = dis_obs_phase[key_phase]
        Y = time_obs_phase[key_phase]

        if(len(X)>20):
            regression[key_phase] = linear_regression(X,Y)
        else:
            print("No enough data: %d, for %s"%(len(X),key_phase))
            regression[key_phase] = [0,0,0]


    # 画图
    fig = plt.figure(figsize=(10,10))
    gridspace = GridSpec(6,6,figure = fig)
    ax2 = fig.add_subplot(gridspace[0:6, 0:6])
    y1 = 99999; y2 = -99999

    # scatter plot
    for iphase in range(len(phase_list)):
        phase = phase_list[iphase]
        color = color_list[iphase]
        ax2.plot(dis_obs_phase[phase],time_obs_phase[phase],'%s.'%(color),markersize=1)

    # linear regression plot
    for iphase in range(len(phase_list)):
        phase = phase_list[iphase]
        color = color_list[iphase]
        (slope,intercept,SEE)= regression[phase]
        ax2.plot([dis_min,dis_max],[dis_min*slope+intercept,dis_max*slope+intercept],'%s-'%(color),linewidth=2,label = "%s: a,b,SEE=%5.2f,%5.2f,%5.2f"%(phase,slope,intercept,SEE))
        y1 = min(y1,intercept-5)
        y2 = max(y2,dis_max*slope+intercept+5)

    ax2.legend(fontsize = 14)
    ax2.tick_params(axis='x',labelsize=18)
    ax2.tick_params(axis='y',labelsize=18)
    ax2.set_xlabel('Distance (km)',fontsize=18)
    ax2.set_ylabel('Traveltime (s)',fontsize=18)
    ax2.set_xlim((dis_min,dis_max))

    
    for iphase in range(len(phase_list)):
        try:
            y1 = min(y1,min(time_obs_phase[phase]))
            y2 = max(y2,max(time_obs_phase[phase]))
        except:
            pass
    ax2.set_ylim((y1,y2))

    title = ""
    for phase in dis_obs_phase:
        title = title + "%s(%d) "%(phase,len(dis_obs_phase[phase]))

    ax2.set_title(title)

    print("a is slope, b is intercept, SEE is standard error of estimate")

# %% [markdown]
# 函数类：结果分析处理
# 
# Functions: results analysis and evaluation

# %%
# 画图：画出初始模型和最终模型的残差直方图 plot: plot the residual histogram of the initial and final model
def fig_residual_histogram(fn_syn_init,fn_syn_final,fn_obs,range_l,range_r,Nbar,tag1 = "initial",tag2 = "final"):  
    
    # 读取初始合成到时数据 read synthetic traveltime data in the initial model
    [ev_info_syn_init, st_info_syn_init] = read_src_rec_file(fn_syn_init)
    time_syn_init = data_dis_time(ev_info_syn_init,st_info_syn_init)[1]

    # 读取最终模型合成到时数据 read synthetic traveltime data in the final model
    [ev_info_syn_final, st_info_syn_final] = read_src_rec_file(fn_syn_final)
    time_syn_final = data_dis_time(ev_info_syn_final,st_info_syn_final)[1]

    # 读取观测到时数据 read observed traveltime data
    [ev_info_obs, st_info_obs] = read_src_rec_file(fn_obs)
    time_obs = data_dis_time(ev_info_obs,st_info_obs)[1] 

    fig = plt.figure(figsize=(6,6))
    gridspace = GridSpec(6,6,figure = fig)

    ax2 = fig.add_subplot(gridspace[0:6, 0:6])

    bins=np.linspace(range_l,range_r,Nbar)
    error_init = time_syn_init - time_obs
    error_final = time_syn_final - time_obs

    hist_init,  _, _ = ax2.hist(error_init,bins=bins,histtype='step', edgecolor = "red", linewidth = 2,
            label = "%s: std = %5.3f s, mean = %5.3f s"%(tag1,np.std(error_init),np.mean(error_init)))
    
    hist_final, _, _ = ax2.hist(error_final,bins=bins,alpha = 0.5, color = "blue",
            label = "%s: std = %5.3f s, mean = %5.3f s"%(tag2,np.std(error_final),np.mean(error_final)))

    print("residual for ",tag1," model is: ","mean: ",np.mean(error_init),"sd: ",np.std(error_init))
    print("residual for ",tag2," model is: ","mean: ",np.mean(error_final),"sd: ",np.std(error_final))
    ax2.legend(fontsize=14)

    ax2.set_xlim(range_l - abs(range_l)*0.1,range_r + abs(range_r)*0.1)
    ax2.set_ylim(0,1.3*max(max(hist_init),max(hist_final)))

    ax2.tick_params(axis='x',labelsize=18)
    ax2.tick_params(axis='y',labelsize=18)
    ax2.set_ylabel('Number of data',fontsize=18)
    ax2.set_xlabel('Traveltime residuals (s)',fontsize=18)
    ax2.set_title("$t_{syn} - t_{obs}$",fontsize=18)
    ax2.grid()


# %%
# 画图：画出两组文件的到时数据差别分布直方图 plot: plot the residual histogram of the initial and final model
def fig_data_difference_histogram(fn_A,fn_B,range_l,range_r,Nbar):  
      
    # 读取数据A read data A
    [ev_info_A, st_info_A] = read_src_rec_file(fn_A)

    # 读取数据B read data B
    [ev_info_B, st_info_B] = read_src_rec_file(fn_B)

    # 绝对到时残差 absolute traveltime residual
    error_t = []
    for key_ev in ev_info_A:
        for key_t in ev_info_A[key_ev].t:
            data_A = ev_info_A[key_ev].t[key_t]
            if (key_ev in ev_info_B and key_t in ev_info_B[key_ev].t):
                data_B = ev_info_B[key_ev].t[key_t]
                error_t.append(data_A[2] - data_B[2])

    # 共源双差到时残差 common-source differential traveltime residual
    error_cs_dt = []
    for key_ev in ev_info_A:
        for key_dt in ev_info_A[key_ev].cs_dt:
            data_A = ev_info_A[key_ev].cs_dt[key_dt]
            if (key_ev in ev_info_B and key_dt in ev_info_B[key_ev].cs_dt):
                data_B = ev_info_B[key_ev].cs_dt[key_dt]
                error_cs_dt.append(data_A[3] - data_B[3])
            else:
                print(key_ev,key_dt)

    # 共台双差到时残差 common-receiver differential traveltime residual
    error_cr_dt = []
    for key_ev in ev_info_A:
        for key_dt in ev_info_A[key_ev].cr_dt:
            data_A = ev_info_A[key_ev].cr_dt[key_dt]
            if (key_ev in ev_info_B and key_dt in ev_info_B[key_ev].cr_dt):
                data_B = ev_info_B[key_ev].cr_dt[key_dt]
                error_cr_dt.append(data_A[3] - data_B[3])

    # 画图
    fig = plt.figure(figsize=(14,6))
    gridspace = GridSpec(6,14,figure = fig)
    

    ax2 = fig.add_subplot(gridspace[0:6, 0:6])
    bins=np.linspace(range_l,range_r,Nbar)
    # hist_t,  _, _ = ax2.hist(error_t,bins=bins,histtype='step', edgecolor = "red", linewidth = 2,
    #         label = "noise: std = %5.3f s, mean = %5.3f s"%(np.std(error_t),np.mean(error_t)))
    hist_t, _, _ = ax2.hist(error_t,bins=bins,alpha = 0.5, color = "blue",
            label = "noise: std = %5.3f s, mean = %5.3f s"%(np.std(error_t),np.mean(error_t)))

    ax2.legend(fontsize=14)
    ax2.set_xlim(range_l - abs(range_l)*0.1,range_r + abs(range_r)*0.1)
    try:
        ax2.set_ylim(0,1.3*max(hist_t))
    except: 
        ax2.set_ylim(0,1.0)
    ax2.tick_params(axis='x',labelsize=18)
    ax2.tick_params(axis='y',labelsize=18)
    ax2.set_ylabel('Number of data',fontsize=18)
    ax2.set_xlabel('Noise (s)',fontsize=18)
    ax2.set_title("Noise of traveltime",fontsize=18)


    ax3 = fig.add_subplot(gridspace[0:6,8:14])
    bins=np.linspace(range_l,range_r,Nbar)
    # hist_t,  _, _ = ax3.hist(error_t,bins=bins,histtype='step', edgecolor = "red", linewidth = 2,
    #         label = "noise: std = %5.3f s, mean = %5.3f s"%(np.std(error_t),np.mean(error_t)))
    hist_cs_dt, _, _ = ax3.hist(error_cs_dt,bins=bins,alpha = 0.5, color = "blue",
            label = "noise: std = %5.3f s, mean = %5.3f s"%(np.std(error_cs_dt),np.mean(error_cs_dt)))

    ax3.legend(fontsize=14)
    ax3.set_xlim(range_l - abs(range_l)*0.1,range_r + abs(range_r)*0.1)
    try:
        ax3.set_ylim(0,1.3*max(hist_cs_dt))
    except: 
        ax3.set_ylim(0,1.0)
    ax3.tick_params(axis='x',labelsize=18)
    ax3.tick_params(axis='y',labelsize=18)
    ax3.set_ylabel('Number of data',fontsize=18)
    ax3.set_xlabel('Noise (s)',fontsize=18)
    ax3.set_title("Noise of differential traveltime",fontsize=18)


