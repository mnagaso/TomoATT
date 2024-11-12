# %%
import pygmt
pygmt.config(FONT="16p", IO_SEGMENT_MARKER="<<<")

import os

# %%
from pytomoatt.model import ATTModel
from pytomoatt.data import ATTData
import numpy as np

# %%
# read model files

Ngrid = [121,121,121]
data_file = '2_models/model_init_N%d_%d_%d.h5'%(Ngrid[0],Ngrid[1],Ngrid[2])
par_file = '3_input_params/input_params_real.yaml'  
model = ATTModel.read(data_file, par_file)
init_model = model.to_xarray()

data_file = 'OUTPUT_FILES/OUTPUT_FILES_real_noshared/final_model.h5'
model = ATTModel.read(data_file, par_file)
inv_model = model.to_xarray()

# %%
# # read earthquakes and stations

# from pytomoatt.src_rec import SrcRec

# # read src_rec_file
# sr = SrcRec.read("1_src_rec_files/src_rec_file.dat")

# # get the coordinates of the stations and earthquakes
# stations = sr.receivers[['stlo','stla','stel']].values.T
# earthquakes = sr.sources[['evlo','evla','evdp']].values.T

# print(stations.shape)
# print(earthquakes.shape)

# %%
# read earthquakes and stations

import sys
sys.path.append('../utils')
import functions_for_data as ffd

ev, st = ffd.read_src_rec_file('1_src_rec_files/src_rec_file.dat')
# earthquake location
ev_lon, ev_lat, ev_dep , _ = ffd.data_lon_lat_dep_wt_ev(ev)
# station location
st_lon, st_lat, _, _ = ffd.data_lon_lat_ele_wt_st(ev,st)

# %%
# load topography
region = [97,106,12,21]
grid_topo = pygmt.datasets.load_earth_relief(resolution="01m", region=region)
grid_gra = pygmt.grdgradient(grid = grid_topo, azimuth = 0)

# %%
def line_read(file):
    doc=open(file,'r')
    file = doc.readlines()
    doc.close()
    lat = []; lon = []; 
    for info in file:
        tmp = info.split()
        lon.append(float(tmp[0]))
        lat.append(float(tmp[1]))
    return((lat,lon))

# %%
# plot imgaing results

fig = pygmt.Figure()
try:
    os.mkdir("img")
except:
    pass


# ------------------ Sub fig 1. topography ------------------
region = [97,106,12,21]
projection = "B101.5/16.5/12/21/10c"
frame = ["xa2+lLongitude", "ya2+lLatitude", "nSWe"]
spacing = [0.1, 0.1]

# topography
pygmt.makecpt(cmap="globe", series=[-4000,4000], background = True)
fig.grdimage(grid=grid_topo, shading = grid_gra, projection=projection, frame=frame,region=region)

# station
fig.plot(x = st_lon, y = st_lat, style = "t0.4c", fill = "blue", pen = "1p,white", label = "Station")


# tectonic setting
(lat,lon) = line_read('tectonics/Khorat_new.txt')   # Khorat Plateau
fig.plot(x = lon, y = lat, pen = "1p,black")    
fig.text(text = "Khorat", x = 103.5, y = 17.5, font = "15p,Helvetica-Bold,black", angle = 0)
fig.text(text = "Plateau", x = 103.5, y = 16.5, font = "15p,Helvetica-Bold,black", angle = 0)
(lat,lon) = line_read('tectonics/WangChaoFault.txt')    # Wang-Chao Fault
fig.plot(x = lon, y = lat, pen = "1p,black,-")
fig.text(text = "WCF", x = 100, y = 16, font = "20p,Helvetica-Bold,white=1p", angle = 315)
(lat,lon) = line_read('tectonics/3PagodasFault.txt')    # Three Pagodas Fault
fig.plot(x = lon, y = lat, pen = "1p,black,-")
fig.text(text = "TPF", x = 98.5, y = 14.5, font = "20p,Helvetica-Bold,white=1p", angle = 315)
(lat,lon) = line_read('tectonics/DBPF.txt') # Dien Bien Phu Fault
fig.plot(x = lon, y = lat, pen = "1p,black,-")  
fig.text(text = "DBPF", x = 102, y = 19.5, font = "20p,Helvetica-Bold,white=1p", angle = 55)
(lat,lon) = line_read('tectonics/SongMaSuture.txt')   # Song Ma Suture
fig.plot(x = lon, y = lat, pen = "1p,black,-")
fig.text(text = "Shan-Thai", x = 99, y = 20, font = "18p,Helvetica-Bold,black=0.5p,white", angle = 0)
fig.text(text = "Block", x = 99, y = 19.3, font = "18p,Helvetica-Bold,black=0.5p,white", angle = 0)
fig.text(text = "Indochina Block", x = 103.5, y = 15, font = "18p,Helvetica-Bold,black=0.5p,white", angle = 315)


fig.shift_origin(xshift= 1, yshift=-1.5)
fig.colorbar(frame = ["a%f"%(4000),"y+lElevation (m)"], position="+w4c/0.3c+h")
fig.shift_origin(xshift=-1, yshift=+1.5)

fig.shift_origin(xshift=0, yshift=-13)


# ------------------ Sub fig 2. earthquakes ------------------

region      = "g"
projection  = "E101/16/90/10c"  # centerlon/centerlat/distnace_range/fig_size
frame       = ["ya180"]
spacing     = [0.1, 5]

fig.basemap(region=region, projection=projection, frame=frame)  
fig.coast(region=region, projection=projection, frame=frame, water="white", land="gray", A=10000)

fig.plot(x=101.5, y=16.5, pen="1p,black,-", style="E-60d")
fig.plot(x=101.5, y=16.5, pen="1p,black,-", style="E-120d")
fig.plot(x=101.5, y=16.5, pen="1p,black,-", style="E-180d")
fig.text(x = [101.5, 101.5, 101.5], y = [-8.0, -38.0, -68], text = ['30', '60', '90'], font="13p")

fig.plot(x=[97,97,106,106,97], y=[11,21,21,11,11], pen="1p,red")
pygmt.makecpt(cmap="jet", series=[0, 100], background=True)
fig.plot(x = ev_lon, y = ev_lat, size = [0.4]*len(ev_lon), fill = ev_dep, style = "a", cmap = True, pen = "0.5p,white")


fig.shift_origin(xshift= 1, yshift=-1)
fig.colorbar(frame = ["a%f"%(50),"y+lFocal depth (km)"], position="+w4c/0.3c+h")
fig.shift_origin(xshift=-1, yshift=+1)

fig.shift_origin(xshift= 13, yshift= 13)


# ------------------ Sub fig 3. depth and vertical profile ------------------

region = [97,106,12,21]
projection = "B101.5/16.5/12/21/10c"
frame = ["xa2+lLongitude", "ya2+lLatitude", "nSWe"]
spacing = [0.1, 0.1]

depth_list = [100,200,300,400]
start_list = [ [97, 19],  [97, 17.2],  [101.8, 12] ]
end_list =   [ [106, 15], [106, 13.2], [101.8, 21] ]

# depth profiles
for idepth, depth in enumerate(depth_list):

    # read models
    vel_init = init_model.interp_dep(depth, field='vel')
    vel_inv = inv_model.interp_dep(depth, field='vel')

    if idepth == 3:
        fig.shift_origin(xshift=-39, yshift=-13)

    pygmt.makecpt(cmap="../utils/svel13_chen.cpt", series=[-2, 2], background=True)
    x = vel_inv[:,0]; y = vel_inv[:,1]; value = (vel_inv[:,2] - vel_init[:,2])/vel_init[:,2]*100
    grid = pygmt.surface(x=x, y=y, z=value, spacing=spacing,region=region)
    fig.grdimage(frame=frame,grid = grid,projection=projection, region=region) # nan_transparent may work

    (lat,lon) = line_read('tectonics/Khorat_new.txt')   # Khorat Plateau
    fig.plot(x = lon, y = lat, pen = "1p,black")    
    
    # vertical profile location
    fig.plot(x = [start_list[0][0],end_list[0][0]], y = [start_list[0][1],end_list[0][1],], pen = "2p,green,-")
    fig.text(text = "A",        x = start_list[0][0] + 0.5, y = start_list[0][1], font = "18p,Helvetica-Bold,black", fill = "lightblue", angle = 0)
    fig.text(text = "A@+'@+",   x = end_list[0][0] - 0.5,   y = end_list[0][1] + 0.5, font = "18p,Helvetica-Bold,black", fill = "lightblue", angle = 0)

    fig.plot(x = [start_list[1][0],end_list[1][0]], y = [start_list[1][1],end_list[1][1],], pen = "2p,green,-")
    fig.text(text = "B",        x = start_list[1][0] + 0.5, y = start_list[1][1], font = "18p,Helvetica-Bold,black", fill = "lightblue", angle = 0)
    fig.text(text = "B@+'@+",   x = end_list[1][0] - 0.5,   y = end_list[1][1] + 0.5, font = "18p,Helvetica-Bold,black", fill = "lightblue", angle = 0)

    fig.plot(x = [start_list[2][0],end_list[2][0]], y = [start_list[2][1],end_list[2][1],], pen = "2p,green,-")
    fig.text(text = "C",        x = start_list[2][0] - 0.5, y = start_list[2][1] + 0.5, font = "18p,Helvetica-Bold,black", fill = "lightblue", angle = 0)
    fig.text(text = "C@+'@+",   x = end_list[2][0] - 0.5,   y = end_list[2][1]   - 0.5, font = "18p,Helvetica-Bold,black", fill = "lightblue", angle = 0)

    # depth label
    fig.text(text="%d km"%(depth), x = 98 , y = 12.5, font = "14p,Helvetica-Bold,black", fill = "white") 


    fig.shift_origin(xshift=13)

fig.shift_origin(yshift=6)

# vertical profiles
for iprof in range(len(start_list)):

    # generate topography data
    Npoints = 100
    point_x = np.linspace(start_list[iprof][0],end_list[iprof][0],Npoints)
    point_y = np.linspace(start_list[iprof][1],end_list[iprof][1],Npoints)
    points = np.hstack((point_x.reshape(-1,1),point_y.reshape(-1,1)))
    topo = np.array(pygmt.grdtrack(points=points, grid=grid_topo)[2])
    topo_dis = [0]
    for ip in range(1, Npoints):
        dis = ffd.cal_dis(point_y[0], point_x[0], point_y[ip], point_x[ip])
        topo_dis.append(dis)
    topo_dis = np.array(topo_dis)
    
    # read models
    vel_init_sec = init_model.interp_sec(start_list[iprof], end_list[iprof], field='vel', val=10)
    vel_inv_sec = inv_model.interp_sec(start_list[iprof], end_list[iprof], field='vel', val=10)
    

    # plot topography
    max_dis = np.max(vel_init_sec[:,2])

    region = [0,max_dis,0,2000]
    projection = "X%f/1c"%(max_dis/400*4)
    frame = ["ya2000+lElevation (m)", "sW"]
    
    fig.shift_origin(yshift=4)
    fig.basemap(region=region, projection=projection, frame=frame)
    fig.plot(x = topo_dis, y = topo, pen = "1p,black", frame = frame, projection = projection, region = region)
    fig.shift_origin(yshift=-4)

    # plot model
    region = [0,max_dis,0,400]
    projection = "X%f/-4c"%(max_dis/400*4)
    frame = ["xa300+lDistance (km)", "ya100+lDepth (km)", "nSWe"]
    spacing = [10, 5]
    
    x_sec = vel_inv_sec[:,2]; y_sec = vel_inv_sec[:,3]; value_sec = (vel_inv_sec[:,4] - vel_init_sec[:,4])/vel_init_sec[:,4]*100
    grid = pygmt.surface(x=x_sec, y=y_sec, z=value_sec, spacing=spacing,region=region)
    fig.grdimage(frame=frame,grid = grid,projection=projection, region=region) # nan_transparent may work

    label_list = ['A', 'B', 'C']
    fig.text(text = "%s"%(label_list[iprof]),        x = 50, y = 50 , font = "18p,Helvetica-Bold,black", fill = "lightblue", angle = 0)
    fig.text(text = "%s@+'@+"%(label_list[iprof]),   x = np.max(x) - 50,   y = 50, font = "18p,Helvetica-Bold,black", fill = "lightblue", angle = 0)


    if (iprof == 0):
        fig.shift_origin(xshift=0, yshift=-6.5)
    elif (iprof == 1):
        fig.shift_origin(xshift=13, yshift=6.5)

fig.shift_origin(xshift= 2, yshift=-2.5)
fig.colorbar(frame = ["a%f"%(2),"y+ldlnVp (%)"], position="+w4c/0.3c+h") # +e,默认是双箭头，f表示forward，b表示background ，w表示长宽，h表示水平
fig.shift_origin(xshift=-2, yshift=+2.5)

fig.show()

fig.savefig("../img/imaging_result.png")

# %%
print(np.max(vel_init_sec[:,2]))


