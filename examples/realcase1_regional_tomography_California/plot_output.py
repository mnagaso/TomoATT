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

Ngrid = [51,89,33]
data_file = '2_models/model_init_N%d_%d_%d.h5'%(Ngrid[0],Ngrid[1],Ngrid[2])
par_file = '3_input_params/input_params_real.yaml'  
model = ATTModel.read(data_file, par_file)
initial_model = model.to_xarray()

data_file = 'OUTPUT_FILES/OUTPUT_FILES_real/final_model.h5'
model = ATTModel.read(data_file, par_file)
inv_model = model.to_xarray()

# %%
# read earthquakes and stations

from pytomoatt.src_rec import SrcRec

# read src_rec_file
sr = SrcRec.read("1_src_rec_files/src_rec_file.dat")

# rotate back to original coordinates
central_lat     =   35.6
central_lon     = -120.45
rotation_angle  = -30
sr.rotate(central_lat, central_lon, rotation_angle, reverse=True)

# get the coordinates of the stations and earthquakes
stations = sr.receivers[['stlo','stla','stel']].values.T
earthquakes = sr.sources[['evlo','evla','evdp']].values.T

print(stations.shape)
print(earthquakes.shape)

# %%
# study region

import sys
sys.path.append('../utils')
import functions_for_data as ffd

lat1 = -1.8;    lat2 = 2.2; 
lon1 = -0.7;    lon2 = 0.7; 

lat_lon_rotate = np.array([[lon1,lat1],[lon1,lat2],[lon2,lat2],[lon2,lat1],[lon1,lat1]])
lat_lon = ffd.rtp_rotation_reverse(lat_lon_rotate[:,1],lat_lon_rotate[:,0],central_lat,central_lon,rotation_angle)
studt_lat = lat_lon[0]
studt_lon = lat_lon[1]

# %%
# load topography
region = [-122.8,-118.5,33.5,38]
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
region = [-122.8,-118.5,33.5,38]
frame  = ["xa1","ya1","nSWe"]
projection = "M10c"

# topography
pygmt.makecpt(cmap="globe", series=[-4000,4000], background = True)
fig.grdimage(grid=grid_topo, shading = grid_gra, projection=projection, frame=frame,region=region)
# study region
fig.plot(x = studt_lon, y = studt_lat, pen = "1.5p,red")
# earthquakes
fig.plot(x = earthquakes[0,:], y = earthquakes[1,:], style = "c0.02c", fill = "red",label = "Earthquake")
# stations
fig.plot(x = stations[0,:], y = stations[1,:], style = "t0.2c", fill = "blue", pen = "white", label = "Station")


fig.basemap(region=[0,1,0,1], frame=["wesn+gwhite"], projection="X4c/2c")
fig.plot(x=0.1, y=0.3, style='c0.2c', fill='red')
fig.text(text="Earthquake", x=0.2, y=0.3, font="16p,Helvetica", justify="LM")
fig.plot(x=0.1, y=0.7, style='t0.4c', fill='blue', pen='black')
fig.text(text="Station", x=0.2, y=0.7, font="16p,Helvetica", justify="LM")


# ------------------ Sub fig 2. colorbar ------------------
fig.shift_origin(xshift= 2, yshift= -2)

pygmt.makecpt(cmap="globe", series=[-4000,4000], background = True)
fig.colorbar(frame = ["a%f"%(4000),"y+lElevation (m)"], position="+e+w4c/0.3c+h") 
fig.shift_origin(yshift=-2)

pygmt.makecpt(cmap="../utils/svel13_chen.cpt", series=[-8, 8], background=True, reverse=False)
fig.colorbar(frame = ["a%f"%(4),"y+ldlnVp (%)"], position="+e+w4c/0.3c+h")
fig.shift_origin(yshift=-2)

pygmt.makecpt(cmap="cool", series=[0, 0.08], background=True, reverse=False)
fig.colorbar(frame = ["a%f"%(0.04),"y+lAnisotropy"], position="+ef+w4c/0.3c+h") 



# ------------------ Sub fig 3. model ------------------
fig.shift_origin(xshift = 10, yshift=8)

region_oblique = [-0.7,0.7,-2.2,1.8]    
projection = "OA%s/%s/%s/4c"%(central_lon,central_lat,rotation_angle-90.0)
perspective = "30/90"
spacing = "1m"

depth_list = [4,8,16]

for idepth, depth in enumerate(depth_list):

    # initial model
    vel_init = initial_model.interp_dep(depth, field='vel')

    # output model
    vel_inv = inv_model.interp_dep(depth, field='vel')  # velocity 
    epsilon_inv = inv_model.interp_dep(depth, field='epsilon')  # magnitude of anisotropy

    # fast velocity directions
    samp_interval = 3
    ani_thd = 0.015
    length = 20
    width = 0.1

    ani_inv_phi = inv_model.interp_dep(depth, field='phi', samp_interval=samp_interval)
    ani_inv_epsilon = inv_model.interp_dep(depth, field='epsilon', samp_interval=samp_interval)
    ani_inv = np.hstack([ani_inv_phi, ani_inv_epsilon[:,2].reshape(-1, 1)*length, np.ones((ani_inv_epsilon.shape[0],1))*width])   # lon, lat, angle, length, width
    idx = np.where(ani_inv_epsilon[:,2] > ani_thd)
    ani = ani_inv[idx[0],:]

    # --------- plot velocity ------------
    if idepth == 0:
        frame  = ["xa100","ya1","nSwE"]
    elif idepth == len(depth_list)-1:
        frame  = ["xa100","ya1","NsWe"]
    else:
        frame  = ["xa100","ya1","nswe"]

    fig.basemap(region=region_oblique, frame=frame, projection=projection, perspective=perspective)
    pygmt.makecpt(cmap="../utils/svel13_chen.cpt", series=[-8, 8], background=True, reverse=False)

    x = vel_init[:,0]; y = vel_init[:,1]; value = (vel_inv[:,2] - vel_init[:,2])/vel_init[:,2] * 100
    y,x = ffd.rtp_rotation_reverse(y,x,central_lat,central_lon,rotation_angle)
    grid = pygmt.surface(x=x, y=y, z=value, spacing=spacing,region=region)
    fig.grdimage(frame=frame,grid = grid,projection=projection, region=region_oblique,perspective=perspective) # nan_transparent may work

    # tectonic setting
    fig.coast(region=region_oblique, frame=frame, projection=projection, perspective=perspective, shorelines="1p,black")    # coastlines
    (SAFy,SAFx) = line_read("tectonics/SAF")
    fig.plot(x = SAFx, y = SAFy, pen = '3.0p,black', perspective = perspective) # SAF
    if idepth == 0:
        fig.text(text = "SMB", x = -120.45 , y = 35.0, font = "16p,Helvetica-Bold,black", angle = 150, fill = "lightblue", perspective = perspective)   # SMB
        fig.text(text = "FT", x = -120.6 , y = 36.50, font = "16p,Helvetica-Bold,black", angle = 150, fill = "lightblue", perspective = perspective)    # Franciscan terrane
        fig.text(text = "ST", x = -121.1 , y = 36.0, font = "16p,Helvetica-Bold,black", angle = 150, fill = "lightblue", perspective = perspective) # Salinian terrane
        fig.text(text = "TR", x = -119.30 , y = 34.70, font = "16p,Helvetica-Bold,black", angle = 150, fill = "lightblue", perspective = perspective)     # Coast Ranges

    # depth label
    fig.text(text="%d km"%(depth), x = -119.8 , y = 34.0, font = "16p,Helvetica-Bold,black", angle = 180, fill = "white", perspective = perspective)     # Coast Ranges

    
    # --------- plot anisotropy ------------
    fig.shift_origin(yshift=-12)

    fig.basemap(region=region_oblique, frame=frame, projection=projection, perspective=perspective)
    pygmt.makecpt(cmap="cool", series=[0, 0.08], background=True)
    value = epsilon_inv[:,2]
    grid = pygmt.surface(x=x, y=y, z=value, spacing=spacing,region=region)
    fig.grdimage(frame=frame,grid = grid,projection=projection, region=region_oblique,perspective=perspective) # nan_transparent may work

    # tectonic setting
    fig.coast(region=region_oblique, frame=frame, projection=projection, perspective=perspective, shorelines="1p,black")    # coastlines
    (line_y,line_x) = line_read("tectonics/SAF_creeping")
    fig.plot(x = line_x, y = line_y, pen = '3.0p,black',perspective = perspective)
    (line_y,line_x) = line_read("tectonics/SAF_transition")
    fig.plot(x = line_x, y = line_y, pen = '3.0p,red',perspective = perspective)
    (line_y,line_x) = line_read("tectonics/SAF_locked")
    fig.plot(x = line_x, y = line_y, pen = '3.0p,blue',perspective = perspective)

    # anisotropy
    if len(ani) > 0:
        # rotate back to original coordinates
        x = ani[:,0]; y = ani[:,1]
        y,x = ffd.rtp_rotation_reverse(y,x,central_lat,central_lon,rotation_angle)
        ani[:,0] = x; ani[:,1] = y; # no need to modify the angle, because the porjection angle and rotate angle are the same
        fig.plot(ani, style='j', fill='yellow1', pen='0.5p,black',perspective=perspective)

    fig.shift_origin(xshift=6,yshift=12)


fig.show()

fig.savefig("img/imaging_result.png")


