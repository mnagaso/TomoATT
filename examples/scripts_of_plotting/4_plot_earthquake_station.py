# %%
import pygmt
pygmt.config(FONT="16p", IO_SEGMENT_MARKER="<<<")


# %%
from pytomoatt.src_rec import SrcRec

# read src_rec_file
sr = SrcRec.read("input_files/src_rec_file.dat")

# get the coordinates of the stations and earthquakes
stations = sr.receivers[['stlo','stla','stel']].values.T
earthquakes = sr.sources[['evlo','evla','evdp']].values.T

print(stations.shape)
print(earthquakes.shape)

# %%
# plot earthquakes and locations

fig = pygmt.Figure()

pygmt.makecpt(cmap="jet", series=[0, 40], background=True, reverse=True)    # colorbar

# -------- horizontal view (x-y) --------
fig.basemap(region=[0,2,0,2], frame=["xa1","ya1","NsWe"], projection="M10c")   # base map
# earthquakes
fig.plot(x = earthquakes[0,:], y = earthquakes[1,:], cmap = True, style = "c0.1c", fill = earthquakes[2,:])
# stations
fig.plot(x = stations[0,:], y = stations[1,:], style = "t0.4c", fill = "blue", pen = "black", label = "Station")


# -------- vertical view (x-z) --------
fig.shift_origin(xshift=0, yshift=-3)  
fig.basemap(region=[0,2,0,40], frame=["xa1","ya20+lDepth (km)","NsWe"], projection="X10c/-2c")   # base map
# earthquakes
fig.plot(x = earthquakes[0,:], y = earthquakes[2,:], cmap = True, style = "c0.1c", fill = earthquakes[2,:])


# -------- vertical view (z-y) --------
fig.shift_origin(xshift=11, yshift=3)  
fig.basemap(region=[0,40,0,2], frame=["xa20+lDepth (km)","ya1","NsWe"], projection="X2c/10c")   # base map
# earthquakes
fig.plot(x = earthquakes[2,:], y = earthquakes[1,:], cmap = True, style = "c0.1c", fill = earthquakes[2,:])


# colorbar
fig.shift_origin(xshift=0, yshift=-1.5)  
fig.colorbar(frame = ["a20","x+lDepth (km)"], position="+e+w2c/0.3c+h") 

fig.savefig("img/4_earthquakes_and_stations.png")


