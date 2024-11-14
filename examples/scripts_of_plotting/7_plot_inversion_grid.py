# %%
import pygmt
pygmt.config(FONT="16p", IO_SEGMENT_MARKER="<<<")


# %%
import sys
sys.path.append('../utils')
import functions_for_data as ffd

# %%
# read inversion grid file

inv_grid_vel, inv_grid_ani = ffd.read_inversion_grid_file("OUTPUT_FILES")

print("inversion grid for velocity: ", inv_grid_vel.shape)
print("inversion grid for anisotropy: ", inv_grid_vel.shape)

Nset    = inv_grid_vel.shape[0]
Ngrid   = inv_grid_vel.shape[1]

colorlist = ["green","blue","red","purple","orange","yellow","black","gray","pink","cyan"]

# %%
# plot earthquakes and locations

fig = pygmt.Figure()

pygmt.makecpt(cmap="jet", series=[0, 40], background=True, reverse=True)    # colorbar

# -------- horizontal view (x-y) --------
fig.basemap(region=[0,2,0,2], frame=["xa1","ya1","NsWe"], projection="M10c")   # base map
# plot inversion grid
for igrid in range(Nset):
    x = inv_grid_vel[igrid,:,0]
    y = inv_grid_vel[igrid,:,1]
    fig.plot(x=x, y=y, style="c0.1c", fill=colorlist[igrid])

# -------- vertical view (x-z) --------
fig.shift_origin(xshift=0, yshift=-3)  
fig.basemap(region=[0,2,0,40], frame=["xa1","ya20+lDepth (km)","NsWe"], projection="X10c/-2c")   # base map
# plot inversion grid
for igrid in range(Nset):
    x = inv_grid_vel[igrid,:,0]
    y = inv_grid_vel[igrid,:,2]
    fig.plot(x=x, y=y, style="c0.1c", fill=colorlist[igrid])


# -------- vertical view (z-y) --------
fig.shift_origin(xshift=11, yshift=3)  
fig.basemap(region=[0,40,0,2], frame=["xa20+lDepth (km)","ya1","NsWe"], projection="X2c/10c")   # base map
# plot inversion grid
for igrid in range(Nset):
    x = inv_grid_vel[igrid,:,2]
    y = inv_grid_vel[igrid,:,1]
    fig.plot(x=x, y=y, style="c0.1c", fill=colorlist[igrid])


fig.savefig("img/7_inversion_grid.png")


