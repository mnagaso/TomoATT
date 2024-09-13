import pygmt
from pygmt.clib import Session
with pygmt.clib.Session() as session:
    session.call_module('gmtset', 'FONT 16p')
pygmt.config(IO_SEGMENT_MARKER="<<<")

import numpy as np

def plot_map(x,y,value,dx,dy,
            contour = False,
            levels = None,
            anisotropy = None,
            station = None,
            earthquake = None,
            fname = "plot_map.png", 
            region = None,
            fig_size = None,
            axis_label = ["X","Y"],
            colorbar = "value", 
            cmap = "seis",
            cpt_range = None,
            y_reverse = False):
    """
    Plot the map of data (X,Y,VALUE) with pygmt
    x: 1D array of x coordinates
    y: 1D array of y coordinates
    value: 1D array of value
    dx: the spacing of x coordinates
    dy: the spacing of y coordinates
    contour: whether to plot contour, default is False
    levels: the levels of contour, default is (max(value) - min(value))/5
    anisotropy: the anisotropy of the data, default is None. If anisotropy is not None, fast velocity direction will be plotted
    station: the station of the data, default is None. If station is not None, station will be plotted
    earthquake: the earthquake of the data, default is None. If earthquake is not None, earthquake will be plotted
    fname: the name of the output file
    region: [x_min, x_max, y_min, y_max], the region of the map, default is [np.min(x),np.max(x),np.min(y),np.max(y)]
    fig_size: [x_size, y_size,], the size of figure, default is 10
    axis_label: the label of x and y axis, default is "X" and "Y"
    colorbar: the label of colorbar, default is "value"
    cpt: the color palette table, default is "seis"
    cpt_range: the range of colorbar, default is [min_value, max_value]
    y_reverse: whether to reverse the y axis, default is False
    """

    if cpt_range is not None:
        cpt_gap = (cpt_range[1] - cpt_range[0])/2
    else:
        min_value = np.min(value)   
        max_value = np.max(value)
        cpt_range = [min_value,max_value]
        cpt_gap = (max_value - min_value)/2 * 0.95
        print(cpt_gap)

    if levels is None:
        min_value = np.min(value)   
        max_value = np.max(value)
        levels = (max_value - min_value)/5

    if station is not None:
        station_x = station[0,:]
        station_y = station[1,:]
    
    if earthquake is not None:
        earthquake_x = earthquake[0,:]
        earthquake_y = earthquake[1,:]

    # set the region and projection
    if region is None:
        region = [np.min(x),np.max(x),np.min(y),np.max(y)]
    else:
        region = region

    if fig_size is None:
        x_size = 10
        y_size = 10* (np.max(y) - np.min(y))/(np.max(x) - np.min(x))
    else:
        x_size = fig_size[0]
        y_size = fig_size[1]
    
    if y_reverse:
        projection = "X%fc/-%fc"%(x_size,y_size)
    else:
        projection = "X%fc/%fc"%(x_size,y_size)

    x_label_gap = (region[1] - region[0])/4
    y_label_gap = (region[3] - region[2])/4
    frame = ["xa%f+l%s"%(x_label_gap, axis_label[0]), "ya%f+l%s"%(y_label_gap, axis_label[1]), "nSWe"]

    # create figure
    fig = pygmt.Figure()
    pygmt.makecpt(cmap=cmap, series=cpt_range, background=True, reverse=False)
    # grid = pygmt.xyz2grd(x=x, y=y, z=value, spacing = "%f/%f"%(dx,dy), region=region,)
    grid = pygmt.surface(x=x, y=y, z=value, spacing = "%f/%f"%(dx,dy), region=region,)

    # plot map
    fig.grdimage(frame=frame,grid = grid,projection=projection, region=region) # nan_transparent may work

    # plot earthquake
    if earthquake is not None:
        fig.plot(x=earthquake_x, y=earthquake_y, style='c0.3c', fill='red', pen='0.1p,white',label='Earthquake')

    # plot station
    if station is not None:
        fig.plot(x=station_x, y=station_y, style='t0.3c', fill='blue', pen='0.5p,white', label='Station')    

    # plot contour
    if contour:
        fig.contour(x=x, y=y, z=value, levels=5, pen="0.5p,black", annotation="5+f12p")

    # plot aniostropy
    if anisotropy is not None:
        fig.plot(anisotropy, style='j', fill='yellow1', pen='0.3p,black')

    fig.shift_origin(xshift= 1, yshift=-2)
    fig.colorbar(frame = ["a%f"%(cpt_gap),"y+l%s"%(colorbar)], position="+w4c/0.3c+h") # +e,默认是双箭头，f表示forward，b表示background ，w表示长宽，h表示水平
    fig.shift_origin(xshift=-1, yshift=+2)

    fig.show()
    fig.savefig(fname)

    return fig



def plot_profile(x, y, value,
                 fname,
                 region,
                 projection,
                 frame,
                 spacing,
                 cpt_range,
                 cpt_gap = None,
                 colorbar = "value",
                 cmap = "seis"
                          ):
    """
    Plot the profile of data (X,Y,VALUE) with pygmt
    x: 1D array of x coordinates
    y: 1D array of y coordinates
    value: 1D array of value
    dx: the spacing of x coordinates
    dy: the spacing of y coordinates
    """

    if cpt_gap is None:
        cpt_gap = (cpt_range[1] - cpt_range[0])/2

    fig = pygmt.Figure()
    pygmt.makecpt(cmap="seis", series=cpt_range, background=True, reverse=False)
    grid = pygmt.surface(x=x, y=y, z=value, spacing=spacing, region=region,)

    fig.grdimage(frame=frame,grid = grid, projection=projection, region=region) # nan_transparent may work

    fig.shift_origin(xshift= 1, yshift=-2)
    fig.colorbar(frame = ["a%f"%(cpt_gap),"y+l%s"%(colorbar)], position="+w4c/0.3c+h") # +e,默认是双箭头，f表示forward，b表示background ，w表示长宽，h表示水平
    fig.shift_origin(xshift=-1, yshift=+2)

    fig.show()

    fig.savefig(fname)

    return fig

