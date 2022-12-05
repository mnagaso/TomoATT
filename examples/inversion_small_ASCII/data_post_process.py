# This script shows how to recompose the output data for further analysis
# As the result data is collectively written out by each process, the data is composed on subgrid basis.
# the data needs to be reindexed and reassembled to a global grid.

import numpy as np
import sys
sys.path.append("../../utils/")

from tomoatt_data_retrieval import get_data_from_ascii


#
# plot calculated travel time field

# number of grids
nx = 50
ny = 50
nz = 10

# number of subgrid
ndiv_x = 2
ndiv_y = 2
ndiv_z = 1

# read data
_src_id = 0
_inv_id = 99

# filling 0 for digit
inv_id = str(_inv_id).zfill(4)
src_id = str(_src_id).zfill(4)

# file name
fname_data = "OUTPUT_FILES/fun_inv_{}.dat".format(inv_id)
fname_grid = "OUTPUT_FILES/out_grid_ptr.dat"

data, r, t, p = get_data_from_ascii(fname_data, fname_grid, nz, ny, nx, ndiv_z, ndiv_y, ndiv_x, verbose=True)

print(data.shape)
print(r.shape)

# plot
import matplotlib.pyplot as plt

fig, axs = plt.subplots(3, 1, figsize=(15, 5))

ir_slice=5
axs[0].pcolormesh(p[ir_slice,0,:], t[ir_slice,:,0], data[ir_slice,:,:], shading='auto')
axs[0].set_title('lon-lat')
axs[0].set_xlabel('longitude')
axs[0].set_ylabel('latitude')
axs[0].set_aspect(1)

it_slice=15
axs[1].pcolormesh(p[0,it_slice,:], r[:,it_slice,0], data[:,it_slice,:], shading='auto')
axs[1].set_title('r-lon')
axs[1].set_xlabel('longitude')
axs[1].set_ylabel('latitude')
axs[1].set_aspect(0.02)

ip_slice=15
axs[2].pcolormesh(t[0,:,ip_slice], r[:,0,ip_slice], data[:,:,ip_slice], shading='auto')
axs[2].set_title('r-lat')
axs[2].set_xlabel('longitude')
axs[2].set_ylabel('latitude')
axs[2].set_aspect(0.02)

plt.tight_layout()
plt.show()