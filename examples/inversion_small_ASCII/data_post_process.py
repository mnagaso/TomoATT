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
_inv_id = 25

# filling 0 for digit
inv_id = str(_inv_id).zfill(4)
src_id = str(_src_id).zfill(4)

# file name
#fname_data = "OUTPUT_FILES/T_res_inv_{}_src_{}.dat".format(inv_id, src_id)
fname_data = "OUTPUT_FILES/fun_inv_{}.dat".format(inv_id)

data = get_data_from_ascii(fname_data, nz, ny, nx, ndiv_z, ndiv_y, ndiv_x, verbose=True)


# plot
import matplotlib.pyplot as plt

fig, axs = plt.subplots(3, 1, figsize=(15, 5))

axs[0].imshow(data[5,:,:], cmap='jet', interpolation='bilinear')
axs[0].set_title('lon-lat')

axs[1].imshow(data[:,15,:], cmap='jet', interpolation='bilinear')
axs[1].set_title('r-lon')

axs[2].imshow(data[:,:,15], cmap='jet', interpolation='bilinear')
axs[2].set_title('r-lat')


plt.show()