# %%
# add ../../utils to the path
import sys
sys.path.append('../../utils')

from src_rec_file_helper import *

# read file
fpath_true = "./src_rec_test_out.dat"
event_list_true = read_src_rec_file(fpath_true)

fpath_mod = "./src_rec_test_out_modified.dat"
event_list_mod = read_src_rec_file(fpath_mod)

fpath_res = "./src_rec_test_out_modified_out.dat"
event_list_res = read_src_rec_file(fpath_res)

# %%
#plot
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(1, 1)
for event_true, event_mod, event_res in zip(event_list_true, event_list_mod, event_list_res):

    # true
    ax.scatter(event_true.lat, event_true.lon, c='b', marker='o', label='true', s=10)
    # modified
    ax.scatter(event_mod.lat, event_mod.lon, c='r', marker='o', label='modified')
    # result
    ax.scatter(event_res.lat, event_res.lon, c='g', marker='o', label='result', alpha=0.5)

ax.set_title("Earthquake Relocation")
ax.set_ylabel("latitude")
ax.set_xlabel("longitude")


# %%


# %%



