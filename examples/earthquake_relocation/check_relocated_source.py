# %%
# add ../../utils to the path
import sys
sys.path.append('../../utils')

from src_rec_file_helper import *

# read file
fpath_true = "./OUTPUT_FILES/src_rec_file_forward.dat"
df_ev_true,_ = read_src_rec_file(fpath_true, no_epi_dist=True)

fpath_mod = "./src_rec_test_out_modified.dat"
df_ev_mod,_ = read_src_rec_file(fpath_mod, no_epi_dist=True)

fpath_res = "./OUTPUT_FILES/src_rec_file_src_reloc_syn.dat"
df_ev_res,_ = read_src_rec_file(fpath_res, no_epi_dist=True)

# %%
#plot
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(1, 1)

for i in range(len(df_ev_true)):
    if (i == 0):
        # true
        ax.scatter(df_ev_true.loc[i, "lat"], df_ev_true.loc[i, "lon"], c='b', marker='o', label='true', s=10)
        # modified
        ax.scatter(df_ev_mod.loc[i, "lat"], df_ev_mod.loc[i, "lon"], c='r', marker='o', label='modified')
        # result
        ax.scatter(df_ev_res.loc[i, "lat"], df_ev_res.loc[i, "lon"], c='g', marker='o', label='result', alpha=0.5)
    else:
        # no legend
        # true
        ax.scatter(df_ev_true.loc[i, "lat"], df_ev_true.loc[i, "lon"], c='b', marker='o', s=10)
        # modified
        ax.scatter(df_ev_mod.loc[i, "lat"], df_ev_mod.loc[i, "lon"], c='r', marker='o', )
        # result
        ax.scatter(df_ev_res.loc[i, "lat"], df_ev_res.loc[i, "lon"], c='g', marker='o', alpha=0.5)


ax.set_title("Earthquake Relocation")
ax.set_ylabel("latitude")
ax.set_xlabel("longitude")

ax.legend()

plt.show()

# %%


# %%



