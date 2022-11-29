# %%
# add ../../utils to the path
import sys
sys.path.append('../../utils')

from src_rec_file_helper import *

# read file
fpath = "./OUTPUT_FILES/src_rec_file_forward.dat"
df_events, df_recs = read_src_rec_file(fpath, no_epi_dist=True)

# %%
df_events

# %%
# copy event list
import copy
df_events_ori = copy.deepcopy(df_events)

# %%
# modify source location
# slightly shift the source location towards the center of the domain
import math

tt1=(30.0-1.5)
tt2=(50.0+1.5)
pp1=(15.0-1.5)
pp2=(40.0+1.5)

center_lon = (pp1+pp2)/2
center_lat = (tt1+tt2)/2

# shift amount is 5 %
factor_shift =  0.05

# shift direction vector
for i in range(len(df_events)):
    v_lon = center_lon - df_events.loc[i, "lon"]
    v_lat = center_lat - df_events.loc[i, "lat"]

    df_events.loc[i, "lon"] += v_lon * factor_shift
    df_events.loc[i, "lat"] += v_lat * factor_shift


# %%
# plot original and modified source location
import matplotlib.pyplot as plt

for i in range(len(df_events_ori)):
    plt.plot(df_events_ori.loc[i, "lon"], df_events_ori.loc[i, "lat"], color="red", marker="o", markersize=1)

for i in range(len(df_events)):
    plt.plot(df_events.loc[i, "lon"], df_events.loc[i, "lat"], color="blue", marker="o", markersize=1)

plt.xlabel("Longitude")
plt.ylabel("Latitude")

plt.show()

# %%
# write out
fpath="src_rec_test_out_modified.dat"
write_src_rec_file(df_events, df_recs, fpath)

# %%



