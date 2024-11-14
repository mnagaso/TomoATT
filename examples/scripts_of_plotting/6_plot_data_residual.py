# %%
import sys
sys.path.append('../utils')
import functions_for_data as ffd


# %%
# synthetic and observational traveltime files in the initial and final models 

file_init_syn = "OUTPUT_FILES/src_rec_file_inv_0000_reloc_0000.dat"         # synthetic traveltime in the initial model
file_init_obs = "input_files/src_rec_file.dat"                              # observational traveltime in the initial model

file_final_syn = "OUTPUT_FILES/src_rec_file_inv_0009_reloc_0009.dat"        # synthetic traveltime in the final model
file_final_obs = "OUTPUT_FILES/src_rec_file_inv_0009_reloc_0009_obs.dat"    # observational traveltime in the final model


# %%
# from pytomoatt.src_rec import SrcRec
# init_syn = SrcRec.read(file_init_syn)
# init_obs = SrcRec.read(file_init_obs)

# final_syn = SrcRec.read(file_final_syn)
# final_obs = SrcRec.read(file_final_obs)

# %%
ev, st = ffd.read_src_rec_file(file_init_syn)
time_init_syn = ffd.data_dis_time(ev, st)[1]    # synthetic traveltime in the initial model

ev, st = ffd.read_src_rec_file(file_init_obs)
time_init_obs = ffd.data_dis_time(ev, st)[1]    # observational traveltime in the initial model

ev, st = ffd.read_src_rec_file(file_final_syn)
time_final_syn = ffd.data_dis_time(ev, st)[1]    # synthetic traveltime in the final model

ev, st = ffd.read_src_rec_file(file_final_obs)
time_final_obs = ffd.data_dis_time(ev, st)[1]    # observational traveltime in the final model

# %%
import os
try:
    os.mkdir("img")
except:
    pass

# %%
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

range_l = -1.5
range_r = 1.5
Nbar = 20

bins=np.linspace(range_l,range_r,Nbar)
error_init = time_init_syn - time_init_obs
error_final = time_final_syn - time_final_obs

tag1 = "initial mode"
tag2 = "final mode"

hist_init,  _, _ = ax.hist(error_init,bins=bins,histtype='step', edgecolor = "red", linewidth = 2,
        label = "%s: std = %5.3f s, mean = %5.3f s"%(tag1,np.std(error_init),np.mean(error_init)))

hist_final, _, _ = ax.hist(error_final,bins=bins,alpha = 0.5, color = "blue",
        label = "%s: std = %5.3f s, mean = %5.3f s"%(tag2,np.std(error_final),np.mean(error_final)))

print("residual for ",tag1," model is: ","mean: ",np.mean(error_init),"sd: ",np.std(error_init))
print("residual for ",tag2," model is: ","mean: ",np.mean(error_final),"sd: ",np.std(error_final))
ax.legend(fontsize=14)

ax.set_xlim(range_l - abs(range_l)*0.1,range_r + abs(range_r)*0.1)
ax.set_ylim(0,1.3*max(max(hist_init),max(hist_final)))

ax.tick_params(axis='x',labelsize=18)
ax.tick_params(axis='y',labelsize=18)
ax.set_ylabel('Count of data',fontsize=18)
ax.set_xlabel('Traveltime residuals (s)',fontsize=18)
ax.set_title("$t_{syn} - t_{obs}$",fontsize=18)
ax.grid()

plt.savefig("img/6_data_residual.png",dpi=300, bbox_inches='tight', edgecolor='w', facecolor='w' )
plt.show()


