# %%
import sys
sys.path.append('../utils')
import functions_for_data as ffd

# %%
# read objective function 

path = "OUTPUT_FILES"
full_curve, location_curve, model_curve = ffd.read_objective_function_file(path)

print("full_curve: ", full_curve.shape, ", the total objective function value during the inversion, including relocation and model update")
print("location_curve: ", location_curve.shape, ", the objective function value during the relocation step")
print("model_curve: ", model_curve.shape, ", the objective function value during the model update step")

print("The first index is iteraion number, the second index is the objective function value vector")


# %%
# (Option 1) objective function value
full_obj = full_curve[:,0]  
location_obj = location_curve[:,0]
model_obj = model_curve[:,0]

# (Option 2) objective function value for only traveltime
full_obj_tt = full_curve[:,1]
location_obj_tt = location_curve[:,1]
model_obj_tt = model_curve[:,1]

# (Option 3) objective function value for only common source differential arrival time
full_obj_cs = full_curve[:,2]
location_obj_cs = location_curve[:,2]
model_obj_cs = model_curve[:,2]

# (Option 4) objective function value for only common receiver differential arrival time
full_obj_cr = full_curve[:,3]
location_obj_cr = location_curve[:,3]
model_obj_cr = model_curve[:,3]

# (Option 5) objective function value for teleseismic differential arrival time
full_obj_tele = full_curve[:,4]
location_obj_tele = location_curve[:,4]
model_obj_tele = model_curve[:,4]

# (Option 6) mean value of all data residual
full_mean = full_curve[:,5]
location_mean = location_curve[:,5]
model_mean = model_curve[:,5]

# (Option 7) standard deviation of all data residual
full_std = full_curve[:,6]
location_std = location_curve[:,6]
model_std = model_curve[:,6]

# (Option 8) mean value of residuals of traveltime
full_mean_tt = full_curve[:,7]
location_mean_tt = location_curve[:,7]
model_mean_tt = model_curve[:,7]

# (Option 9) standard deviation of residuals of traveltime
full_std_tt = full_curve[:,8]
location_std_tt = location_curve[:,8]
model_std_tt = model_curve[:,8]

# (Option 10) mean value of residuals of common source differential arrival time
full_mean_cs = full_curve[:,9]
location_mean_cs = location_curve[:,9]
model_mean_cs = model_curve[:,9]

# (Option 11) standard deviation of residuals of common source differential arrival time
full_std_cs = full_curve[:,10]
location_std_cs = location_curve[:,10]
model_std_cs = model_curve[:,10]

# (Option 12) mean value of residuals of common receiver differential arrival time
full_mean_cr = full_curve[:,11]
location_mean_cr = location_curve[:,11]
model_mean_cr = model_curve[:,11]

# (Option 13) standard deviation of residuals of common receiver differential arrival time
full_std_cr = full_curve[:,12]
location_std_cr = location_curve[:,12]
model_std_cr = model_curve[:,12]

# (Option 14) mean value of residuals of teleseismic differential arrival time
full_mean_tele = full_curve[:,13]
location_mean_tele = location_curve[:,13]
model_mean_tele = model_curve[:,13]

# (Option 15) standard deviation of residuals of teleseismic differential arrival time
full_std_tele = full_curve[:,14]
location_std_tele = location_curve[:,14]
model_std_tele = model_curve[:,14]


# %%
import os
try:
    os.mkdir("img")
except:
    pass

# %%
# plot objective functin reduction

import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(10, 6))
ax = plt.subplot(1, 1, 1)

ax.plot(model_obj/np.max(model_obj), label='objective function', linewidth=2)
ax.set_xlim([-0.2, len(model_obj)-0.8])
ax.set_ylim([0, 1.1])
ax.grid()
ax.set_xlabel('Iteration number',fontsize=14)
ax.set_ylabel('Normalized value',fontsize=14)
ax.tick_params(axis='x', labelsize=14)  
ax.tick_params(axis='y', labelsize=14) 
ax.legend(fontsize=14)

plt.savefig('img/5_objective_function_reduction.png', dpi=300, bbox_inches='tight', edgecolor='w', facecolor='w')


