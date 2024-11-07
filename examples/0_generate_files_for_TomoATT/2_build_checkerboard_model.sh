#!/bin/bash

# Create a checkerboard model using the checkerboard model generator "pta create_checkerboard" from the PyTomoATT package
# The checkerboard model is created with the following parameters:
# -n4/4/2: number of anomalies in x/y/z direction for velocity
# -a4/4/2/60: number of anomalies in x/y/z direction for anisotropy, 60 is the fast velocity direction
# -p0.05/0.05: pertubation of anomalies for velocity and anisotropy, respectively
# -z0/40: depth of the model, limit the depth of the anormlies from 0 to 40 km
# -i2_models/model_init_N61_61_61.h5: initial model
# -o2_models/model_ckb_N61_61_61.h5: output model
# ./3_input_params/0_input_params_forward_simulation.yaml: input parameter file for forward simulation
pta create_checkerboard -n4/4/2 -a4/4/2/60 -p0.05/0.05 -z0/40 -i2_models/model_init_N61_61_61.h5 -o2_models/model_ckb_N61_61_61.h5 ./3_input_params/0_input_params_forward_simulation.yaml