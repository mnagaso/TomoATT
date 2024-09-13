#!/bin/bash

python make_test_model.py

# run for preparing true travel times
mpirun -n 2 ../../build/bin/TOMOATT -i input_params_pre.yml

# run for inversion
mpirun -n 2 ../../build/bin/TOMOATT -i input_params.yml
#mpirun -n 8 ../../build/bin/TOMOATT -i input_params_flexible_inv_grid.yml

# then final model can be plot by e.g. check_3d_out.ipynb
#paraview OUTPUT_FILES/out_data_sim.xmf
