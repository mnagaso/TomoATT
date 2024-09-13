#!/bin/bash

# python make_test_model.py

# run for forward simulation in the mesh of 41*41*41
mpirun -n 1 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_N41.yaml

# run for forward simulation in the mesh of 61*61*61
mpirun -n 1 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_N61.yaml

# run for forward simulation in the mesh of 81*81*81
mpirun -n 1 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_N81.yaml

# run for forward simulation in the mesh of 121*121*121
mpirun -n 1 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_N121.yaml

# run for forward simulation in the mesh of 161*161*161
mpirun -n 1 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_N161.yaml

# this should be called from run_test.sh
#python compare_src_rec.py

# then final model can be plot by e.g. check_3d_out.ipynb
#paraview OUTPUT_FILES/out_data_sim.xmf
