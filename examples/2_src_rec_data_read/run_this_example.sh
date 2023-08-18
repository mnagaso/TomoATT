#!/bin/bash

python make_test_model.py

# run for no swap serial
mpirun -n 1 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_no_swap_serial.yml

# run for no swap parallel
mpirun -n 2 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_no_swap_parallel.yml

# run for swap serial
mpirun -n 1 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_swap_serial.yml

# run for swap parallel
mpirun -n 2 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_swap_parallel.yml

#python compare_src_rec.py

# then final model can be plot by e.g. check_3d_out.ipynb
#paraview OUTPUT_FILES/out_data_sim.xmf
