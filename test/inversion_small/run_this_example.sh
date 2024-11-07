#!/bin/bash

python make_test_model.py

# get number of sweep parallel processes from input variable
nproc_sweep=$1
# get number of processes for domain decomposition
nproc_dd=$2

# modify the input_params_pre.yml file
#  nproc_sub: (integer) # number of processors for sweep parallelization (parallel the fast sweep method)
sed -i "s/nproc_sub: [0-9]\+/nproc_sub: $nproc_sweep/g" input_params_pre.yml
#   ndiv_rtp: [1, 1, 1] # number of subdivision on each direction (parallel the computional domain)
sed -i "s/ndiv_rtp: \[[0-9]\+, [0-9]\+, [0-9]\+\]/ndiv_rtp: \[$nproc_dd, $nproc_dd, $nproc_dd\]/g" input_params_pre.yml

#  nproc_sub: (integer) # number of processors for sweep parallelization (parallel the fast sweep method)
sed -i "s/nproc_sub: [0-9]\+/nproc_sub: $nproc_sweep/g" input_params.yml
#   ndiv_rtp: [1, 1, 1] # number of subdivision on each direction (parallel the computional domain)
sed -i "s/ndiv_rtp: \[[0-9]\+, [0-9]\+, [0-9]\+\]/ndiv_rtp: \[$nproc_dd, $nproc_dd, $nproc_dd\]/g" input_params.yml

# run for preparing true travel times
mpirun --oversubscribe -n 2 ../../build/bin/TOMOATT -i input_params_pre.yml

# run for inversion
mpirun --oversubscribe -n 2 ../../build/bin/TOMOATT -i input_params.yml
#mpirun -n 8 ../../build/bin/TOMOATT -i input_params_flexible_inv_grid.yml

# then final model can be plot by e.g. check_3d_out.ipynb
#paraview OUTPUT_FILES/out_data_sim.xmf
