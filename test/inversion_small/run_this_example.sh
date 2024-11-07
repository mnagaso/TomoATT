#!/bin/bash

#python make_test_model.py

# get number of sweep parallel processes from input variable
nproc_sweep=$1
# get number of processes for domain decomposition
nproc_dd=$2
# calculate the total number of processes
nproc_total=$((nproc_sweep*nproc_dd*nproc_dd*nproc_dd))

# echo the number of processes
echo "Number of sweep parallel processes: $nproc_sweep"
echo "Number of processes for domain decomposition: $nproc_dd"
echo "Total number of processes: $nproc_total"

# modify the input_params_pre.yml file
#  nproc_sub: (integer) # number of processors for sweep parallelization (parallel the fast sweep method)
sed -i "s/^\([[:space:]]*nproc_sub: *\)[0-9]*/\1${nproc_sweep}/" input_params_pre.yml
#   ndiv_rtp: [1, 1, 1] # number of subdivision on each direction (parallel the computional domain)
sed -i "s/^\([[:space:]]*ndiv_rtp: *\)\[.*\]/\1[${nproc_dd}, ${nproc_dd}, ${nproc_dd}]/" input_params_pre.yml

#  nproc_sub: (integer) # number of processors for sweep parallelization (parallel the fast sweep method)
sed -i "s/^\([[:space:]]*nproc_sub: *\)[0-9]*/\1${nproc_sweep}/" input_params.yml
#   ndiv_rtp: [1, 1, 1] # number of subdivision on each direction (parallel the computional domain)
sed -i "s/^\([[:space:]]*ndiv_rtp: *\)\[.*\]/\1[${nproc_dd}, ${nproc_dd}, ${nproc_dd}]/" input_params.yml

# run for preparing true travel times
mpirun --oversubscribe -n $nproc_total ../../build/bin/TOMOATT -i input_params_pre.yml

# run for inversion
mpirun --oversubscribe -n $nproc_total ../../build/bin/TOMOATT -i input_params.yml
#mpirun -n 8 ../../build/bin/TOMOATT -i input_params_flexible_inv_grid.yml

# then final model can be plot by e.g. check_3d_out.ipynb
#paraview OUTPUT_FILES/out_data_sim.xmf
