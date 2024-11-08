#!/bin/bash

python make_test_model.py

# get number of sweep parallel processes from input variable (default 1)
nproc_sweep=${1:-1}
# get number of processes for domain decomposition
nproc_dd=${2:-1}
# calculate the total number of processes
nproc_total=$((nproc_sweep*nproc_dd*nproc_dd*nproc_dd))

# echo the number of processes
echo "Number of sweep parallel processes: $nproc_sweep"
echo "Number of processes for domain decomposition: $nproc_dd"
echo "Total number of processes: $nproc_total"

# modify the input_params_pre.yml file
#  nproc_sub: (integer) # number of processors for sweep parallelization (parallel the fast sweep method)
pta setpar input_params_pre.yml parallel.nproc_sub ${nproc_sweep}
#   ndiv_rtp: [1, 1, 1] # number of subdivision on each direction (parallel the computional domain)
pta setpar input_params_pre.yml parallel.ndiv_rtp ${nproc_dd},${nproc_dd},${nproc_dd}

#  nproc_sub: (integer) # number of processors for sweep parallelization (parallel the fast sweep method)
pta setpar input_params.yml parallel.nproc_sub ${nproc_sweep}
#   ndiv_rtp: [1, 1, 1] # number of subdivision on each direction (parallel the computional domain)
pta setpar input_params.yml parallel.ndiv_rtp ${nproc_dd},${nproc_dd},${nproc_dd}

# run for preparing true travel times
mpirun --oversubscribe -n $nproc_total ../../build/bin/TOMOATT -i input_params_pre.yml

# run for inversion
mpirun --oversubscribe -n $nproc_total ../../build/bin/TOMOATT -i input_params.yml
#mpirun -n 8 ../../build/bin/TOMOATT -i input_params_flexible_inv_grid.yml

# then final model can be plot by e.g. check_3d_out.ipynb
#paraview OUTPUT_FILES/out_data_sim.xmf
