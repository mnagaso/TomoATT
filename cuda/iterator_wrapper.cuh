#ifndef ITERATOR_WRAPPER_CUH
#define ITERATOR_WRAPPER_CUH

#include <cuda_runtime.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cooperative_groups.h>
//#include <cooperative_groups/reduce.h>

#include <stdio.h>
//#include "mpi_funcs.h"
#include "grid_wrapper.cuh"
#include "cuda_utils.cuh"

//namespace cg = cooperative_groups;


void cuda_run_iterate_forward(Grid_on_device*, CUSTOMREAL, int, int);

void cuda_copy_T_loc_tau_loc_to_host(CUSTOMREAL* h_T_loc, CUSTOMREAL* h_tau_loc, Grid_on_device* grid);



#endif // ITERATOR_WRAPPER_CUH