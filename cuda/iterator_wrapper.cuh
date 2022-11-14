#ifndef ITERATOR_WRAPPER_CUH
#define ITERATOR_WRAPPER_CUH

#include <cuda_runtime.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

#include <stdio.h>
#include "grid_wrapper.cuh"
#include "cuda_utils.cuh"


//void cuda_do_sweep_level_kernel_3rd();
//void cuda_do_sweep_level_kernel_1st();

void initialize_sweep_params(Grid_on_device*);
void cuda_run_iterate_forward(Grid_on_device*, int const&);


#endif // ITERATOR_WRAPPER_CUH