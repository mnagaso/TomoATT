#ifndef GRID_WRAPPER_CUH
#define GRID_WRAPPER_CUH

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cuda.h>

#include <vector>

//#include "config.h"
#include "cuda_constants.cuh"
#include "cuda_utils.cuh"

// structure for storing grid information on device
typedef struct Grid_on_device {

} Grid_on_device;


void cuda_allocate_memory_for_grid_1st(Grid_on_device* grid_dv, int const& loc_I, int const& loc_J, int const& loc_K,
                CUSTOMREAL const& dp, CUSTOMREAL const& dt, CUSTOMREAL const& dr, \
                std::vector<std::vector<int*>>        const& vv_i__j__k__, \
                std::vector<std::vector<int*>>        const& vv_ip1j__k__, \
                std::vector<std::vector<int*>>        const& vv_im1j__k__, \
                std::vector<std::vector<int*>>        const& vv_i__jp1k__, \
                std::vector<std::vector<int*>>        const& vv_i__jm1k__, \
                std::vector<std::vector<int*>>        const& vv_i__j__kp1, \
                std::vector<std::vector<int*>>        const& vv_i__j__km1, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_fac_a, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_fac_b, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_fac_c, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_fac_f, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_T0v, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_T0r, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_T0t, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_T0p, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_fun, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_change);

void cuda_allocate_memory_for_grid_3rd(Grid_on_device* grid_dv, int const& loc_I, int const& loc_J, int const& loc_K,
                CUSTOMREAL const& dp, CUSTOMREAL const& dt, CUSTOMREAL const& dr, \
                std::vector<std::vector<int*>>        const& vv_i__j__k__, \
                std::vector<std::vector<int*>>        const& vv_ip1j__k__, \
                std::vector<std::vector<int*>>        const& vv_im1j__k__, \
                std::vector<std::vector<int*>>        const& vv_i__jp1k__, \
                std::vector<std::vector<int*>>        const& vv_i__jm1k__, \
                std::vector<std::vector<int*>>        const& vv_i__j__kp1, \
                std::vector<std::vector<int*>>        const& vv_i__j__km1, \
                std::vector<std::vector<int*>>        const& vv_ip2j__k__, \
                std::vector<std::vector<int*>>        const& vv_im2j__k__, \
                std::vector<std::vector<int*>>        const& vv_i__jp2k__, \
                std::vector<std::vector<int*>>        const& vv_i__jm2k__, \
                std::vector<std::vector<int*>>        const& vv_i__j__kp2, \
                std::vector<std::vector<int*>>        const& vv_i__j__km2, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_fac_a, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_fac_b, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_fac_c, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_fac_f, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_T0v, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_T0r, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_T0t, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_T0p, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_fun, \
                std::vector<std::vector<CUSTOMREAL*>> const& vv_change);


#endif // GRID_WRAPPER_CUH