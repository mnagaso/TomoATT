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

    // parameters
    int loc_I_host, loc_J_host, loc_K_host;
    int n_nodes_total_host;
    int n_levels_host;

    // index storage
    int* n_nodes_on_levels, *n_nodes_on_levels_host;

    int* vv_i__j__k__0, *vv_i__j__k__1, *vv_i__j__k__2, *vv_i__j__k__3, *vv_i__j__k__4, *vv_i__j__k__5, *vv_i__j__k__6, *vv_i__j__k__7;
    int* vv_ip1j__k__0, *vv_ip1j__k__1, *vv_ip1j__k__2, *vv_ip1j__k__3, *vv_ip1j__k__4, *vv_ip1j__k__5, *vv_ip1j__k__6, *vv_ip1j__k__7;
    int* vv_im1j__k__0, *vv_im1j__k__1, *vv_im1j__k__2, *vv_im1j__k__3, *vv_im1j__k__4, *vv_im1j__k__5, *vv_im1j__k__6, *vv_im1j__k__7;
    int* vv_i__jp1k__0, *vv_i__jp1k__1, *vv_i__jp1k__2, *vv_i__jp1k__3, *vv_i__jp1k__4, *vv_i__jp1k__5, *vv_i__jp1k__6, *vv_i__jp1k__7;
    int* vv_i__jm1k__0, *vv_i__jm1k__1, *vv_i__jm1k__2, *vv_i__jm1k__3, *vv_i__jm1k__4, *vv_i__jm1k__5, *vv_i__jm1k__6, *vv_i__jm1k__7;
    int* vv_i__j__kp10, *vv_i__j__kp11, *vv_i__j__kp12, *vv_i__j__kp13, *vv_i__j__kp14, *vv_i__j__kp15, *vv_i__j__kp16, *vv_i__j__kp17;
    int* vv_i__j__km10, *vv_i__j__km11, *vv_i__j__km12, *vv_i__j__km13, *vv_i__j__km14, *vv_i__j__km15, *vv_i__j__km16, *vv_i__j__km17;
    int* vv_ip2j__k__0, *vv_ip2j__k__1, *vv_ip2j__k__2, *vv_ip2j__k__3, *vv_ip2j__k__4, *vv_ip2j__k__5, *vv_ip2j__k__6, *vv_ip2j__k__7;
    int* vv_im2j__k__0, *vv_im2j__k__1, *vv_im2j__k__2, *vv_im2j__k__3, *vv_im2j__k__4, *vv_im2j__k__5, *vv_im2j__k__6, *vv_im2j__k__7;
    int* vv_i__jp2k__0, *vv_i__jp2k__1, *vv_i__jp2k__2, *vv_i__jp2k__3, *vv_i__jp2k__4, *vv_i__jp2k__5, *vv_i__jp2k__6, *vv_i__jp2k__7;
    int* vv_i__jm2k__0, *vv_i__jm2k__1, *vv_i__jm2k__2, *vv_i__jm2k__3, *vv_i__jm2k__4, *vv_i__jm2k__5, *vv_i__jm2k__6, *vv_i__jm2k__7;
    int* vv_i__j__kp20, *vv_i__j__kp21, *vv_i__j__kp22, *vv_i__j__kp23, *vv_i__j__kp24, *vv_i__j__kp25, *vv_i__j__kp26, *vv_i__j__kp27;
    int* vv_i__j__km20, *vv_i__j__km21, *vv_i__j__km22, *vv_i__j__km23, *vv_i__j__km24, *vv_i__j__km25, *vv_i__j__km26, *vv_i__j__km27;

    // constants
    CUSTOMREAL* vv_fac_a_0, *vv_fac_a_1, *vv_fac_a_2, *vv_fac_a_3, *vv_fac_a_4, *vv_fac_a_5, *vv_fac_a_6, *vv_fac_a_7;
    CUSTOMREAL* vv_fac_b_0, *vv_fac_b_1, *vv_fac_b_2, *vv_fac_b_3, *vv_fac_b_4, *vv_fac_b_5, *vv_fac_b_6, *vv_fac_b_7;
    CUSTOMREAL* vv_fac_c_0, *vv_fac_c_1, *vv_fac_c_2, *vv_fac_c_3, *vv_fac_c_4, *vv_fac_c_5, *vv_fac_c_6, *vv_fac_c_7;
    CUSTOMREAL* vv_fac_f_0, *vv_fac_f_1, *vv_fac_f_2, *vv_fac_f_3, *vv_fac_f_4, *vv_fac_f_5, *vv_fac_f_6, *vv_fac_f_7;
    CUSTOMREAL* vv_T0v_0, *vv_T0v_1, *vv_T0v_2, *vv_T0v_3, *vv_T0v_4, *vv_T0v_5, *vv_T0v_6, *vv_T0v_7;
    CUSTOMREAL* vv_T0r_0, *vv_T0r_1, *vv_T0r_2, *vv_T0r_3, *vv_T0r_4, *vv_T0r_5, *vv_T0r_6, *vv_T0r_7;
    CUSTOMREAL* vv_T0t_0, *vv_T0t_1, *vv_T0t_2, *vv_T0t_3, *vv_T0t_4, *vv_T0t_5, *vv_T0t_6, *vv_T0t_7;
    CUSTOMREAL* vv_T0p_0, *vv_T0p_1, *vv_T0p_2, *vv_T0p_3, *vv_T0p_4, *vv_T0p_5, *vv_T0p_6, *vv_T0p_7;
    CUSTOMREAL* vv_fun_0, *vv_fun_1, *vv_fun_2, *vv_fun_3, *vv_fun_4, *vv_fun_5, *vv_fun_6, *vv_fun_7;
    CUSTOMREAL* vv_change_0, *vv_change_1, *vv_change_2, *vv_change_3, *vv_change_4, *vv_change_5, *vv_change_6, *vv_change_7;

    // temporary variables
    CUSTOMREAL* tau;

} Grid_on_device;


void cuda_initialize_grid_1st(std::vector< std::vector<int> >& ijk, Grid_on_device* grid_dv, int const& loc_I, int const& loc_J, int const& loc_K,
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

void cuda_initialize_grid_3rd(std::vector< std::vector<int> >& ijk, Grid_on_device* grid_dv, int const& loc_I, int const& loc_J, int const& loc_K,
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