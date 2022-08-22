#ifndef GRID_WRAPPER_CUH
#define GRID_WRAPPER_CUH

#include <mpi.h>

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cuda.h>

#include <vector>

//#include "config.h"
#include "cuda_constants.cuh"
#include "cuda_utils.cuh"

//#include "mpi_funcs.h"

#pragma once
#include "grid.fwd.h"

//#include "grid.h"

// structure for storing grid information on device
typedef struct Grid_on_device {

    //
    // grid info
    //

    int ngrid_i_host, ngrid_j_host, ngrid_k_host;

    int loc_I_host, loc_J_host, loc_K_host;

    int loc_I_excl_ghost_host, loc_J_excl_ghost_host, loc_K_excl_ghost_host;
    int n_points_ijk;
    int n_points_ijk_excl_ghost;

    int n_grid_bound_i_host, n_grid_bound_j_host, n_grid_bound_k_host;
    int n_ghost_layers_host;
    CUSTOMREAL dr_host, dt_host, dp_host;

    int i_start_loc_host, j_start_loc_host, k_start_loc_host;

    int i_end_loc_host, j_end_loc_host, k_end_loc_host;

    // neighbors_id_dev
    int* neighbors_id_host;

    // flag if inversion mode
    bool inverse_flag_host;

    // number of subdomains
    int n_subdomains_host;

    CUSTOMREAL L1_dev, L1_host;


    // tau_old_loc
    CUSTOMREAL* tau_old_loc_dev;
    // tau_loc
    CUSTOMREAL* tau_loc_dev;
    // fac_a_loc
    CUSTOMREAL* fac_a_loc_dev;
    // fac_b_loc
    CUSTOMREAL* fac_b_loc_dev;
    // fac_c_loc
    CUSTOMREAL* fac_c_loc_dev;
    // fac_f_loc
    CUSTOMREAL* fac_f_loc_dev;
    // T0v_loc
    CUSTOMREAL* T0v_loc_dev;
    // T0r_loc
    CUSTOMREAL* T0r_loc_dev;
    // T0t_loc
    CUSTOMREAL* T0t_loc_dev;
    // T0p_loc
    CUSTOMREAL* T0p_loc_dev;
    // fun_loc
    CUSTOMREAL* fun_loc_dev;

    //
    // arrays for mpi boundary communication
    // memcopy between device and host is done only for those arrays
    // to minimise the amount of memory transfer

    // bin_s
    CUSTOMREAL* bin_s_dev;
    // bin_r
    CUSTOMREAL* bin_r_dev;
    // bip_s
    CUSTOMREAL* bip_s_dev;
    // bip_r
    CUSTOMREAL* bip_r_dev;
    // bjn_s
    CUSTOMREAL* bjn_s_dev;
    // bjn_r
    CUSTOMREAL* bjn_r_dev;
    // bjp_s
    CUSTOMREAL* bjp_s_dev;
    // bjp_r
    CUSTOMREAL* bjp_r_dev;
    // bkn_s
    CUSTOMREAL* bkn_s_dev;
    // bkn_r
    CUSTOMREAL* bkn_r_dev;
    // bkp_s
    CUSTOMREAL* bkp_s_dev;
    // bkp_r
    CUSTOMREAL* bkp_r_dev;

    //
    // for adjoint run
    //
    // Tadj_loc
    CUSTOMREAL* Tadj_loc_dev;
    // t_loc_1d
    CUSTOMREAL* t_loc_1d_dev;
    // r_loc_1d
    CUSTOMREAL* r_loc_1d_dev;
    // zeta_loc
    CUSTOMREAL* zeta_loc_dev;
    // xi_loc
    CUSTOMREAL* xi_loc_dev;
    // eta_loc
    CUSTOMREAL* eta_loc_dev;
    // T_loc
    CUSTOMREAL* T_loc_dev;

    // is changed
    bool* is_changed_dev;

    // tmp array unified
    CUSTOMREAL* L1_arr_unified;
    CUSTOMREAL* L1_tmp_dev;

    // streams
    // use two streams for computation and memory transfer
    // for overlapping those processes

    // stream for compute
    cudaStream_t compute_stream;
    // stream for memory copy
    cudaStream_t copy_stream;
    // array of streams
    cudaStream_t* level_streams;

    // store MPI_comm for communicating between subdomains
    MPI_Comm inter_sub_comm_host;
    // mpi rank
    int myrank_host;

    // level set collection
    int* i_id_all_level;
    int* j_id_all_level;
    int* k_id_all_level;
    int* n_points_per_level_dev;
    int* n_points_per_level_host; // host copy of n_points_per_level

    // other constants
    int st_level_host, ed_level_host;
    int n_max_points_level_host; // max number of points among all levels
    int n_total_points_host; // total number of points in the grid

    // stencil order
    int stencil_order_host;

    // thead and grid for sweeping
    dim3 grid_sweep_host, threads_sweep_host;
    // if use minimum unified kernel for all level
    bool use_unified_kernel_host = true;

    // thread and grid for computing L1
    dim3 grid_L1_host, threads_L1_host;
    int n_blocks_L1_host;

    // thread and grid for 3d full increment
    dim3 grid_3d_full_incr_host, threads_3d_full_incr_host;


} Grid_on_device;

// allocate memory for grid arrays on device
void cuda_allocate_memory_for_grid(Grid_on_device*, \
    int ngrid_i_host, int ngrid_j_host, int ngrid_k_host, \
    int loc_I_host, \
    int loc_J_host, \
    int loc_K_host, \
    int loc_I_excl_ghost_host, \
    int loc_J_excl_ghost_host, \
    int loc_K_excl_ghost_host, \
    int n_grid_bound_i_host, \
    int n_grid_bound_j_host, \
    int n_grid_bound_k_host, \
    int n_ghost_layers_host, \
    bool inverse_flag, \
    int i_start_loc_host, \
    int j_start_loc_host, \
    int k_start_loc_host, \
    int i_end_loc_host, \
    int j_end_loc_host, \
    int k_end_loc_host, \
    int n_subdomains_tmp, \
    int myrank_host);

// copy data from host to device
void cuda_data_transfer_to_device(Grid_on_device*, \
    int          loc_I_host, \
    int          loc_J_host, \
    int          loc_K_host, \
    int n_grid_bound_i_host, \
    int n_grid_bound_j_host, \
    int n_grid_bound_k_host, \
    int n_ghost_layers_host, \
    CUSTOMREAL      dr_host, \
    CUSTOMREAL      dt_host, \
    CUSTOMREAL      dp_host, \
    std::vector<int> neighbors_id_host, \
    CUSTOMREAL* tau_old_loc_host, \
    CUSTOMREAL*   tau_loc_host, \
    CUSTOMREAL* fac_a_loc_host, \
    CUSTOMREAL* fac_b_loc_host, \
    CUSTOMREAL* fac_c_loc_host, \
    CUSTOMREAL* fac_f_loc_host, \
    CUSTOMREAL*   T0v_loc_host, \
    CUSTOMREAL*   T0r_loc_host, \
    CUSTOMREAL*   T0t_loc_host, \
    CUSTOMREAL*   T0p_loc_host, \
    CUSTOMREAL*   fun_loc_host, \
    CUSTOMREAL*  Tadj_loc_host, \
    CUSTOMREAL*  t_loc_1d_host, \
    CUSTOMREAL*  r_loc_1d_host, \
    CUSTOMREAL*  zeta_loc_host, \
    CUSTOMREAL*    xi_loc_host, \
    CUSTOMREAL*   eta_loc_host, \
    CUSTOMREAL*     T_loc_host, \
    bool*      is_changed_host);


// free memory for grid arrays on device
void cuda_free_memory_for_grid(Grid_on_device*);

// prepare grid informations on gpu
void cuda_init_gpu_params(Grid_on_device*, \
                          MPI_Comm, \
                          std::vector<std::vector<std::vector<int>>> &, \
                          CUSTOMREAL, CUSTOMREAL, CUSTOMREAL);

// resert grid when target source changed
void cuda_reinitialize_grid(Grid_on_device*, \
                            CUSTOMREAL*, \
                            CUSTOMREAL*, \
                            CUSTOMREAL*, \
                            CUSTOMREAL*, \
                            CUSTOMREAL*, \
                            CUSTOMREAL*, \
                            CUSTOMREAL*, \
                            CUSTOMREAL*, \
                            CUSTOMREAL*, \
                            CUSTOMREAL*, \
                            bool*);

// send and receiv boundary data between each mpi
void cuda_send_recev_boundary_data(Grid_on_device*);
void cuda_prepare_boundary_data_to_send(Grid_on_device*);
void cuda_assign_received_data_to_ghost(Grid_on_device*);


#endif // GRID_WRAPPER_CUH