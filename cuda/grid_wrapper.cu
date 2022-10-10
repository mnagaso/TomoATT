#include "grid_wrapper.cuh"

void cuda_allocate_memory_for_grid(Grid_on_device* grid_on_dv, \
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
    int myrank_host){

    // copy parameters to device
    printf("cuda_allocate_memory_for_grid\n");;

    grid_on_dv->myrank_host = myrank_host;

    grid_on_dv->ngrid_i_host = ngrid_i_host;
    grid_on_dv->ngrid_j_host = ngrid_j_host;
    grid_on_dv->ngrid_k_host = ngrid_k_host;

    grid_on_dv->loc_I_host = loc_I_host;
    grid_on_dv->loc_J_host = loc_J_host;
    grid_on_dv->loc_K_host = loc_K_host;
    grid_on_dv->loc_I_excl_ghost_host = loc_I_excl_ghost_host;
    grid_on_dv->loc_J_excl_ghost_host = loc_J_excl_ghost_host;
    grid_on_dv->loc_K_excl_ghost_host = loc_K_excl_ghost_host;
    grid_on_dv->n_grid_bound_i_host = n_grid_bound_i_host;
    grid_on_dv->n_grid_bound_j_host = n_grid_bound_j_host;
    grid_on_dv->n_grid_bound_k_host = n_grid_bound_k_host;
    grid_on_dv->n_ghost_layers_host = n_ghost_layers_host;
    grid_on_dv->st_level_host = 6;
    grid_on_dv->ed_level_host = loc_I_host+loc_J_host+loc_K_host-3;
    grid_on_dv->n_points_ijk = loc_I_host*loc_J_host*loc_K_host;
    grid_on_dv->n_points_ijk_excl_ghost = loc_I_excl_ghost_host*loc_J_excl_ghost_host*loc_K_excl_ghost_host;
    grid_on_dv->i_start_loc_host = i_start_loc_host;
    grid_on_dv->j_start_loc_host = j_start_loc_host;
    grid_on_dv->k_start_loc_host = k_start_loc_host;
    grid_on_dv->i_end_loc_host = i_end_loc_host;
    grid_on_dv->j_end_loc_host = j_end_loc_host;
    grid_on_dv->k_end_loc_host = k_end_loc_host;
    grid_on_dv->n_subdomains_host = n_subdomains_tmp; // in cuda_constants.cuh

    int n_total_points = loc_I_host * loc_J_host * loc_K_host;

    grid_on_dv->inverse_flag_host = inverse_flag;

    // allocate memory for grid on device
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->tau_old_loc_dev), n_total_points), 15);
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->tau_loc_dev),     n_total_points), 16);
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->fac_a_loc_dev),   n_total_points), 17);
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->fac_b_loc_dev),   n_total_points), 18);
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->fac_c_loc_dev),   n_total_points), 19);
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->fac_f_loc_dev),   n_total_points), 20);
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->T0v_loc_dev),     n_total_points), 21);
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->T0r_loc_dev),     n_total_points), 22);
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->T0t_loc_dev),     n_total_points), 23);
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->T0p_loc_dev),     n_total_points), 24);
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->fun_loc_dev),     n_total_points), 25);
    print_CUDA_error_if_any(allocate_memory_on_device_cv_pinned((void**) &(grid_on_dv->bin_s_dev),       n_ghost_layers_host*n_grid_bound_i_host), 26);
    print_CUDA_error_if_any(allocate_memory_on_device_cv_pinned((void**) &(grid_on_dv->bin_r_dev),       n_ghost_layers_host*n_grid_bound_i_host), 27);
    print_CUDA_error_if_any(allocate_memory_on_device_cv_pinned((void**) &(grid_on_dv->bip_s_dev),       n_ghost_layers_host*n_grid_bound_i_host), 28);
    print_CUDA_error_if_any(allocate_memory_on_device_cv_pinned((void**) &(grid_on_dv->bip_r_dev),       n_ghost_layers_host*n_grid_bound_i_host), 29);
    print_CUDA_error_if_any(allocate_memory_on_device_cv_pinned((void**) &(grid_on_dv->bjn_s_dev),       n_ghost_layers_host*n_grid_bound_j_host), 30);
    print_CUDA_error_if_any(allocate_memory_on_device_cv_pinned((void**) &(grid_on_dv->bjn_r_dev),       n_ghost_layers_host*n_grid_bound_j_host), 31);
    print_CUDA_error_if_any(allocate_memory_on_device_cv_pinned((void**) &(grid_on_dv->bjp_s_dev),       n_ghost_layers_host*n_grid_bound_j_host), 32);
    print_CUDA_error_if_any(allocate_memory_on_device_cv_pinned((void**) &(grid_on_dv->bjp_r_dev),       n_ghost_layers_host*n_grid_bound_j_host), 33);
    print_CUDA_error_if_any(allocate_memory_on_device_cv_pinned((void**) &(grid_on_dv->bkn_s_dev),       n_ghost_layers_host*n_grid_bound_k_host), 34);
    print_CUDA_error_if_any(allocate_memory_on_device_cv_pinned((void**) &(grid_on_dv->bkn_r_dev),       n_ghost_layers_host*n_grid_bound_k_host), 35);
    print_CUDA_error_if_any(allocate_memory_on_device_cv_pinned((void**) &(grid_on_dv->bkp_s_dev),       n_ghost_layers_host*n_grid_bound_k_host), 36);
    print_CUDA_error_if_any(allocate_memory_on_device_cv_pinned((void**) &(grid_on_dv->bkp_r_dev),       n_ghost_layers_host*n_grid_bound_k_host), 37);
//    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->bin_s_dev),       n_ghost_layers_host*n_grid_bound_i_host), 26);
//    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->bin_r_dev),       n_ghost_layers_host*n_grid_bound_i_host), 27);
//    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->bip_s_dev),       n_ghost_layers_host*n_grid_bound_i_host), 28);
//    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->bip_r_dev),       n_ghost_layers_host*n_grid_bound_i_host), 29);
//    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->bjn_s_dev),       n_ghost_layers_host*n_grid_bound_j_host), 30);
//    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->bjn_r_dev),       n_ghost_layers_host*n_grid_bound_j_host), 31);
//    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->bjp_s_dev),       n_ghost_layers_host*n_grid_bound_j_host), 32);
//    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->bjp_r_dev),       n_ghost_layers_host*n_grid_bound_j_host), 33);
//    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->bkn_s_dev),       n_ghost_layers_host*n_grid_bound_k_host), 34);
//    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->bkn_r_dev),       n_ghost_layers_host*n_grid_bound_k_host), 35);
//    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->bkp_s_dev),       n_ghost_layers_host*n_grid_bound_k_host), 36);
//    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->bkp_r_dev),       n_ghost_layers_host*n_grid_bound_k_host), 37);
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->t_loc_1d_dev),    loc_J_host), 38);
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->r_loc_1d_dev),    loc_K_host), 39);
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->zeta_loc_dev),    n_total_points), 40);
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->xi_loc_dev),      n_total_points), 41);
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->eta_loc_dev),     n_total_points), 42);
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**) &(grid_on_dv->T_loc_dev),       n_total_points), 43);
    print_CUDA_error_if_any(allocate_memory_on_device_bl((void**) &(grid_on_dv->is_changed_dev), n_total_points), 44);

    //print_CUDA_error_if_any(allocate_memory_on_device_i((void**)&(grid_on_dv->neighbors_id_dev), 6), 44);
    grid_on_dv->neighbors_id_host = new int[6];

    // array for inversion
    if (inverse_flag){
        print_CUDA_error_if_any(allocate_memory_on_device_cv((void**)&(grid_on_dv->Tadj_loc_dev),    n_total_points), 44);
    }

    // single value for storing L1
    grid_on_dv->L1_host = 0.0;
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_cv(&(grid_on_dv->L1_dev), &(grid_on_dv->L1_host), 1), 45);
    // array for reducing L1 values in each block
    int n_blocks_x, n_blocks_y;
    int n_total_pt = grid_on_dv->loc_I_host*grid_on_dv->loc_J_host*grid_on_dv->loc_K_host;
    int n_theads_per_block = CUDA_L1_BLOCK_SIZE;
    get_block_xy(ceil(n_total_pt/n_theads_per_block+0.5), &n_blocks_x, &n_blocks_y);
    int nblocks = n_blocks_x*n_blocks_y;
    printf("nblocks = %d\n", nblocks);

    print_CUDA_error_if_any(cudaMallocManaged((void**)&(grid_on_dv->L1_arr_unified), nblocks*sizeof(CUSTOMREAL)),191920);
    print_CUDA_error_if_any(cudaMalloc((void**)&(grid_on_dv->L1_tmp_dev), n_total_points*sizeof(CUSTOMREAL)),191921);

}

void cuda_free_memory_for_grid(Grid_on_device* grid_on_dv){
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->tau_old_loc_dev ), 10001);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->tau_loc_dev     ), 10002);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->fac_a_loc_dev   ), 10003);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->fac_b_loc_dev   ), 10004);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->fac_c_loc_dev   ), 10005);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->fac_f_loc_dev   ), 10006);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->T0v_loc_dev     ), 10007);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->T0r_loc_dev     ), 10008);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->T0t_loc_dev     ), 10009);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->T0p_loc_dev     ), 10010);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->fun_loc_dev     ), 10011);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->bin_s_dev       ), 10012);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->bin_r_dev       ), 10013);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->bip_s_dev       ), 10014);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->bip_r_dev       ), 10015);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->bjn_s_dev       ), 10016);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->bjn_r_dev       ), 10017);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->bjp_s_dev       ), 10018);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->bjp_r_dev       ), 10019);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->bkn_s_dev       ), 10020);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->bkn_r_dev       ), 10021);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->bkp_s_dev       ), 10022);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->bkp_r_dev       ), 10023);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->t_loc_1d_dev    ), 10024);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->r_loc_1d_dev    ), 10025);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->zeta_loc_dev    ), 10026);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->xi_loc_dev      ), 10027);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->eta_loc_dev     ), 10028);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->T_loc_dev       ), 10029);
    print_CUDA_error_if_any(deallocate_memory_on_device_bl(grid_on_dv->is_changed_dev  ), 10030);
    delete [] grid_on_dv->neighbors_id_host;

    if (grid_on_dv->inverse_flag_host){
        print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->Tadj_loc_dev),10032);
    }

    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_on_dv->i_id_all_level),10033);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_on_dv->j_id_all_level),10034);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_on_dv->k_id_all_level),10035);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_on_dv->n_points_per_level_dev),10036);
    delete [] grid_on_dv->n_points_per_level_host;

    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->L1_arr_unified),10038);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_on_dv->L1_tmp_dev),10039);
}


void cuda_data_transfer_to_device(Grid_on_device* grid_on_dv, \
    int loc_I_host, \
    int loc_J_host, \
    int loc_K_host, \
    int n_grid_bound_i_host, \
    int n_grid_bound_j_host, \
    int n_grid_bound_k_host, \
    int n_ghost_layers_host, \
    CUSTOMREAL dr_host, \
    CUSTOMREAL dt_host, \
    CUSTOMREAL dp_host, \
    std::vector<int> neighbors_id_host, \
    CUSTOMREAL* tau_old_loc_host, \
    CUSTOMREAL* tau_loc_host, \
    CUSTOMREAL* fac_a_loc_host, \
    CUSTOMREAL* fac_b_loc_host, \
    CUSTOMREAL* fac_c_loc_host, \
    CUSTOMREAL* fac_f_loc_host, \
    CUSTOMREAL* T0v_loc_host, \
    CUSTOMREAL* T0r_loc_host, \
    CUSTOMREAL* T0t_loc_host, \
    CUSTOMREAL* T0p_loc_host, \
    CUSTOMREAL* fun_loc_host, \
    CUSTOMREAL* Tadj_loc_host, \
    CUSTOMREAL* t_loc_1d_host, \
    CUSTOMREAL* r_loc_1d_host, \
    CUSTOMREAL* zeta_loc_host, \
    CUSTOMREAL* xi_loc_host, \
    CUSTOMREAL* eta_loc_host, \
    CUSTOMREAL* T_loc_host, \
    bool* is_changed_host){

    //printf("cuda_data_transfer_to_device\n");

    int n_total_points = loc_I_host * loc_J_host * loc_K_host;

    // copy host memory to device memory
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->tau_old_loc_dev, tau_old_loc_host, n_total_points), 46);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->tau_loc_dev,     tau_loc_host,     n_total_points), 47);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->fac_a_loc_dev,   fac_a_loc_host,   n_total_points), 48);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->fac_b_loc_dev,   fac_b_loc_host,   n_total_points), 49);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->fac_c_loc_dev,   fac_c_loc_host,   n_total_points), 50);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->fac_f_loc_dev,   fac_f_loc_host,   n_total_points), 51);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->T0v_loc_dev,     T0v_loc_host,     n_total_points), 52);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->T0r_loc_dev,     T0r_loc_host,     n_total_points), 53);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->T0t_loc_dev,     T0t_loc_host,     n_total_points), 54);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->T0p_loc_dev,     T0p_loc_host,     n_total_points), 55);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->fun_loc_dev,     fun_loc_host,     n_total_points), 56);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->T_loc_dev,       T_loc_host,       n_total_points), 57);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->r_loc_1d_dev,    r_loc_1d_host,    loc_K_host    ), 58);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->t_loc_1d_dev,    t_loc_1d_host,    loc_J_host    ), 59);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->xi_loc_dev,      xi_loc_host,      n_total_points), 60);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->eta_loc_dev,     eta_loc_host,     n_total_points), 61);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->zeta_loc_dev,    zeta_loc_host,    n_total_points), 62);
    print_CUDA_error_if_any(copy_host_to_device_bl(grid_on_dv->is_changed_dev,  is_changed_host,  n_total_points), 62);


    // neighboring mpi ranks
    int i=0;
    for (auto& neighbor_id : neighbors_id_host){
        grid_on_dv->neighbors_id_host[i] = neighbor_id;
        i++;
    }

    if (grid_on_dv->inverse_flag_host){
        print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->Tadj_loc_dev,    Tadj_loc_host,    n_total_points), 64);
    }

    // calculate total memory in Mb allocated on device
    size_t total_memory_in_bytes = 0;
    size_t free_memory_in_bytes = 0;
    cudaMemGetInfo(&free_memory_in_bytes, &total_memory_in_bytes);

}


void cuda_init_gpu_params(Grid_on_device* grid_on_dv, \
                          MPI_Comm inter_sub_comm, \
                          std::vector<std::vector<int>>& ijk_levels, \
                          CUSTOMREAL dr, CUSTOMREAL dt, CUSTOMREAL dp) {

    //printf("cuda_init_gpu_params\n");

    // prepare interprocess exchange information
    grid_on_dv->inter_sub_comm_host = inter_sub_comm;

    // memory allocation
    int n_levels = ijk_levels.size();

    int total_points = 0;
    int max_points = 0;
    for (auto &level : ijk_levels) {
        int n_points = level.size();
        total_points += n_points;
        if (n_points > max_points) {
            max_points = n_points;
        }
    }
    // store maximum number of points in a level
    grid_on_dv->n_max_points_level_host = max_points;
    // store total number of points in all levels
    grid_on_dv->n_total_points_host = total_points;

    grid_on_dv->n_points_per_level_host = new int[n_levels];
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**)&(grid_on_dv->i_id_all_level), total_points), 66);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**)&(grid_on_dv->j_id_all_level), total_points), 67);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**)&(grid_on_dv->k_id_all_level), total_points), 68);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**)&(grid_on_dv->n_points_per_level_dev), n_levels), 69);

    int* i_id_all_level_tmp = new int[total_points];
    int* j_id_all_level_tmp = new int[total_points];
    int* k_id_all_level_tmp = new int[total_points];

    // store the points in temporary arrays
    int current_pt = 0;
    int iip,jjt,kkr;
    for (int ilevel = 0; ilevel < n_levels; ilevel++){
        int n_points = ijk_levels[ilevel].size();
        for (int ipoint = 0; ipoint < n_points; ipoint++){
            V2I(ijk_levels[ilevel][ipoint], iip, jjt, kkr, grid_on_dv->loc_I_host, grid_on_dv->loc_J_host, grid_on_dv->loc_K_host);
            i_id_all_level_tmp[current_pt] = iip;
            j_id_all_level_tmp[current_pt] = jjt;
            k_id_all_level_tmp[current_pt] = kkr;
            current_pt++;
        }
    }

    print_CUDA_error_if_any(copy_host_to_device_i(grid_on_dv->i_id_all_level, i_id_all_level_tmp, total_points), 69);
    print_CUDA_error_if_any(copy_host_to_device_i(grid_on_dv->j_id_all_level, j_id_all_level_tmp, total_points), 70);
    print_CUDA_error_if_any(copy_host_to_device_i(grid_on_dv->k_id_all_level, k_id_all_level_tmp, total_points), 71);

    delete [] i_id_all_level_tmp;
    delete [] j_id_all_level_tmp;
    delete [] k_id_all_level_tmp;


    //int current_pt = 0;
    for (int ilevel = 0; ilevel < n_levels; ilevel++){
        int n_points = ijk_levels[ilevel].size();
        // store the number of points on each level
        grid_on_dv->n_points_per_level_host[ilevel] = n_points;
    }
    print_CUDA_error_if_any(copy_host_to_device_i(grid_on_dv->n_points_per_level_dev, grid_on_dv->n_points_per_level_host, n_levels), 70);

    grid_on_dv->dr_host = dr;
    grid_on_dv->dt_host = dt;
    grid_on_dv->dp_host = dp;

}


void cuda_reinitialize_grid(Grid_on_device* grid_on_dv, \
                            CUSTOMREAL* fac_a_loc_host, \
                            CUSTOMREAL* fac_b_loc_host, \
                            CUSTOMREAL* fac_c_loc_host, \
                            CUSTOMREAL* fac_f_loc_host, \
                            CUSTOMREAL* T0v_loc_host, \
                            CUSTOMREAL* T0r_loc_host, \
                            CUSTOMREAL* T0t_loc_host, \
                            CUSTOMREAL* T0p_loc_host, \
                            CUSTOMREAL* tau_loc_host, \
                            CUSTOMREAL* tau_old_loc_host, \
                            bool* is_changed_host){

    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->fac_a_loc_dev, fac_a_loc_host, grid_on_dv->n_total_points_host), 72);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->fac_b_loc_dev, fac_b_loc_host, grid_on_dv->n_total_points_host), 73);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->fac_c_loc_dev, fac_c_loc_host, grid_on_dv->n_total_points_host), 74);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->fac_f_loc_dev, fac_f_loc_host, grid_on_dv->n_total_points_host), 75);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->T0v_loc_dev, T0v_loc_host, grid_on_dv->n_total_points_host), 76);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->T0r_loc_dev, T0r_loc_host, grid_on_dv->n_total_points_host), 77);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->T0t_loc_dev, T0t_loc_host, grid_on_dv->n_total_points_host), 78);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->T0p_loc_dev, T0p_loc_host, grid_on_dv->n_total_points_host), 79);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->tau_loc_dev, tau_loc_host, grid_on_dv->n_total_points_host), 80);
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_on_dv->tau_old_loc_dev, tau_old_loc_host, grid_on_dv->n_total_points_host), 81);
    print_CUDA_error_if_any(copy_host_to_device_bl(grid_on_dv->is_changed_dev, is_changed_host, grid_on_dv->n_total_points_host), 82);

}



__global__ void cuda_prepare_boundary_data_to_send_jk_plane_negative( \
    int i_start_loc, int j_start_loc, int k_start_loc, \
    int i_end_loc, int j_end_loc, int k_end_loc, \
    int loc_I, int loc_J, int loc_K, \
    int n_ghost_layers, \
    int loc_I_excl_ghost, int loc_J_excl_ghost, int loc_K_excl_ghost, \
    CUSTOMREAL* bin_s, \
    CUSTOMREAL* tau_loc){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y + j_start_loc;
    int k = blockIdx.z * blockDim.z + threadIdx.z + k_start_loc;

    // return if ijk are out of bounds
    if (i > n_ghost_layers || j > j_end_loc || k > k_end_loc || i < 0 || j < j_start_loc || k < k_start_loc) {
        return;
    }

    int ijk = I2V_cuda(i+n_ghost_layers, j, k, loc_I, loc_J);
    bin_s[i*loc_J_excl_ghost*loc_K_excl_ghost + (k-k_start_loc)*loc_J_excl_ghost \
                                                    + (j-j_start_loc)] = tau_loc[ijk];

}


__global__ void cuda_prepare_boundary_data_to_send_jk_plane_positive(
    int i_start_loc, int j_start_loc, int k_start_loc, \
    int i_end_loc, int j_end_loc, int k_end_loc, \
    int loc_I, int loc_J, int loc_K, \
    int n_ghost_layers, \
    int loc_I_excl_ghost, int loc_J_excl_ghost, int loc_K_excl_ghost, \
    CUSTOMREAL* bip_s, \
    CUSTOMREAL* tau_loc){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y + j_start_loc;
    int k = blockIdx.z * blockDim.z + threadIdx.z + k_start_loc;

    // return if ijk are out of bounds
    if (i > n_ghost_layers || j > j_end_loc || k > k_end_loc || i < 0 || j < j_start_loc || k < k_start_loc) {
        return;
    }

    int ijk = I2V_cuda(loc_I-n_ghost_layers*2+i, j, k, loc_I, loc_J);
    bip_s[i*loc_J_excl_ghost*loc_K_excl_ghost + (k-k_start_loc)*loc_J_excl_ghost \
                                                        + (j-j_start_loc)] = tau_loc[ijk];

}

__global__ void cuda_prepare_boundary_data_to_send_ik_plane_negative( \
    int i_start_loc, int j_start_loc, int k_start_loc, \
    int i_end_loc, int j_end_loc, int k_end_loc, \
    int loc_I, int loc_J, int loc_K, \
    int n_ghost_layers, \
    int loc_I_excl_ghost, int loc_J_excl_ghost, int loc_K_excl_ghost, \
    CUSTOMREAL* bjn_s, \
    CUSTOMREAL* tau_loc){

    int i = blockIdx.x * blockDim.x + threadIdx.x + i_start_loc;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z + k_start_loc;

    // return if ijk are out of bounds
    if (i > i_end_loc || j > n_ghost_layers || k > k_end_loc || i < i_start_loc || j < 0 || k < k_start_loc) {
        return;
    }

    int ijk = I2V_cuda(i, j+n_ghost_layers, k, loc_I, loc_J);
    bjn_s[j*loc_I_excl_ghost*loc_K_excl_ghost \
                 + (k-k_start_loc)*loc_I_excl_ghost + (i-i_start_loc)] = tau_loc[ijk];

}


__global__ void cuda_prepare_boundary_data_to_send_ik_plane_positive( \
    int i_start_loc, int j_start_loc, int k_start_loc, \
    int i_end_loc, int j_end_loc, int k_end_loc, \
    int loc_I, int loc_J, int loc_K, \
    int n_ghost_layers, \
    int loc_I_excl_ghost, int loc_J_excl_ghost, int loc_K_excl_ghost, \
    CUSTOMREAL* bjp_s, \
    CUSTOMREAL* tau_loc){

    int i = blockIdx.x * blockDim.x + threadIdx.x + i_start_loc;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z + k_start_loc;

    // return if ijk are out of bounds
    if (i > i_end_loc || j > n_ghost_layers || k > k_end_loc || i < i_start_loc || j < 0 || k < k_start_loc) {
        return;
    }

    int ijk = I2V_cuda(i, loc_J-n_ghost_layers*2+j, k, loc_I, loc_J);
    bjp_s[j*loc_I_excl_ghost*loc_K_excl_ghost \
                 + (k-k_start_loc)*loc_I_excl_ghost + (i-i_start_loc)] =tau_loc[ijk];

}


__global__ void cuda_prepare_boundary_data_to_send_ij_plane_negative( \
    int i_start_loc, int j_start_loc, int k_start_loc, \
    int i_end_loc, int j_end_loc, int k_end_loc, \
    int loc_I, int loc_J, int loc_K, \
    int n_ghost_layers, \
    int loc_I_excl_ghost, int loc_J_excl_ghost, int loc_K_excl_ghost, \
    CUSTOMREAL* bkn_s, \
    CUSTOMREAL* tau_loc){

    int i = blockIdx.x * blockDim.x + threadIdx.x + i_start_loc;
    int j = blockIdx.y * blockDim.y + threadIdx.y + j_start_loc;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    // return if ijk are out of bounds
    if (i > i_end_loc || j > j_end_loc || k > n_ghost_layers || i < i_start_loc || j < j_start_loc || k < 0) {
        return;
    }

    int ijk = I2V_cuda(i, j, k+n_ghost_layers, loc_I, loc_J);
    bkn_s[k*loc_I_excl_ghost*loc_J_excl_ghost \
                 + (j-j_start_loc)*loc_I_excl_ghost + (i-i_start_loc)] = tau_loc[ijk];

}


__global__ void cuda_prepare_boundary_data_to_send_ij_plane_positive( \
    int i_start_loc, int j_start_loc, int k_start_loc, \
    int i_end_loc, int j_end_loc, int k_end_loc, \
    int loc_I, int loc_J, int loc_K, \
    int n_ghost_layers, \
    int loc_I_excl_ghost, int loc_J_excl_ghost, int loc_K_excl_ghost, \
    CUSTOMREAL* bkp_s, \
    CUSTOMREAL* tau_loc){

    int i = blockIdx.x * blockDim.x + threadIdx.x + i_start_loc;
    int j = blockIdx.y * blockDim.y + threadIdx.y + j_start_loc;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    // return if ijk are out of bounds
    if (i > i_end_loc || j > j_end_loc || k > n_ghost_layers || i < i_start_loc || j < j_start_loc || k < 0) {
        return;
    }

    int ijk = I2V_cuda(i, j, loc_K-n_ghost_layers*2+k, loc_I, loc_J);
    bkp_s[k*loc_I_excl_ghost*loc_J_excl_ghost \
                 + (j-j_start_loc)*loc_I_excl_ghost + (i-i_start_loc)] = tau_loc[ijk];

}

void cuda_prepare_boundary_data_to_send(Grid_on_device* grid){
    // define threads and blocks
    dim3 blocks, threads;

    // j-k plane negative
    if (grid->neighbors_id_host[0] != -1) {
        get_thread_block_for_ibound(grid->n_ghost_layers_host, grid->loc_J_excl_ghost_host, grid->loc_K_excl_ghost_host, &threads, &blocks);
        cuda_prepare_boundary_data_to_send_jk_plane_negative<<<blocks, threads>>>(\
            grid->i_start_loc_host, \
            grid->j_start_loc_host, \
            grid->k_start_loc_host, \
            grid->i_end_loc_host, \
            grid->j_end_loc_host, \
            grid->k_end_loc_host, \
            grid->loc_I_host, \
            grid->loc_J_host, \
            grid->loc_K_host, \
            grid->n_ghost_layers_host, \
            grid->loc_I_excl_ghost_host, \
            grid->loc_J_excl_ghost_host, \
            grid->loc_K_excl_ghost_host, \
            grid->bin_s_dev, \
            grid->tau_loc_dev);
    }

    // j-k plane positive
    if (grid->neighbors_id_host[1] != -1) {
        get_thread_block_for_ibound(grid->n_ghost_layers_host, grid->loc_J_excl_ghost_host, grid->loc_K_excl_ghost_host, &threads, &blocks);
        cuda_prepare_boundary_data_to_send_jk_plane_positive<<<blocks, threads>>>(\
            grid->i_start_loc_host, \
            grid->j_start_loc_host, \
            grid->k_start_loc_host, \
            grid->i_end_loc_host, \
            grid->j_end_loc_host, \
            grid->k_end_loc_host, \
            grid->loc_I_host, \
            grid->loc_J_host, \
            grid->loc_K_host, \
            grid->n_ghost_layers_host, \
            grid->loc_I_excl_ghost_host, \
            grid->loc_J_excl_ghost_host, \
            grid->loc_K_excl_ghost_host, \
            grid->bip_s_dev, \
            grid->tau_loc_dev);
    }

    // i-k plane negative
    if (grid->neighbors_id_host[2] != -1) {
        get_thread_block_for_jbound(grid->loc_I_excl_ghost_host, grid->n_ghost_layers_host, grid->loc_K_excl_ghost_host, &threads, &blocks);
        cuda_prepare_boundary_data_to_send_ik_plane_negative<<<blocks, threads>>>( \
            grid->i_start_loc_host, \
            grid->j_start_loc_host, \
            grid->k_start_loc_host, \
            grid->i_end_loc_host, \
            grid->j_end_loc_host, \
            grid->k_end_loc_host, \
            grid->loc_I_host, \
            grid->loc_J_host, \
            grid->loc_K_host, \
            grid->n_ghost_layers_host, \
            grid->loc_I_excl_ghost_host, \
            grid->loc_J_excl_ghost_host, \
            grid->loc_K_excl_ghost_host, \
            grid->bjn_s_dev, \
            grid->tau_loc_dev);
    }

    // i-k plane positive
    if (grid->neighbors_id_host[3] != -1) {
        get_thread_block_for_jbound(grid->loc_I_excl_ghost_host, grid->n_ghost_layers_host, grid->loc_K_excl_ghost_host, &threads, &blocks);
        cuda_prepare_boundary_data_to_send_ik_plane_positive<<<blocks, threads>>>( \
            grid->i_start_loc_host, \
            grid->j_start_loc_host, \
            grid->k_start_loc_host, \
            grid->i_end_loc_host, \
            grid->j_end_loc_host, \
            grid->k_end_loc_host, \
            grid->loc_I_host, \
            grid->loc_J_host, \
            grid->loc_K_host, \
            grid->n_ghost_layers_host, \
            grid->loc_I_excl_ghost_host, \
            grid->loc_J_excl_ghost_host, \
            grid->loc_K_excl_ghost_host, \
            grid->bjp_s_dev, \
            grid->tau_loc_dev);
    }

    // i-j plane negative
    if (grid->neighbors_id_host[4] != -1) {
        get_thread_block_for_kbound(grid->loc_I_excl_ghost_host, grid->loc_J_excl_ghost_host, grid->n_ghost_layers_host, &threads, &blocks);
        cuda_prepare_boundary_data_to_send_ij_plane_negative<<<blocks, threads>>>( \
            grid->i_start_loc_host, \
            grid->j_start_loc_host, \
            grid->k_start_loc_host, \
            grid->i_end_loc_host, \
            grid->j_end_loc_host, \
            grid->k_end_loc_host, \
            grid->loc_I_host, \
            grid->loc_J_host, \
            grid->loc_K_host, \
            grid->n_ghost_layers_host, \
            grid->loc_I_excl_ghost_host, \
            grid->loc_J_excl_ghost_host, \
            grid->loc_K_excl_ghost_host, \
            grid->bkn_s_dev, \
            grid->tau_loc_dev);
    }

    // i-j plane positive
    if (grid->neighbors_id_host[5] != -1) {
        get_thread_block_for_kbound(grid->loc_I_excl_ghost_host, grid->loc_J_excl_ghost_host, grid->n_ghost_layers_host, &threads, &blocks);
        cuda_prepare_boundary_data_to_send_ij_plane_positive<<<blocks, threads>>>( \
            grid->i_start_loc_host, \
            grid->j_start_loc_host, \
            grid->k_start_loc_host, \
            grid->i_end_loc_host, \
            grid->j_end_loc_host, \
            grid->k_end_loc_host, \
            grid->loc_I_host, \
            grid->loc_J_host, \
            grid->loc_K_host, \
            grid->n_ghost_layers_host, \
            grid->loc_I_excl_ghost_host, \
            grid->loc_J_excl_ghost_host, \
            grid->loc_K_excl_ghost_host, \
            grid->bkp_s_dev, \
            grid->tau_loc_dev);
    }

}


__global__ void cuda_assign_received_data_to_ghost_jk_plane_negative( \
    int i_start_loc, int j_start_loc, int k_start_loc, \
    int i_end_loc, int j_end_loc, int k_end_loc, \
    int loc_I, int loc_J, int loc_K, \
    int n_ghost_layers, \
    int loc_I_excl_ghost, int loc_J_excl_ghost, int loc_K_excl_ghost, \
    CUSTOMREAL* bin_r, \
    CUSTOMREAL* tau_loc){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y + j_start_loc;
    int k = blockIdx.z * blockDim.z + threadIdx.z + k_start_loc;

    // return if ijk are out of bounds
    if (i > n_ghost_layers || j > j_end_loc || k > k_end_loc || i < 0 || j < j_start_loc || k < k_start_loc) {
        return;
    }

    int ijk = I2V_cuda(i, j, k, loc_I, loc_J);
    tau_loc[ijk] = bin_r[i*loc_J_excl_ghost*loc_K_excl_ghost \
                                         + (k-k_start_loc)*loc_J_excl_ghost + (j-j_start_loc)];

}


__global__ void cuda_assign_received_data_to_ghost_jk_plane_positive( \
    int i_start_loc, int j_start_loc, int k_start_loc, \
    int i_end_loc, int j_end_loc, int k_end_loc, \
    int loc_I, int loc_J, int loc_K, \
    int n_ghost_layers, \
    int loc_I_excl_ghost, int loc_J_excl_ghost, int loc_K_excl_ghost, \
    CUSTOMREAL* bip_r, \
    CUSTOMREAL* tau_loc){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y + j_start_loc;
    int k = blockIdx.z * blockDim.z + threadIdx.z + k_start_loc;

    // return if ijk are out of bounds
    if (i > n_ghost_layers || j > j_end_loc || k > k_end_loc || i < 0 || j < j_start_loc || k < k_start_loc) {
    return;
    }

    int ijk = I2V_cuda(loc_I-n_ghost_layers+i, j, k, loc_I, loc_J);
    tau_loc[ijk] = bip_r[i*loc_J_excl_ghost*loc_K_excl_ghost \
                                          + (k-k_start_loc)*loc_J_excl_ghost + (j-j_start_loc)];

}


__global__ void cuda_assign_received_data_to_ghost_ik_plane_negative( \
    int i_start_loc, int j_start_loc, int k_start_loc, \
    int i_end_loc, int j_end_loc, int k_end_loc, \
    int loc_I, int loc_J, int loc_K, \
    int n_ghost_layers, \
    int loc_I_excl_ghost, int loc_J_excl_ghost, int loc_K_excl_ghost, \
    CUSTOMREAL* bjn_r, \
    CUSTOMREAL* tau_loc){

    int i = blockIdx.x * blockDim.x + threadIdx.x + i_start_loc;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z + k_start_loc;

    // return if ijk are out of bounds
    if (i > i_end_loc || j > n_ghost_layers || k > k_end_loc || i < i_start_loc || j < 0 || k < k_start_loc) {
        return;
    }

    int ijk = I2V_cuda(i, j, k, loc_I, loc_J);
    tau_loc[ijk] = bjn_r[j*loc_I_excl_ghost*loc_K_excl_ghost + \
                                            (k-k_start_loc)*loc_I_excl_ghost + (i-i_start_loc)];

}


__global__ void cuda_assign_received_data_to_ghost_ik_plane_positive( \
    int i_start_loc, int j_start_loc, int k_start_loc, \
    int i_end_loc, int j_end_loc, int k_end_loc, \
    int loc_I, int loc_J, int loc_K, \
    int n_ghost_layers, \
    int loc_I_excl_ghost, int loc_J_excl_ghost, int loc_K_excl_ghost, \
    CUSTOMREAL* bjp_r, \
    CUSTOMREAL* tau_loc){

    int i = blockIdx.x * blockDim.x + threadIdx.x + i_start_loc;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z + k_start_loc;

    // return if ijk are out of bounds
    if (i > i_end_loc || j > n_ghost_layers || k > k_end_loc || i < i_start_loc || j < 0 || k < k_start_loc) {
        return;
    }

    int ijk = I2V_cuda(i, loc_J-n_ghost_layers+j, k, loc_I, loc_J);
    tau_loc[ijk] = bjp_r[j*loc_I_excl_ghost*loc_K_excl_ghost + \
                                            (k-k_start_loc)*loc_I_excl_ghost + (i-i_start_loc)];

}


__global__ void cuda_assign_received_data_to_ghost_ij_plane_negative( \
    int i_start_loc, int j_start_loc, int k_start_loc, \
    int i_end_loc, int j_end_loc, int k_end_loc, \
    int loc_I, int loc_J, int loc_K, \
    int n_ghost_layers, \
    int loc_I_excl_ghost, int loc_J_excl_ghost, int loc_K_excl_ghost, \
    CUSTOMREAL* bkn_r, \
    CUSTOMREAL* tau_loc){

    int i = blockIdx.x * blockDim.x + threadIdx.x + i_start_loc;
    int j = blockIdx.y * blockDim.y + threadIdx.y + j_start_loc;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    // return if ijk are out of bounds
    if (i > i_end_loc || j > j_end_loc || k > n_ghost_layers || i < i_start_loc || j < j_start_loc || k < 0) {
        return;
    }

    int ijk = I2V_cuda(i, j, k, loc_I, loc_J);
    tau_loc[ijk] = bkn_r[k*loc_I_excl_ghost*loc_J_excl_ghost + \
                                            (j-j_start_loc)*loc_I_excl_ghost + (i-i_start_loc)];

}


__global__ void cuda_assign_received_data_to_ghost_ij_plane_positive( \
    int i_start_loc, int j_start_loc, int k_start_loc, \
    int i_end_loc, int j_end_loc, int k_end_loc, \
    int loc_I, int loc_J, int loc_K, \
    int n_ghost_layers, \
    int loc_I_excl_ghost, int loc_J_excl_ghost, int loc_K_excl_ghost, \
    CUSTOMREAL* bkp_r, \
    CUSTOMREAL* tau_loc){

    int i = blockIdx.x * blockDim.x + threadIdx.x + i_start_loc;
    int j = blockIdx.y * blockDim.y + threadIdx.y + j_start_loc;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    // return if ijk are out of bounds
    if (i > i_end_loc || j > j_end_loc || k > n_ghost_layers || i < i_start_loc || j < j_start_loc || k < 0) {
        return;
    }

    int ijk = I2V_cuda(i, j, loc_K-n_ghost_layers+k, loc_I, loc_J);
    tau_loc[ijk] = bkp_r[k*loc_I_excl_ghost*loc_J_excl_ghost + \
                                            (j-j_start_loc)*loc_I_excl_ghost + (i-i_start_loc)];
}


void cuda_assign_received_data_to_ghost(Grid_on_device* grid){

    // define threads and blocks
    dim3 blocks, threads;

    // j-k plane negative
    if (grid->neighbors_id_host[0] != -1) {
        get_thread_block_for_ibound(grid->n_ghost_layers_host, grid->loc_J_excl_ghost_host, grid->loc_K_excl_ghost_host, &threads, &blocks);
        cuda_assign_received_data_to_ghost_jk_plane_negative<<<blocks, threads>>>(\
        grid->i_start_loc_host, \
        grid->j_start_loc_host, \
        grid->k_start_loc_host, \
        grid->i_end_loc_host, \
        grid->j_end_loc_host, \
        grid->k_end_loc_host, \
        grid->loc_I_host, \
        grid->loc_J_host, \
        grid->loc_K_host, \
        grid->n_ghost_layers_host, \
        grid->loc_I_excl_ghost_host, \
        grid->loc_J_excl_ghost_host, \
        grid->loc_K_excl_ghost_host, \
        grid->bin_r_dev, \
        grid->tau_loc_dev);
    }

    // j-k plane positive
    if (grid->neighbors_id_host[1] != -1) {
        get_thread_block_for_ibound(grid->n_ghost_layers_host, grid->loc_J_excl_ghost_host, grid->loc_K_excl_ghost_host, &threads, &blocks);
        cuda_assign_received_data_to_ghost_jk_plane_positive<<<blocks, threads>>>(\
        grid->i_start_loc_host, \
        grid->j_start_loc_host, \
        grid->k_start_loc_host, \
        grid->i_end_loc_host, \
        grid->j_end_loc_host, \
        grid->k_end_loc_host, \
        grid->loc_I_host, \
        grid->loc_J_host, \
        grid->loc_K_host, \
        grid->n_ghost_layers_host, \
        grid->loc_I_excl_ghost_host, \
        grid->loc_J_excl_ghost_host, \
        grid->loc_K_excl_ghost_host, \
        grid->bip_r_dev, \
        grid->tau_loc_dev);
    }

    // i-k plane negative
    if (grid->neighbors_id_host[2] != -1) {
        get_thread_block_for_jbound(grid->loc_I_excl_ghost_host, grid->n_ghost_layers_host, grid->loc_K_excl_ghost_host, &threads, &blocks);
        cuda_assign_received_data_to_ghost_ik_plane_negative<<<blocks, threads>>>( \
        grid->i_start_loc_host, \
        grid->j_start_loc_host, \
        grid->k_start_loc_host, \
        grid->i_end_loc_host, \
        grid->j_end_loc_host, \
        grid->k_end_loc_host, \
        grid->loc_I_host, \
        grid->loc_J_host, \
        grid->loc_K_host, \
        grid->n_ghost_layers_host, \
        grid->loc_I_excl_ghost_host, \
        grid->loc_J_excl_ghost_host, \
        grid->loc_K_excl_ghost_host, \
        grid->bjn_r_dev, \
        grid->tau_loc_dev);
    }

    // i-k plane positive
    if (grid->neighbors_id_host[3] != -1) {
        get_thread_block_for_jbound(grid->loc_I_excl_ghost_host, grid->n_ghost_layers_host, grid->loc_K_excl_ghost_host, &threads, &blocks);
        cuda_assign_received_data_to_ghost_ik_plane_positive<<<blocks, threads>>>( \
        grid->i_start_loc_host, \
        grid->j_start_loc_host, \
        grid->k_start_loc_host, \
        grid->i_end_loc_host, \
        grid->j_end_loc_host, \
        grid->k_end_loc_host, \
        grid->loc_I_host, \
        grid->loc_J_host, \
        grid->loc_K_host, \
        grid->n_ghost_layers_host, \
        grid->loc_I_excl_ghost_host, \
        grid->loc_J_excl_ghost_host, \
        grid->loc_K_excl_ghost_host, \
        grid->bjp_r_dev, \
        grid->tau_loc_dev);
    }

    // i-j plane negative
    if (grid->neighbors_id_host[4] != -1) {
        get_thread_block_for_kbound(grid->loc_I_excl_ghost_host, grid->loc_J_excl_ghost_host, grid->n_ghost_layers_host, &threads, &blocks);
        cuda_assign_received_data_to_ghost_ij_plane_negative<<<blocks, threads>>>( \
        grid->i_start_loc_host, \
        grid->j_start_loc_host, \
        grid->k_start_loc_host, \
        grid->i_end_loc_host, \
        grid->j_end_loc_host, \
        grid->k_end_loc_host, \
        grid->loc_I_host, \
        grid->loc_J_host, \
        grid->loc_K_host, \
        grid->n_ghost_layers_host, \
        grid->loc_I_excl_ghost_host, \
        grid->loc_J_excl_ghost_host, \
        grid->loc_K_excl_ghost_host, \
        grid->bkn_r_dev, \
        grid->tau_loc_dev);
    }

    // i-j plane positive
    if (grid->neighbors_id_host[5] != -1) {
        get_thread_block_for_kbound(grid->loc_I_excl_ghost_host, grid->loc_J_excl_ghost_host, grid->n_ghost_layers_host, &threads, &blocks);
        cuda_assign_received_data_to_ghost_ij_plane_positive<<<blocks, threads>>>(\
        grid->i_start_loc_host, \
        grid->j_start_loc_host, \
        grid->k_start_loc_host, \
        grid->i_end_loc_host, \
        grid->j_end_loc_host, \
        grid->k_end_loc_host, \
        grid->loc_I_host, \
        grid->loc_J_host, \
        grid->loc_K_host, \
        grid->n_ghost_layers_host, \
        grid->loc_I_excl_ghost_host, \
        grid->loc_J_excl_ghost_host, \
        grid->loc_K_excl_ghost_host, \
        grid->bkp_r_dev, \
        grid->tau_loc_dev);
    }

}


void cuda_send_recev_boundary_data(Grid_on_device* grid){

    if (grid->n_subdomains_host <= 1) return;

    cuda_synchronize_all_sub(grid->inter_sub_comm_host);

    cuda_prepare_boundary_data_to_send(grid);

    MPI_Request* mpi_reqs = new MPI_Request[6];
    // initialize mpi_reqs
    for (int i = 0; i < 6; i++) {
        mpi_reqs[i] = MPI_REQUEST_NULL;
    }

    // i-direction negative
    if (grid->neighbors_id_host[0] != -1) {
//        // send boundary layer to neighbor
//        cuda_isend_cr(grid->bin_s_dev, grid->n_grid_bound_i_host*grid->n_ghost_layers_host, grid->neighbors_id_host[0], grid->inter_sub_comm_host, mpi_reqs[0]);
//        // receive boundary layer from neighbor
//        cuda_irecv_cr(grid->bin_r_dev, grid->n_grid_bound_i_host*grid->n_ghost_layers_host, grid->neighbors_id_host[0], grid->inter_sub_comm_host, mpi_reqs[0]);

        // send boundary layer to neighbor
        cuda_send_cr(grid->bin_s_dev, grid->n_grid_bound_i_host*grid->n_ghost_layers_host, grid->neighbors_id_host[0], grid->inter_sub_comm_host);
        // receive boundary layer from neighbor
        cuda_recv_cr(grid->bin_r_dev, grid->n_grid_bound_i_host*grid->n_ghost_layers_host, grid->neighbors_id_host[0], grid->inter_sub_comm_host);

    }
    // i-direction positive
    if (grid->neighbors_id_host[1] != -1) {
//        // send boundary layer to neighbor
//        cuda_isend_cr(grid->bip_s_dev, grid->n_grid_bound_i_host*grid->n_ghost_layers_host, grid->neighbors_id_host[1], grid->inter_sub_comm_host, mpi_reqs[1]);
//        // receive boundary layer from neighbor
//        cuda_irecv_cr(grid->bip_r_dev, grid->n_grid_bound_i_host*grid->n_ghost_layers_host, grid->neighbors_id_host[1], grid->inter_sub_comm_host, mpi_reqs[1]);
        // send boundary layer to neighbor
        cuda_send_cr(grid->bip_s_dev, grid->n_grid_bound_i_host*grid->n_ghost_layers_host, grid->neighbors_id_host[1], grid->inter_sub_comm_host);
        // receive boundary layer from neighbor
        cuda_recv_cr(grid->bip_r_dev, grid->n_grid_bound_i_host*grid->n_ghost_layers_host, grid->neighbors_id_host[1], grid->inter_sub_comm_host);
    }
    // j-direction negative
    if (grid->neighbors_id_host[2] != -1) {
//        // send boundary layer to neighbor
//        cuda_isend_cr(grid->bjn_s_dev, grid->n_grid_bound_j_host*grid->n_ghost_layers_host, grid->neighbors_id_host[2], grid->inter_sub_comm_host, mpi_reqs[2]);
//        // receive boundary layer from neighbor
//        cuda_irecv_cr(grid->bjn_r_dev, grid->n_grid_bound_j_host*grid->n_ghost_layers_host, grid->neighbors_id_host[2], grid->inter_sub_comm_host, mpi_reqs[2]);
        // send boundary layer to neighbor
        cuda_send_cr(grid->bjn_s_dev, grid->n_grid_bound_j_host*grid->n_ghost_layers_host, grid->neighbors_id_host[2], grid->inter_sub_comm_host);
        // receive boundary layer from neighbor
        cuda_recv_cr(grid->bjn_r_dev, grid->n_grid_bound_j_host*grid->n_ghost_layers_host, grid->neighbors_id_host[2], grid->inter_sub_comm_host);
    }
    // j-direction positive
    if (grid->neighbors_id_host[3] != -1) {
//        // send boundary layer to neighbor
//        cuda_isend_cr(grid->bjp_s_dev, grid->n_grid_bound_j_host*grid->n_ghost_layers_host, grid->neighbors_id_host[3], grid->inter_sub_comm_host, mpi_reqs[3]);
//        // receive boundary layer from neighbor
//        cuda_irecv_cr(grid->bjp_r_dev, grid->n_grid_bound_j_host*grid->n_ghost_layers_host, grid->neighbors_id_host[3], grid->inter_sub_comm_host, mpi_reqs[3]);
        // send boundary layer to neighbor
        cuda_send_cr(grid->bjp_s_dev, grid->n_grid_bound_j_host*grid->n_ghost_layers_host, grid->neighbors_id_host[3], grid->inter_sub_comm_host);
        // receive boundary layer from neighbor
        cuda_recv_cr(grid->bjp_r_dev, grid->n_grid_bound_j_host*grid->n_ghost_layers_host, grid->neighbors_id_host[3], grid->inter_sub_comm_host);
    }
    // k-direction negative
    if (grid->neighbors_id_host[4] != -1) {
//        // send boundary layer to neighbor
//        cuda_isend_cr(grid->bkn_s_dev, grid->n_grid_bound_k_host*grid->n_ghost_layers_host, grid->neighbors_id_host[4], grid->inter_sub_comm_host, mpi_reqs[4]);
//        // receive boundary layer from neighbor
//        cuda_irecv_cr(grid->bkn_r_dev, grid->n_grid_bound_k_host*grid->n_ghost_layers_host, grid->neighbors_id_host[4], grid->inter_sub_comm_host, mpi_reqs[4]);
        // send boundary layer to neighbor
        cuda_send_cr(grid->bkn_s_dev, grid->n_grid_bound_k_host*grid->n_ghost_layers_host, grid->neighbors_id_host[4], grid->inter_sub_comm_host);
        // receive boundary layer from neighbor
        cuda_recv_cr(grid->bkn_r_dev, grid->n_grid_bound_k_host*grid->n_ghost_layers_host, grid->neighbors_id_host[4], grid->inter_sub_comm_host);
    }
    // k-direction positive
    if (grid->neighbors_id_host[5] != -1) {
//        // send boundary layer to neighbor
//        cuda_isend_cr(grid->bkp_s_dev, grid->n_grid_bound_k_host*grid->n_ghost_layers_host, grid->neighbors_id_host[5], grid->inter_sub_comm_host, mpi_reqs[5]);
//        // receive boundary layer from neighbor
//        cuda_irecv_cr(grid->bkp_r_dev, grid->n_grid_bound_k_host*grid->n_ghost_layers_host, grid->neighbors_id_host[5], grid->inter_sub_comm_host, mpi_reqs[5]);
        // send boundary layer to neighbor
        cuda_send_cr(grid->bkp_s_dev, grid->n_grid_bound_k_host*grid->n_ghost_layers_host, grid->neighbors_id_host[5], grid->inter_sub_comm_host);
        // receive boundary layer from neighbor
        cuda_recv_cr(grid->bkp_r_dev, grid->n_grid_bound_k_host*grid->n_ghost_layers_host, grid->neighbors_id_host[5], grid->inter_sub_comm_host);
    }

    // wait for finishing communication
//    for (int i = 0; i < 6; i++) {
//        if (grid->neighbors_id_host[i] != -1) {
//            cuda_wait_req(mpi_reqs[i]);
//        }
//    }

    delete [] mpi_reqs;

    cuda_assign_received_data_to_ghost(grid);

}

