#include "grid.h"

Grid::Grid(InputParams& IP, IO_utils& io) {
    stdout_by_main("--- grid object initialization starts. ---");

    // initialize grid parameters are done by only the main process of each subdomain
    if (subdom_main) {
        // domain decomposition
        init_decomposition(IP);

        // memory allocation
        memory_allocation();
    }

    // allocate memory for shm arrays if necessary
    if (n_subprocs > 1)
        shm_memory_allocation();

    synchronize_all_world();

    // setup grid parameters
    if (subdom_main)
        setup_grid_params(IP, io);

}

Grid::~Grid() {
    // deallocate memory from all arrays
    if (subdom_main) memory_deallocation();

    if (n_subprocs > 1) {
        // deallocate memory from shm arrays
        shm_memory_deallocation();
    }
}

void Grid::init_decomposition(InputParams& IP) {

    int tmp_rank = 0;

    for (int kk = 0; kk < ndiv_k; kk++)
        for (int jj = 0; jj < ndiv_j; jj++)
            for (int ii = 0; ii < ndiv_i; ii++) {
                tmp_rank = kk * ndiv_i * ndiv_j + jj * ndiv_i + ii;
                if (tmp_rank == myrank) {
                    // check neighbors on i axis
                    if (ndiv_i > 1) {
                        if (ii == 0) {                 // neighbor on the bound +i
                            neighbors_id[1] = tmp_rank + 1;
                        } else if (ii == ndiv_i - 1) { // neighbor on the bound -i
                            neighbors_id[0] = tmp_rank - 1;
                        } else {                       // middle layers
                            neighbors_id[0] = tmp_rank - 1;
                            neighbors_id[1] = tmp_rank + 1;
                        }
                    }

                    // check neighbors on j axis
                    if (ndiv_j > 1) {
                        if (jj == 0) {                 // neighbor on the bound +j
                            neighbors_id[3] = tmp_rank + ndiv_i;
                        } else if (jj == ndiv_j - 1) { // neighbor on the bound -j
                            neighbors_id[2] = tmp_rank - ndiv_i;
                        } else {                       // middle layers
                            neighbors_id[2] = tmp_rank - ndiv_i;
                            neighbors_id[3] = tmp_rank + ndiv_i;
                        }
                    }

                    // check neighbors on k axis
                    if (ndiv_k > 1) {
                        if (kk == 0) {                 // neighbor on the bound +k
                            neighbors_id[5] = tmp_rank + ndiv_i * ndiv_j;
                        } else if (kk == ndiv_k - 1) { // neighbor on the bound -k
                            neighbors_id[4] = tmp_rank - ndiv_i * ndiv_j;
                        } else {                       // middle layers
                            neighbors_id[4] = tmp_rank - ndiv_i * ndiv_j;
                            neighbors_id[5] = tmp_rank + ndiv_i * ndiv_j;
                        }
                    }

                    // store the layer id
                    domain_i = ii;
                    domain_j = jj;
                    domain_k = kk;

                    // stop searching immediately if this rank found the local ids
                    goto end_pos_detected;
                }
            }

    end_pos_detected:

    // calculate the number of grid points for this subdomain
    loc_I = (int)ngrid_i/ndiv_i;
    loc_J = (int)ngrid_j/ndiv_j;
    loc_K = (int)ngrid_k/ndiv_k;

    // calculate the number of grid points for this subdomain
    offset_i = domain_i * loc_I;
    offset_j = domain_j * loc_J;
    offset_k = domain_k * loc_K;

    // add modulus to the last rank
    if (i_last()) loc_I += ngrid_i % ndiv_i;
    if (j_last()) loc_J += ngrid_j % ndiv_j;
    if (k_last()) loc_K += ngrid_k % ndiv_k;

    // number of grid points on the boundary surfaces
    n_grid_bound_i = loc_J * loc_K;
    n_grid_bound_j = loc_I * loc_K;
    n_grid_bound_k = loc_I * loc_J;

    // store the number of grids excluding the ghost grids
    loc_I_excl_ghost = loc_I;
    loc_J_excl_ghost = loc_J;
    loc_K_excl_ghost = loc_K;

    // initialize iteration start/end
    i_start_loc = 0;
    j_start_loc = 0;
    k_start_loc = 0;
    i_end_loc   = loc_I-1; // end id needs to be included for the iterations
    j_end_loc   = loc_J-1; // thus please pay attaintion to the for loop
    k_end_loc   = loc_K-1; // to be e.g. for (i = i_start_loc; i <= i_end_loc; i++)

    // add ghost nodes to comunicate with neighbor domains
    // i axis
    if (ndiv_i > 1) {
        if (i_first() || i_last()) {
            loc_I += n_ghost_layers;
        } else {
            loc_I += n_ghost_layers * 2;
        }
        // store the start/end ids of iteration for axis i
        if (!i_first()) i_start_loc = n_ghost_layers;
        if (!i_last())  i_end_loc   = loc_I - n_ghost_layers - 1;
        if (i_last())   i_end_loc   = loc_I - 1;
    }

    // j axis
    if (ndiv_j > 1) {
        if (j_first() || j_last()) {
            loc_J += n_ghost_layers;
        } else {
            loc_J += n_ghost_layers * 2;
        }
        // store the start/end ids of iteration for axis j
        if (!j_first()) j_start_loc = n_ghost_layers;
        if (!j_last())  j_end_loc   = loc_J - n_ghost_layers - 1;
        if (j_last())   j_end_loc   = loc_J - 1;
    }
    // k axis
    if (ndiv_k > 1) {
        if (k_first() || k_last()) {
            loc_K += n_ghost_layers;
        } else {
            loc_K += n_ghost_layers * 2;
        }
        // store the start/end ids of iteration for axis k
        if (!k_first()) k_start_loc = n_ghost_layers;
        if (!k_last())  k_end_loc   = loc_K - n_ghost_layers - 1;
        if (k_last())   k_end_loc   = loc_K - 1;
    }

    // calculate the nodes' and elms' offset for data output
    loc_nnodes =  loc_I_excl_ghost      *  loc_J_excl_ghost      *  loc_K_excl_ghost;
    loc_nelms  = (loc_I_excl_ghost - 1) * (loc_J_excl_ghost - 1) * (loc_K_excl_ghost - 1);

    nnodes = new int[nprocs];
    nelms  = new int[nprocs];

    // gather the number of nodes and elements
    allgather_i_single(&loc_nnodes, nnodes);
    allgather_i_single(&loc_nelms, nelms);

    if (myrank > 0){
        // sum nnodes and nelms
        for (int i = 0; i < myrank; i++) {
            offset_nnodes += nnodes[i];
            offset_nelms  += nelms[i];
        }
    }

    // setup number of elements and nodes for visualization (additinoal 1 layer for connection)
    loc_I_vis = loc_I_excl_ghost;
    loc_J_vis = loc_J_excl_ghost;
    loc_K_vis = loc_K_excl_ghost;
    if (!i_last()) loc_I_vis++;
    if (!j_last()) loc_J_vis++;
    if (!k_last()) loc_K_vis++;

    loc_nnodes_vis = loc_I_vis * loc_J_vis * loc_K_vis;
    loc_nelms_vis  = (loc_I_vis - 1) * (loc_J_vis - 1) * (loc_K_vis - 1);
    nnodes_vis = new int[nprocs];
    nelms_vis  = new int[nprocs];
    allgather_i_single(&loc_nnodes_vis, nnodes_vis);
    allgather_i_single(&loc_nelms_vis, nelms_vis);
    if (myrank > 0){
        // sum nnodes and nelms
        for (int i = 0; i < myrank; i++) {
            offset_nnodes_vis += nnodes_vis[i];
            offset_nelms_vis  += nelms_vis[i];
        }
    }
    i_start_vis = i_start_loc;
    j_start_vis = j_start_loc;
    k_start_vis = k_start_loc;
    i_end_vis = i_end_loc;
    j_end_vis = j_end_loc;
    k_end_vis = k_end_loc;
    // add one layer for connection to the next subdomain
    if (!i_last()) i_end_vis++;
    if (!j_last()) j_end_vis++;
    if (!k_last()) k_end_vis++;

    // inversion setup
    // check if inversion grids are needed
    if (IP.get_run_mode()==1){
        inverse_flag = true;
        setup_inversion_grids(IP);
    } else {
        inverse_flag = false;
    }

    // check subdomains contacting lines and points
    // ij nn
    if (neighbors_id[0] != -1 && neighbors_id[2] != -1)
        neighbors_id_ij[0] = neighbors_id[2]-1;
    // ij np
    if (neighbors_id[0] != -1 && neighbors_id[3] != -1)
        neighbors_id_ij[1] = neighbors_id[3]-1;
    // ij pn
    if (neighbors_id[1] != -1 && neighbors_id[2] != -1)
        neighbors_id_ij[2] = neighbors_id[2]+1;
    // ij pp
    if (neighbors_id[1] != -1 && neighbors_id[3] != -1)
        neighbors_id_ij[3] = neighbors_id[3]+1;
    // jk nn
    if (neighbors_id[2] != -1 && neighbors_id[4] != -1)
        neighbors_id_jk[0] = neighbors_id[4]-ndiv_i;
    // jk np
    if (neighbors_id[2] != -1 && neighbors_id[5] != -1)
        neighbors_id_jk[1] = neighbors_id[5]-ndiv_i;
    // jk pn
    if (neighbors_id[3] != -1 && neighbors_id[4] != -1)
        neighbors_id_jk[2] = neighbors_id[4]+ndiv_i;
    // jk pp
    if (neighbors_id[3] != -1 && neighbors_id[5] != -1)
        neighbors_id_jk[3] = neighbors_id[5]+ndiv_i;
    // ik nn
    if (neighbors_id[0] != -1 && neighbors_id[4] != -1)
        neighbors_id_ik[0] = neighbors_id[4]-1;
    // ik np
    if (neighbors_id[0] != -1 && neighbors_id[5] != -1)
        neighbors_id_ik[1] = neighbors_id[5]-1;
    // ik pn
    if (neighbors_id[1] != -1 && neighbors_id[4] != -1)
        neighbors_id_ik[2] = neighbors_id[4]+1;
    // ik pp
    if (neighbors_id[1] != -1 && neighbors_id[5] != -1)
        neighbors_id_ik[3] = neighbors_id[5]+1;
    // ijk nnn
    if (neighbors_id[0] != -1 && neighbors_id[2] != -1 && neighbors_id[4] != -1)
        neighbors_id_ijk[0] = neighbors_id[4]-ndiv_i;
    // ijk nnp
    if (neighbors_id[0] != -1 && neighbors_id[2] != -1 && neighbors_id[5] != -1)
        neighbors_id_ijk[1] = neighbors_id[5]-ndiv_i;
    // ijk npn
    if (neighbors_id[0] != -1 && neighbors_id[3] != -1 && neighbors_id[4] != -1)
        neighbors_id_ijk[2] = neighbors_id[4]+ndiv_i;
    // ijk npp
    if (neighbors_id[0] != -1 && neighbors_id[3] != -1 && neighbors_id[5] != -1)
        neighbors_id_ijk[3] = neighbors_id[5]+ndiv_i;
    // ijk pnn
    if (neighbors_id[1] != -1 && neighbors_id[2] != -1 && neighbors_id[4] != -1)
        neighbors_id_ijk[4] = neighbors_id[4]-ndiv_i;
    // ijk pnp
    if (neighbors_id[1] != -1 && neighbors_id[2] != -1 && neighbors_id[5] != -1)
        neighbors_id_ijk[5] = neighbors_id[5]-ndiv_i;
    // ijk ppn
    if (neighbors_id[1] != -1 && neighbors_id[3] != -1 && neighbors_id[4] != -1)
        neighbors_id_ijk[6] = neighbors_id[4]+ndiv_i;
    // ijk ppp
    if (neighbors_id[1] != -1 && neighbors_id[3] != -1 && neighbors_id[5] != -1)
        neighbors_id_ijk[7] = neighbors_id[5]+ndiv_i;

    // debug output
    stdout_by_main("domain decomposition initialization end.");
    for (int i = 0; i < nprocs; i++) {
        synchronize_all_inter();

        if (i == myrank) {
            std::cout << "--- myrank: ---" << myrank << std::endl;

            if (i == 0) {
                std::cout << "ndiv_i j k : "  << ndiv_i  << " " << ndiv_j  << " " << ndiv_k  << std::endl;
                std::cout << "ngrid_i j k : " << ngrid_i << " " << ngrid_j << " " << ngrid_k << std::endl;
            }

            std::cout << "domain_i j k: "          << domain_i << " " << domain_j << " " << domain_k << std::endl;
            std::cout << "loc_I J K: "             << loc_I << " " << loc_J << " " << loc_K << std::endl;
            std::cout << "loc_I_excl_ghost J_excl_ghost K_excl_ghost: " << loc_I_excl_ghost << " " << loc_J_excl_ghost << " " << loc_K_excl_ghost << std::endl;
            std::cout << "i_start_loc i_end_loc: " << i_start_loc << " " << i_end_loc << std::endl;
            std::cout << "j_start_loc j_end_loc: " << j_start_loc << " " << j_end_loc << std::endl;
            std::cout << "k_start_loc k_end_loc: " << k_start_loc << " " << k_end_loc << std::endl;
            std::cout << "offset_nnodes offset_nelms: " << offset_nnodes << " " << offset_nelms << std::endl;
            std::cout << "n total local grids: "   << loc_I*loc_J*loc_K << "   max vector elm id:  " << I2V(loc_I-1,loc_J-1,loc_K-1) << std::endl;
            // print neighbors_id
            std::cout << "neighbors_id: ";
            for (int i = 0; i < 6; i++) {
                std::cout << neighbors_id[i] << " ";
            }
            std::cout << std::endl;
        }
    }

    synchronize_all_inter();
}


void Grid::setup_inversion_grids(InputParams& IP) {

    n_inv_grids = IP.get_n_inversion_grid();

    if(IP.get_type_dep_inv() == 0){
        ngrid_k_inv = IP.get_n_inv_r();
    } else if (IP.get_type_dep_inv() == 1) {
        ngrid_k_inv = IP.get_n_inv_r_flex();
    } else {
        std::cout << "unknown type of inversion grid" << std::endl;
        exit(1);
    }

    if(IP.get_type_lat_inv() == 0){
        ngrid_j_inv = IP.get_n_inv_t();
    } else if (IP.get_type_lat_inv() == 1) {
        ngrid_j_inv = IP.get_n_inv_t_flex();
    } else {
        std::cout << "unknown type of inversion grid" << std::endl;
        exit(1);
    }

    if(IP.get_type_lon_inv() == 0){
        ngrid_i_inv = IP.get_n_inv_p();
    } else if (IP.get_type_lon_inv() == 1) {
        ngrid_i_inv = IP.get_n_inv_p_flex();
    } else {
        std::cout << "unknown type of inversion grid" << std::endl;
        exit(1);
    }

    // +1 grid is added here for placing the inversion grid homogeneousily at the end of the domain
    // TODO : check if new inv grid type accept this +1
    n_inv_I_loc = ngrid_i_inv+1;
    n_inv_J_loc = ngrid_j_inv+1;
    n_inv_K_loc = ngrid_k_inv+1;

}


// allocate memory for arrays, called only for subdom_main
void Grid::memory_allocation() {
    // 3d arrays
    int n_total_loc_grid_points = loc_I * loc_J * loc_K;
    // skip allocation of shm arrays if not using shm
    if (sub_nprocs <= 1) {
#ifdef USE_CUDA
        if(use_gpu)
            cudaMallocHost((void**)&tau_loc, n_total_loc_grid_points * sizeof(CUSTOMREAL));
        else
            tau_loc     = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
#else
            tau_loc     = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
#endif

        T_loc       = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
        tau_old_loc = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);

        T0r_loc     = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
        T0t_loc     = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
        T0p_loc     = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
        T0v_loc     = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
        fac_a_loc   = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
        fac_b_loc   = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
        fac_c_loc   = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
        fac_f_loc   = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
        fun_loc     = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
        is_changed  = (bool *)       malloc(sizeof(bool)       * n_total_loc_grid_points);

        // 1d arrays
        p_loc_1d = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * loc_I);
        t_loc_1d = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * loc_J);
        r_loc_1d = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * loc_K);
    }

    if (if_test)
        u_loc   = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);

    // ghost layer arrays
    if (neighbors_id[0] != -1) {
        bin_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *n_ghost_layers*n_grid_bound_i);
        bin_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *n_ghost_layers*n_grid_bound_i);
    }
    if (neighbors_id[1] != -1) {
        bip_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *n_ghost_layers*n_grid_bound_i);
        bip_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *n_ghost_layers*n_grid_bound_i);
    }
    if (neighbors_id[2] != -1) {
        bjn_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *n_ghost_layers*n_grid_bound_j);
        bjn_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *n_ghost_layers*n_grid_bound_j);
    }
    if (neighbors_id[3] != -1) {
        bjp_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *n_ghost_layers*n_grid_bound_j);
        bjp_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *n_ghost_layers*n_grid_bound_j);
    }
    if (neighbors_id[4] != -1) {
        bkn_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *n_ghost_layers*n_grid_bound_k);
        bkn_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *n_ghost_layers*n_grid_bound_k);
    }
    if (neighbors_id[5] != -1) {
        bkp_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *n_ghost_layers*n_grid_bound_k);
        bkp_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *n_ghost_layers*n_grid_bound_k);
    }
    if (neighbors_id_ij[0] != -1) {
        bij_nn_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_K_vis);
        bij_nn_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_K_vis);
    }
    if (neighbors_id_ij[1] != -1) {
        bij_np_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_K_vis);
        bij_np_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_K_vis);
    }
    if (neighbors_id_ij[2] != -1) {
        bij_pn_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_K_vis);
        bij_pn_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_K_vis);
    }
    if (neighbors_id_ij[3] != -1) {
        bij_pp_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_K_vis);
        bij_pp_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_K_vis);
    }
    if (neighbors_id_jk[0] != -1) {
        bjk_nn_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_I_vis);
        bjk_nn_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_I_vis);
    }
    if (neighbors_id_jk[1] != -1) {
        bjk_np_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_I_vis);
        bjk_np_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_I_vis);
    }
    if (neighbors_id_jk[2] != -1) {
        bjk_pn_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_I_vis);
        bjk_pn_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_I_vis);
    }
    if (neighbors_id_jk[3] != -1) {
        bjk_pp_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_I_vis);
        bjk_pp_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_I_vis);
    }
    if (neighbors_id_ik[0] != -1) {
        bik_nn_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_J_vis);
        bik_nn_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_J_vis);
    }
    if (neighbors_id_ik[1] != -1) {
        bik_np_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_J_vis);
        bik_np_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_J_vis);
    }
    if (neighbors_id_ik[2] != -1) {
        bik_pn_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_J_vis);
        bik_pn_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_J_vis);
    }
    if (neighbors_id_ik[3] != -1) {
        bik_pp_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_J_vis);
        bik_pp_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *loc_J_vis);
    }
    if (neighbors_id_ijk[0] != -1) {
        bijk_nnn_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *1);
        bijk_nnn_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *1);
    }
    if (neighbors_id_ijk[1] != -1) {
        bijk_nnp_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *1);
        bijk_nnp_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *1);
    }
    if (neighbors_id_ijk[2] != -1) {
        bijk_npn_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *1);
        bijk_npn_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *1);
    }
    if (neighbors_id_ijk[3] != -1) {
        bijk_npp_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *1);
        bijk_npp_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *1);
    }
    if (neighbors_id_ijk[4] != -1) {
        bijk_pnn_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *1);
        bijk_pnn_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *1);
    }
    if (neighbors_id_ijk[5] != -1) {
        bijk_pnp_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *1);
        bijk_pnp_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *1);
    }
    if (neighbors_id_ijk[6] != -1) {
        bijk_ppn_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *1);
        bijk_ppn_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *1);
    }
    if (neighbors_id_ijk[7] != -1) {
        bijk_ppp_s = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *1);
        bijk_ppp_r = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *1);
    }

    // array for mpi request
    mpi_reqs = (MPI_Request *) malloc(sizeof(MPI_Request) * 6);
    for (int i = 0; i < 6; i++)
        mpi_reqs[i] = MPI_REQUEST_NULL;
    mpi_reqs_kosumi = (MPI_Request *) malloc(sizeof(MPI_Request) * 20);
    for (int i = 0; i < 20; i++)
        mpi_reqs_kosumi[i] = MPI_REQUEST_NULL;

    // velocity model
    //velo_loc = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);

    // arrays for data output
    int nnodes_loc_vis =  loc_I_vis    *  loc_J_vis    *  loc_K_vis;
    int nelms_loc_vis  = (loc_I_vis-1) * (loc_J_vis-1) * (loc_K_vis-1);
    x_loc_3d     = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL *) *nnodes_loc_vis);
    y_loc_3d     = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL *) *nnodes_loc_vis);
    z_loc_3d     = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL *) *nnodes_loc_vis);
    p_loc_3d     = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL *) *nnodes_loc_vis);
    t_loc_3d     = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL *) *nnodes_loc_vis);
    r_loc_3d     = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL *) *nnodes_loc_vis);
    elms_conn    = (int *) malloc(sizeof(int *) *nelms_loc_vis*9);
    my_proc_dump = (int *) malloc(sizeof(int *) *nnodes_loc_vis); // DEBUG
    vis_data     = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) *nnodes_loc_vis); // temp array for visuzulization

    // arrays for inversion
    if (inverse_flag) {
        int n_total_loc_inv_grid = n_inv_I_loc * n_inv_J_loc * n_inv_K_loc;
        r_loc_inv       = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_inv_K_loc * n_inv_grids);
        t_loc_inv       = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_inv_J_loc * n_inv_grids);
        p_loc_inv       = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_inv_I_loc * n_inv_grids);
        Ks_loc          = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
        Ks_inv_loc      = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_inv_grid);
        Ks_update_loc   = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
        if (sub_nprocs <= 1)
            Tadj_loc = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);

        if (optim_method==HALVE_STEPPING_MODE) {
            fac_b_loc_back = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
            fac_c_loc_back = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
            fac_f_loc_back = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
            fun_loc_back = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
        }

        if (optim_method==LBFGS_MODE) {
            fac_b_loc_back = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
            fac_c_loc_back = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
            fac_f_loc_back = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);
            fun_loc_back = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_grid_points);

            int n_total_loc_lbfgs = n_total_loc_grid_points * Mbfgs;
            Ks_grad_store_loc    = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_lbfgs);
            Ks_model_store_loc   = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_lbfgs);
            Ks_descent_dir_loc   = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_lbfgs);
            Ks_regularization_penalty_loc = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_lbfgs);
            fun_regularization_penalty_loc = (CUSTOMREAL *) malloc(sizeof(CUSTOMREAL) * n_total_loc_lbfgs);

            // initialize
            for (int i = 0; i < n_total_loc_lbfgs; i++) {
                Ks_grad_store_loc[i]    = _0_CR;
                Ks_model_store_loc[i]   = _0_CR;
                Ks_descent_dir_loc[i]   = _0_CR;
                fun_regularization_penalty_loc[i] = _0_CR;
                Ks_regularization_penalty_loc[i]   = _0_CR;
            }
        }
    } // end of if inverse_flag

    stdout_by_main("Memory allocation done.");
}


void Grid::shm_memory_allocation() {

    synchronize_all_sub();

    int n_total_loc_grid_points = loc_I * loc_J * loc_K;
    if (sub_rank != 0)
        n_total_loc_grid_points = 0;

    prepare_shm_array_cr(n_total_loc_grid_points, T_loc, win_T_loc);
    prepare_shm_array_cr(n_total_loc_grid_points, tau_old_loc, win_tau_old_loc);

    prepare_shm_array_cr(n_total_loc_grid_points, T0r_loc, win_T0r_loc);
    prepare_shm_array_cr(n_total_loc_grid_points, T0t_loc, win_T0t_loc);
    prepare_shm_array_cr(n_total_loc_grid_points, T0p_loc, win_T0p_loc);
    prepare_shm_array_cr(n_total_loc_grid_points, T0v_loc, win_T0v_loc);
    prepare_shm_array_cr(n_total_loc_grid_points, tau_loc, win_tau_loc);
    prepare_shm_array_cr(n_total_loc_grid_points, fac_a_loc, win_fac_a_loc);
    prepare_shm_array_cr(n_total_loc_grid_points, fac_b_loc, win_fac_b_loc);
    prepare_shm_array_cr(n_total_loc_grid_points, fac_c_loc, win_fac_c_loc);
    prepare_shm_array_cr(n_total_loc_grid_points, fac_f_loc, win_fac_f_loc);
    prepare_shm_array_cr(n_total_loc_grid_points, fun_loc, win_fun_loc);
    prepare_shm_array_bool(n_total_loc_grid_points, is_changed, win_is_changed);

    prepare_shm_array_cr(loc_I, p_loc_1d, win_p_loc_1d);
    prepare_shm_array_cr(loc_J, t_loc_1d, win_t_loc_1d);
    prepare_shm_array_cr(loc_K, r_loc_1d, win_r_loc_1d);

    // inversion
    prepare_shm_array_cr(n_total_loc_grid_points, Tadj_loc, win_Tadj_loc);

}


void Grid::shm_memory_deallocation() {
    // add if necessary
}


// function for memory allocation, called only for subdomain.
void Grid::memory_deallocation() {
   // shm array deallocation is done in shm_memory_deallocation() if shm is used
    if (sub_nprocs <= 1) {

#ifdef USE_CUDA
        if(use_gpu)
            cudaFree(tau_loc);
        else
            free(tau_loc);
#else
        free(tau_loc);
#endif
        free(T_loc);
        free(tau_old_loc);

        free(T0v_loc);
        free(T0r_loc);
        free(T0t_loc);
        free(T0p_loc);
        free(fac_a_loc);
        free(fac_b_loc);
        free(fac_c_loc);
        free(fac_f_loc);
        free(fun_loc);
        free(is_changed);

        free(t_loc_1d);
        free(p_loc_1d);
        free(r_loc_1d);
    }
    if(if_test) free(u_loc);

    if (neighbors_id[0] != -1) {
        free(bin_s);
        free(bin_r);
    }
    if (neighbors_id[1] != -1) {
        free(bip_s);
        free(bip_r);
    }

    if (neighbors_id[2] != -1) {
        free(bjn_s);
        free(bjn_r);
    }
    if (neighbors_id[3] != -1) {
        free(bjp_s);
        free(bjp_r);
    }
    if (neighbors_id[4] != -1) {
        free(bkn_s);
        free(bkn_r);
    }
    if (neighbors_id[5] != -1) {
        free(bkp_s);
        free(bkp_r);
    }
    if (neighbors_id_ij[0] != -1) {
        free(bij_nn_s);
        free(bij_nn_r);
    }
    if (neighbors_id_ij[1] != -1) {
        free(bij_np_s);
        free(bij_np_r);
    }
    if (neighbors_id_ij[2] != -1) {
        free(bij_pn_s);
        free(bij_pn_r);
    }
    if (neighbors_id_ij[3] != -1) {
        free(bij_pp_s);
        free(bij_pp_r);
    }
    if (neighbors_id_jk[0] != -1) {
        free(bjk_nn_s);
        free(bjk_nn_r);
    }
    if (neighbors_id_jk[1] != -1) {
        free(bjk_np_s);
        free(bjk_np_r);
    }
    if (neighbors_id_jk[2] != -1) {
        free(bjk_pn_s);
        free(bjk_pn_r);
    }
    if (neighbors_id_jk[3] != -1) {
        free(bjk_pp_s);
        free(bjk_pp_r);
    }
    if (neighbors_id_ik[0] != -1) {
        free(bik_nn_s);
        free(bik_nn_r);
    }
    if (neighbors_id_ik[1] != -1) {
        free(bik_np_s);
        free(bik_np_r);
    }
    if (neighbors_id_ik[2] != -1) {
        free(bik_pn_s);
        free(bik_pn_r);
    }
    if (neighbors_id_ik[3] != -1) {
        free(bik_pp_s);
        free(bik_pp_r);
    }
    if (neighbors_id_ijk[0] != -1) {
        free(bijk_nnn_s);
        free(bijk_nnn_r);
    }
    if (neighbors_id_ijk[1] != -1) {
        free(bijk_nnp_s);
        free(bijk_nnp_r);
    }
    if (neighbors_id_ijk[2] != -1) {
        free(bijk_npn_s);
        free(bijk_npn_r);
    }
    if (neighbors_id_ijk[3] != -1) {
        free(bijk_npp_s);
        free(bijk_npp_r);
    }
    if (neighbors_id_ijk[4] != -1) {
        free(bijk_pnn_s);
        free(bijk_pnn_r);
    }
    if (neighbors_id_ijk[5] != -1) {
        free(bijk_pnp_s);
        free(bijk_pnp_r);
    }
    if (neighbors_id_ijk[6] != -1) {
        free(bijk_ppn_s);
        free(bijk_ppn_r);
    }
    if (neighbors_id_ijk[7] != -1) {
        free(bijk_ppp_s);
        free(bijk_ppp_r);
    }

    free(mpi_reqs);
    free(mpi_reqs_kosumi);

    // delete arrays
    delete [] nnodes;
    delete [] nelms;
    delete [] nnodes_vis;
    delete [] nelms_vis;

    //free(velo_loc);

    free(x_loc_3d);
    free(y_loc_3d);
    free(z_loc_3d);
    free(p_loc_3d);
    free(t_loc_3d);
    free(r_loc_3d);
    free(elms_conn);
    free(my_proc_dump);
    free(vis_data);

    // inversion arrays
    if (inverse_flag){
        free(r_loc_inv);
        free(p_loc_inv);
        free(t_loc_inv);
        free(Ks_loc);
        free(Ks_inv_loc);
        free(Ks_update_loc);

        if (sub_nprocs <= 1)
            free(Tadj_loc);

        if (optim_method==HALVE_STEPPING_MODE) {
            free(fac_b_loc_back);
            free(fac_c_loc_back);
            free(fac_f_loc_back);
            free(fun_loc_back);
            free(xi_loc_back);
            free(eta_loc_back);
        }

        if (optim_method==LBFGS_MODE) {
            free(fun_loc_back);
            free(xi_loc_back);
            free(eta_loc_back);

            free(Ks_grad_store_loc   );
            free(Ks_model_store_loc  );
            free(Ks_descent_dir_loc);
            free(fun_regularization_penalty_loc);
            free(Ks_regularization_penalty_loc);
        }
    } // end if inverse_flag

    stdout_by_main("Memory deallocation done.");
}


// setput the grid parameters (material and source etc.)
void Grid::setup_grid_params(InputParams &IP, IO_utils& io) {

    // setup coordinates
    r_min   = depth2radius(IP.get_max_dep()); // convert from depth to radius
    r_max   = depth2radius(IP.get_min_dep()); // convert from depth to radius
    lat_min = IP.get_min_lat();
    lat_max = IP.get_max_lat();
    lon_min = IP.get_min_lon();
    lon_max = IP.get_max_lon();

    dr   = (r_max - r_min)     / (ngrid_k - 1);
    dlat = (lat_max - lat_min) / (ngrid_j - 1);
    dlon = (lon_max - lon_min) / (ngrid_i - 1);

    dt = dlat; dp = dlon;

    // starting position in global indices including ghost layers
    int tmp_offset_k = get_offset_k();
    int tmp_offset_j = get_offset_j();
    int tmp_offset_i = get_offset_i(); // offset_i - i_start_loc

    // output global grid parameters
    if (inter_sub_rank == 0 && myrank == 0) {
        std::cout << "\n\n Global grid information: \n\n" << std::endl;
        std::cout << "  r_min r_max dr : "                << r_min   << " " << r_max << " " << dr << std::endl;
        std::cout << "  depth_min depth_max dz : "        << IP.get_min_dep() << " " << IP.get_max_dep() << std::endl;
        std::cout << "  lat_min lat_max dlat [degree] : " << lat_min*RAD2DEG << " " << lat_max*RAD2DEG << " " << dlat*RAD2DEG << std::endl;
        std::cout << "  lon_min lon_max dlon [degree] : " << lon_min*RAD2DEG << " " << lon_max*RAD2DEG << " " << dlon*RAD2DEG << std::endl;
        std::cout << "  ngrid_i ngrid_j ngrid_k : "      << ngrid_i << " " << ngrid_j << " " << ngrid_k << std::endl;
        std::cout << "\n\n" << std::endl;
    }

    // assign coordinates to 1d arrays
    // r
    for (int k = 0; k < loc_K; k++) {
        r_loc_1d[k] = r_min + (tmp_offset_k + k)*dr;
    }
    // lat
    for (int j = 0; j < loc_J; j++) {
        t_loc_1d[j] = lat_min + (tmp_offset_j + j)*dlat;
    }
    // lon
    for (int i = 0; i < loc_I; i++) {
        p_loc_1d[i] = lon_min + (tmp_offset_i + i)*dlon;
    }

    for (int i = 0; i < nprocs; i++) {
        synchronize_all_inter();
        if (i == myrank) {
            std::cout << "\n--- subdomain info ---- rank : " << i << "\n\n" << std::endl;
            std::cout << "  p_loc_1d min max in deg.: "  << p_loc_1d[0]*RAD2DEG << " " << p_loc_1d[loc_I-1]*RAD2DEG << std::endl;
            std::cout << "  t_loc_1d min max in deg.: "  << t_loc_1d[0]*RAD2DEG << " " << t_loc_1d[loc_J-1]*RAD2DEG << std::endl;
            std::cout << "  r_loc_1d min max in km  : "  << r_loc_1d[0]         << " " << r_loc_1d[loc_K-1]         << std::endl;

            std::cout << "  p_loc_1d min max in rad.: "  << p_loc_1d[0] << " " << p_loc_1d[loc_I-1] << std::endl;
            std::cout << "  t_loc_1d min max in rad.: "  << t_loc_1d[0] << " " << t_loc_1d[loc_J-1] << std::endl;

            std::cout << "  r_min, get_offset_k(), dr in km: "     << r_min           << " " << tmp_offset_k << " " << dr           << std::endl;
            std::cout << "  lat_min, get_offset_j(), dlat in deg.: " << lat_min*RAD2DEG << " " << tmp_offset_j << " " << dlat*RAD2DEG << std::endl;
            std::cout << "  lon_min, get_offset_i(), dlon in deg.: " << lon_min*RAD2DEG << " " << tmp_offset_i << " " << dlon*RAD2DEG << std::endl;
            std::cout << "\n" << std::endl;
        }
    }

    // independent read model data
    if (id_sim == 0){
        // read init model
        std::string f_model_path = IP.get_init_model_path();
        io.read_model(f_model_path,"vel",  fun_loc,   tmp_offset_i, tmp_offset_j, tmp_offset_k); // use slowness array temprarily
        if(if_test) {
            // solver test
            io.read_model(f_model_path, "u", u_loc, tmp_offset_i, tmp_offset_j, tmp_offset_k);
        }
    }

    // broadcast
    int n_total_loc_grid_points = loc_I * loc_J * loc_K;
    broadcast_cr_inter_sim(fun_loc,   n_total_loc_grid_points, 0); // here passing velocity array
    if(if_test) {
        broadcast_cr_inter_sim(u_loc, n_total_loc_grid_points, 0);
    }

    // center of the domain
    CUSTOMREAL lon_center = (lon_min + lon_max) / 2.0;
    CUSTOMREAL lat_center = (lat_min + lat_max) / 2.0;

    // initialize other necessary arrays
    for (int k_r = 0; k_r < loc_K; k_r++) {
        for (int j_lat = 0; j_lat < loc_J; j_lat++) {
            for (int i_lon = 0; i_lon < loc_I; i_lon++) {

                // convert velocity to slowness
                fun_loc[I2V(i_lon, j_lat, k_r)] = _1_CR / fun_loc[I2V(i_lon, j_lat, k_r)];

                // calculate fac_a_loc, fac_b_loc, fac_c_loc, fac_f_loc
                fac_a_loc[I2V(i_lon, j_lat, k_r)] = _1_CR;
                fac_b_loc[I2V(i_lon, j_lat, k_r)] = _1_CR;
                fac_c_loc[I2V(i_lon, j_lat, k_r)] = _1_CR;
                fac_f_loc[I2V(i_lon, j_lat, k_r)] = _0_CR;

                // construct 3d coordinate arrays and node connectivity for visualization
                // exclude the ghost nodes
                if (i_start_vis <= i_lon && i_lon <= i_end_vis \
                 && j_start_vis <= j_lat && j_lat <= j_end_vis \
                 && k_start_vis <= k_r   && k_r   <= k_end_vis) {
                    int i_loc_tmp = i_lon - i_start_vis;
                    int j_loc_tmp = j_lat - j_start_vis;
                    int k_loc_tmp = k_r   - k_start_vis;

                    // get the coordinates of the point
                    CUSTOMREAL x, y, z;
                    //RLonLat2xyz(p_loc_1d[i_lon], t_loc_1d[j_lat], r_loc_1d[k_r], x, y, z); // complete sphere
                    WGS84ToCartesian(p_loc_1d[i_lon], t_loc_1d[j_lat], r_loc_1d[k_r], x, y, z, lon_center, lat_center); // WGS84

                    x_loc_3d[I2V_VIS(i_loc_tmp,j_loc_tmp,k_loc_tmp)] = x;
                    y_loc_3d[I2V_VIS(i_loc_tmp,j_loc_tmp,k_loc_tmp)] = y;
                    z_loc_3d[I2V_VIS(i_loc_tmp,j_loc_tmp,k_loc_tmp)] = z;
                    p_loc_3d[I2V_VIS(i_loc_tmp,j_loc_tmp,k_loc_tmp)] = p_loc_1d[i_lon]*RAD2DEG;
                    t_loc_3d[I2V_VIS(i_loc_tmp,j_loc_tmp,k_loc_tmp)] = t_loc_1d[j_lat]*RAD2DEG;
                    // r_loc_3d[I2V_VIS(i_loc_tmp,j_loc_tmp,k_loc_tmp)] = r_loc_1d[k_r]*M2KM;
                    r_loc_3d[I2V_VIS(i_loc_tmp,j_loc_tmp,k_loc_tmp)] = r_loc_1d[k_r];

                    // dump proc id
                    my_proc_dump[I2V_VIS(i_loc_tmp,j_loc_tmp,k_loc_tmp)] = myrank;

                    // node connectivity
                    // exclude the last layers
                    if (i_lon < i_end_vis \
                     && j_lat < j_end_vis \
                     && k_r   < k_end_vis) {
                        elms_conn[I2V_ELM_CONN(i_loc_tmp,j_loc_tmp,k_loc_tmp)*9] = cell_type;
                        elms_conn[I2V_ELM_CONN(i_loc_tmp,j_loc_tmp,k_loc_tmp)*9+1] = offset_nnodes_vis+I2V_VIS(i_loc_tmp,  j_loc_tmp,  k_loc_tmp);
                        elms_conn[I2V_ELM_CONN(i_loc_tmp,j_loc_tmp,k_loc_tmp)*9+2] = offset_nnodes_vis+I2V_VIS(i_loc_tmp+1,j_loc_tmp,  k_loc_tmp);
                        elms_conn[I2V_ELM_CONN(i_loc_tmp,j_loc_tmp,k_loc_tmp)*9+3] = offset_nnodes_vis+I2V_VIS(i_loc_tmp+1,j_loc_tmp+1,k_loc_tmp);
                        elms_conn[I2V_ELM_CONN(i_loc_tmp,j_loc_tmp,k_loc_tmp)*9+4] = offset_nnodes_vis+I2V_VIS(i_loc_tmp,  j_loc_tmp+1,k_loc_tmp);
                        elms_conn[I2V_ELM_CONN(i_loc_tmp,j_loc_tmp,k_loc_tmp)*9+5] = offset_nnodes_vis+I2V_VIS(i_loc_tmp,  j_loc_tmp,  k_loc_tmp+1);
                        elms_conn[I2V_ELM_CONN(i_loc_tmp,j_loc_tmp,k_loc_tmp)*9+6] = offset_nnodes_vis+I2V_VIS(i_loc_tmp+1,j_loc_tmp,  k_loc_tmp+1);
                        elms_conn[I2V_ELM_CONN(i_loc_tmp,j_loc_tmp,k_loc_tmp)*9+7] = offset_nnodes_vis+I2V_VIS(i_loc_tmp+1,j_loc_tmp+1,k_loc_tmp+1);
                        elms_conn[I2V_ELM_CONN(i_loc_tmp,j_loc_tmp,k_loc_tmp)*9+8] = offset_nnodes_vis+I2V_VIS(i_loc_tmp,  j_loc_tmp+1,k_loc_tmp+1);
                    }
                }
            }
        }
    } // end of for loop

    // setup parameters for inversion grids
    if (inverse_flag)
        setup_inv_grid_params(IP);

}


void Grid::setup_inv_grid_params(InputParams& IP) {

    // inversion grid for depth
    if(IP.get_type_dep_inv() == 0){     // uniform inversion grid for depth
        CUSTOMREAL r_min_inv   = depth2radius(IP.get_max_dep_inv()); // convert from depth to radius
        CUSTOMREAL r_max_inv   = depth2radius(IP.get_min_dep_inv()); // convert from depth to radius
        // inversion grid is defined for all processes which covers entire domain
        dinv_r = (r_max_inv   - r_min_inv)   / (ngrid_k_inv-1);
        // shift of each set of inversion grid
        dinv_lr = dinv_r/n_inv_grids;

        for (int l = 0; l < n_inv_grids; l++) {
            for (int k = 0; k < n_inv_K_loc; k++)
                r_loc_inv[I2V_INV_GRIDS_1DK(k,l)] = r_min_inv   + k*dinv_r - l*dinv_lr;
        }
    } else {            // flexibly designed inversion grid for depth
        CUSTOMREAL* dep_inv = IP.get_dep_inv();

        for (int k = 0; k < n_inv_K_loc; k++){
            if(k < n_inv_K_loc - 1)
                dinv_lr = (depth2radius(dep_inv[k+1]) - depth2radius(dep_inv[k]))/n_inv_grids;
            else
                dinv_lr = (depth2radius(dep_inv[n_inv_K_loc-1]) - depth2radius(dep_inv[n_inv_K_loc-2]))/n_inv_grids;

            for (int l = 0; l < n_inv_grids; l++)
                r_loc_inv[I2V_INV_GRIDS_1DK(k,l)] = depth2radius(dep_inv[k]) + l*dinv_lr;
        }
    }

    // inversion grid for latitude
    if(IP.get_type_lat_inv() == 0){     // uniform inversion grid for latitude
        CUSTOMREAL lat_min_inv = IP.get_min_lat_inv();
        CUSTOMREAL lat_max_inv = IP.get_max_lat_inv();
        // inversion grid is defined for all processes which covers entire domain
        dinv_t = (lat_max_inv - lat_min_inv) / (ngrid_j_inv-1);
        // shift of each set of inversion grid
        dinv_lt = dinv_t/n_inv_grids;

        for (int l = 0; l < n_inv_grids; l++) {
            for (int j = 0; j < n_inv_J_loc; j++)
                t_loc_inv[I2V_INV_GRIDS_1DJ(j,l)] = lat_min_inv + j*dinv_t - l*dinv_lt;
        }
    } else {                // flexibly designed inversion grid for latitude
        CUSTOMREAL* lat_inv = IP.get_lat_inv();

        for (int j = 0; j < n_inv_J_loc; j++){
            if(j < n_inv_J_loc - 1)
                dinv_lt = (lat_inv[j+1] - lat_inv[j])*DEG2RAD/n_inv_grids;
            else
                dinv_lt = (lat_inv[n_inv_J_loc-1] - lat_inv[n_inv_J_loc-2])*DEG2RAD/n_inv_grids;

            for (int l = 0; l < n_inv_grids; l++)
                t_loc_inv[I2V_INV_GRIDS_1DJ(j,l)] = lat_inv[j]*DEG2RAD - l*dinv_lt;
        }
    }

    // inversion grid for longitude
    if(IP.get_type_lon_inv() == 0){     // uniform inversion grid for longitude
        CUSTOMREAL lon_min_inv = IP.get_min_lon_inv();
        CUSTOMREAL lon_max_inv = IP.get_max_lon_inv();
        // inversion grid is defined for all processes which covers entire domain
        dinv_p = (lon_max_inv - lon_min_inv) / (ngrid_i_inv-1);
        // shift of each set of inversion grid
        dinv_lp = dinv_p/n_inv_grids;

        for (int l = 0; l < n_inv_grids; l++) {
            for (int i = 0; i < n_inv_I_loc; i++)
                p_loc_inv[I2V_INV_GRIDS_1DI(i,l)] = lon_min_inv + i*dinv_p - l*dinv_lp;
        }
    } else {
        CUSTOMREAL* lon_inv = IP.get_lon_inv();

        for (int i = 0; i < n_inv_I_loc; i++){
            if(i < n_inv_I_loc - 1)
                dinv_lp = (lon_inv[i+1] - lon_inv[i])*DEG2RAD/n_inv_grids;
            else
                dinv_lp = (lon_inv[n_inv_I_loc-1] - lon_inv[n_inv_I_loc-2])*DEG2RAD/n_inv_grids;

            for (int l = 0; l < n_inv_grids; l++)
                p_loc_inv[I2V_INV_GRIDS_1DI(i,l)] = lon_inv[i]*DEG2RAD - l*dinv_lp;
        }
    }


}


void Grid::initialize_kernels(){
    // initialize kernels
    if (subdom_main){
        for (int k = 0; k < loc_K; k++) {
            for (int j = 0; j < loc_J; j++) {
                for (int i = 0; i < loc_I; i++) {
                    Ks_loc[  I2V(i,j,k)] = _0_CR;
                }
            }
        }
    }
}


// get a part of pointers from the requested array for visualization
CUSTOMREAL* Grid::get_array_for_vis(CUSTOMREAL* arr) {

    //send_recev_boundary_data(arr);
    // add a routine for communication the boundary value
    // with the neighbors with line/point contact
    send_recev_boundary_data_kosumi(arr);

    for (int k_r = 0; k_r < loc_K; k_r++) {
        for (int j_lat = 0; j_lat < loc_J; j_lat++) {
            for (int i_lon = 0; i_lon < loc_I; i_lon++) {
                 if (i_start_vis <= i_lon && i_lon <= i_end_vis \
                  && j_start_vis <= j_lat && j_lat <= j_end_vis \
                  && k_start_vis <= k_r   && k_r   <= k_end_vis) {
                    int i_loc_tmp = i_lon - i_start_vis;
                    int j_loc_tmp = j_lat - j_start_vis;
                    int k_loc_tmp = k_r   - k_start_vis;

                    vis_data[I2V_VIS(i_loc_tmp,j_loc_tmp,k_loc_tmp)] = arr[I2V(i_lon,j_lat,k_r)];
                 }
            }
        }
    }

    return vis_data;
}


// set visualization array to local array
void Grid::set_array_from_vis(CUSTOMREAL* arr) {

    for (int k_r = 0; k_r < loc_K; k_r++) {
        for (int j_lat = 0; j_lat < loc_J; j_lat++) {
            for (int i_lon = 0; i_lon < loc_I; i_lon++) {
                 if (i_start_vis <= i_lon && i_lon <= i_end_vis \
                  && j_start_vis <= j_lat && j_lat <= j_end_vis \
                  && k_start_vis <= k_r   && k_r   <= k_end_vis) {
                    int i_loc_tmp = i_lon - i_start_vis;
                    int j_loc_tmp = j_lat - j_start_vis;
                    int k_loc_tmp = k_r   - k_start_vis;

                    arr[I2V(i_lon,j_lat,k_r)] = vis_data[I2V_VIS(i_loc_tmp,j_loc_tmp,k_loc_tmp)];
                 }
            }
        }
    }

    // set values to ghost layer
    send_recev_boundary_data(arr);

}


void Grid::reinitialize_abcf(){
    if (subdom_main) {
        for (int k_r = 0; k_r < loc_K; k_r++) {
            for (int j_lat = 0; j_lat < loc_J; j_lat++) {
                for (int i_lon = 0; i_lon < loc_I; i_lon++) {
                    // initialize arrays
                    fac_a_loc[I2V(i_lon,j_lat,k_r)] = fac_a_loc[I2V(i_lon,j_lat,k_r)];
                    fac_b_loc[I2V(i_lon,j_lat,k_r)] = fac_b_loc[I2V(i_lon,j_lat,k_r)]/ my_square(r_loc_1d[k_r]);
                    fac_c_loc[I2V(i_lon,j_lat,k_r)] = fac_c_loc[I2V(i_lon,j_lat,k_r)]/(my_square(r_loc_1d[k_r])*my_square(std::cos(t_loc_1d[j_lat])));
                    fac_f_loc[I2V(i_lon,j_lat,k_r)] = fac_f_loc[I2V(i_lon,j_lat,k_r)]/(my_square(r_loc_1d[k_r])*          std::cos(t_loc_1d[j_lat]));
                }
            }
        }
    }
}


void Grid::setup_factors(Source &src){

    // calculate factors for the source
    a0   = src.get_fac_at_source(fac_a_loc);
    b0   = src.get_fac_at_source(fac_b_loc);
    c0   = src.get_fac_at_source(fac_c_loc);
    f0   = src.get_fac_at_source(fac_f_loc);
    fun0 = src.get_fac_at_source(fun_loc);

}


void Grid::initialize_fields(Source& src, InputParams& IP){

    // get source position
    CUSTOMREAL src_r = src.get_src_r();
    CUSTOMREAL src_t = src.get_src_t();
    CUSTOMREAL src_p = src.get_src_p();

    // std out src positions
    CUSTOMREAL c0b0_minus_f0f0 = c0*b0 - f0*f0;

    // debug
    int n_source_node = 0;

    // std::cout << a0 << ' ' << b0 << ' ' << c0 << ' ' << f0 << ' ' << std::endl;



    for (int k_r = 0; k_r < loc_K; k_r++) {
        for (int j_lat = 0; j_lat < loc_J; j_lat++) {
            for (int i_lon = 0; i_lon < loc_I; i_lon++) {
                CUSTOMREAL dr_from_src = r_loc_1d[k_r]   - src_r;
                CUSTOMREAL dt_from_src = t_loc_1d[j_lat] - src_t;
                CUSTOMREAL dp_from_src = p_loc_1d[i_lon] - src_p;

                T0v_loc[I2V(i_lon,j_lat,k_r)] = fun0 * std::sqrt( _1_CR/a0                  *my_square(dr_from_src) \
                                                                + c0/(c0b0_minus_f0f0)      *my_square(dt_from_src) \
                                                                + b0/(c0b0_minus_f0f0)      *my_square(dp_from_src) \
                                                                + _2_CR*f0/(c0b0_minus_f0f0)*dt_from_src*dp_from_src);

                is_changed[I2V(i_lon,j_lat,k_r)] = true;

                if (isZero(T0v_loc[I2V(i_lon,j_lat,k_r)])) {
                    T0r_loc[I2V(i_lon,j_lat,k_r)] = _0_CR;
                    T0t_loc[I2V(i_lon,j_lat,k_r)] = _0_CR;
                    T0p_loc[I2V(i_lon,j_lat,k_r)] = _0_CR;
                } else {
                    T0r_loc[I2V(i_lon,j_lat,k_r)] = my_square(fun0)*(_1_CR/a0*dr_from_src)/T0v_loc[I2V(i_lon,j_lat,k_r)];
                    T0t_loc[I2V(i_lon,j_lat,k_r)] = my_square(fun0)*(c0/(c0b0_minus_f0f0)*dt_from_src+f0/(c0b0_minus_f0f0)*dp_from_src)/T0v_loc[I2V(i_lon,j_lat,k_r)];
                    T0p_loc[I2V(i_lon,j_lat,k_r)] = my_square(fun0)*(b0/(c0b0_minus_f0f0)*dp_from_src+f0/(c0b0_minus_f0f0)*dt_from_src)/T0v_loc[I2V(i_lon,j_lat,k_r)];
                }

                if (IP.get_stencil_order() == 1){
                    source_width = _1_CR-0.1;
                } else {
                    source_width = _2_CR;
                }

                if (std::abs(dr_from_src/dr) <= source_width \
                 && std::abs(dt_from_src/dt) <= source_width \
                 && std::abs(dp_from_src/dp) <= source_width) {

                    tau_loc[I2V(i_lon,j_lat,k_r)] = TAU_INITIAL_VAL;
                    is_changed[I2V(i_lon,j_lat,k_r)] = false;

                    n_source_node++;

                    if (if_verbose) {
                        if ( (k_r == 0 && k_first()) || (k_r == loc_K-1 && k_last()) )
                            std::cout << "Warning: source is on the boundary k of the grid.\n";
                        if ( (j_lat == 0 && j_first()) || (j_lat == loc_J-1 && j_last()) )
                            std::cout << "Warning: source is on the boundary j of the grid.\n";
                        if ( (i_lon == 0 && i_first()) || (i_lon == loc_I-1 && i_last()) )
                            std::cout << "Warning: source is on the boundary i of the grid.\n";
                    }

                } else {
                    if (IP.get_stencil_type()==1)   // upwind scheme, initial tau should be large enough
                        tau_loc[I2V(i_lon,j_lat,k_r)] = TAU_INF_VAL;
                    else
                        tau_loc[I2V(i_lon,j_lat,k_r)] = TAU_INITIAL_VAL;
                    is_changed[I2V(i_lon,j_lat,k_r)] = true;

                }

                tau_old_loc[I2V(i_lon,j_lat,k_r)] = _0_CR;

            } // end loop i
        } // end loop j
    } // end loop k

    // warning if source node is not found
    if( n_source_node > 0 && if_verbose )
        std::cout << "rank  n_source_node: " << myrank << "  " << n_source_node << std::endl;

}


// copy the tau values to tau_old
void Grid::tau2tau_old() {
    std::memcpy(tau_old_loc, tau_loc, sizeof(CUSTOMREAL)*loc_I*loc_J*loc_K);
}


// copy the T values to tau_old
void Grid::T2tau_old() {
    std::memcpy(tau_old_loc, T_loc, sizeof(CUSTOMREAL)*loc_I*loc_J*loc_K);
}


// copy the tau values to Tadj
void Grid::update_Tadj() {
    std::memcpy(Tadj_loc, tau_loc, sizeof(CUSTOMREAL)*loc_I*loc_J*loc_K);
}


void Grid::back_up_fun_bcf() {
    std::memcpy(fun_loc_back, fun_loc, sizeof(CUSTOMREAL)*loc_I*loc_J*loc_K);
    std::memcpy(fac_b_loc_back, fac_b_loc, sizeof(CUSTOMREAL)*loc_I*loc_J*loc_K);
    std::memcpy(fac_c_loc_back, fac_c_loc, sizeof(CUSTOMREAL)*loc_I*loc_J*loc_K);
    std::memcpy(fac_f_loc_back, fac_f_loc, sizeof(CUSTOMREAL)*loc_I*loc_J*loc_K);
}


void Grid::restore_fun_bcf() {
    std::memcpy(fun_loc, fun_loc_back, sizeof(CUSTOMREAL)*loc_I*loc_J*loc_K);
    std::memcpy(fac_b_loc, fac_b_loc_back, sizeof(CUSTOMREAL)*loc_I*loc_J*loc_K);
    std::memcpy(fac_c_loc, fac_c_loc_back, sizeof(CUSTOMREAL)*loc_I*loc_J*loc_K);
    std::memcpy(fac_f_loc, fac_f_loc_back, sizeof(CUSTOMREAL)*loc_I*loc_J*loc_K);
}


void Grid::calc_L1_and_Linf_diff(CUSTOMREAL& L1_diff, CUSTOMREAL& Linf_diff) {
    if (subdom_main) {
        L1_diff   = 0.0;
        Linf_diff = 0.0;
        CUSTOMREAL obj_func_glob = 0.0;

        // calculate L1 error
        for (int k_r = k_start_loc; k_r <= k_end_loc; k_r++) {
            for (int j_lat = j_start_loc; j_lat <= j_end_loc; j_lat++) {
                #pragma omp simd reduction(+:L1_diff)
                for (int i_lon = i_start_loc; i_lon <= i_end_loc; i_lon++) {
                    L1_diff   +=                    std::abs(tau_loc[I2V(i_lon,j_lat,k_r)] - tau_old_loc[I2V(i_lon,j_lat,k_r)]) * T0v_loc[I2V(i_lon,j_lat,k_r)];
                    Linf_diff  = std::max(Linf_diff,std::abs(tau_loc[I2V(i_lon,j_lat,k_r)] - tau_old_loc[I2V(i_lon,j_lat,k_r)]) * T0v_loc[I2V(i_lon,j_lat,k_r)]);
                    // std::cout << R_earth-r_loc_3d[k_r] << ' ' << t_loc_3d[j_lat] << ' ' << p_loc_3d[i_lon] << ' ' << T0v_loc[I2V(i_lon,j_lat,k_r)] << ' ' << std::endl;
                }
            }
        }

        // sum up the values of all processes
        obj_func_glob=0.0;
        allreduce_cr_single(L1_diff, obj_func_glob);
        L1_diff = obj_func_glob/((ngrid_i-2)*(ngrid_j-2)*(ngrid_k-2));
        obj_func_glob=0.0;
        allreduce_cr_single_max(Linf_diff, obj_func_glob);
        Linf_diff = obj_func_glob;///(ngrid_i*ngrid_j*ngrid_k);
    }
}


void Grid::calc_L1_and_Linf_diff_adj(CUSTOMREAL& L1_diff, CUSTOMREAL& Linf_diff) {

    if (subdom_main) {
        //L1_diff   = 0.0;
        Linf_diff = 0.0;
        CUSTOMREAL obj_func_glob = 0.0;
        CUSTOMREAL adj_factor = 1.0;

        // calculate L1 error
        for (int k_r = k_start_loc; k_r <= k_end_loc; k_r++) {
            for (int j_lat = j_start_loc; j_lat <= j_end_loc; j_lat++) {
                #pragma omp simd
                for (int i_lon = i_start_loc; i_lon <= i_end_loc; i_lon++) {
                    // Adjoint simulation use only Linf
                    Linf_diff  = std::max(Linf_diff,std::abs(tau_loc[I2V(i_lon,j_lat,k_r)] - Tadj_loc[I2V(i_lon,j_lat,k_r)]));
                }
            }
        }

        // sum up the values of all processes
        obj_func_glob=0.0;
        allreduce_cr_single_max(Linf_diff, obj_func_glob);
        Linf_diff = obj_func_glob*adj_factor;///(ngrid_i*ngrid_j*ngrid_k);
    }
}


void Grid::calc_L1_and_Linf_error(CUSTOMREAL& L1_error, CUSTOMREAL& Linf_error) {
    L1_error   = 0.0;
    Linf_error = 0.0;
    CUSTOMREAL obj_func_glob = 0.0;

    // calculate L1 error
    for (int k_r = k_start_loc; k_r <= k_end_loc; k_r++) {
        for (int j_lat = j_start_loc; j_lat <= j_end_loc; j_lat++) {
            #pragma omp simd reduction(+:L1_error)
            for (int i_lon = i_start_loc; i_lon <= i_end_loc; i_lon++) {
                L1_error   +=                     std::abs(u_loc[I2V(i_lon,j_lat,k_r)] - T0v_loc[I2V(i_lon,j_lat,k_r)] * tau_loc[I2V(i_lon,j_lat,k_r)]);
                Linf_error  = std::max(Linf_error,std::abs(u_loc[I2V(i_lon,j_lat,k_r)] - T0v_loc[I2V(i_lon,j_lat,k_r)] * tau_loc[I2V(i_lon,j_lat,k_r)]));
            }
        }
    }

    // sum up the values of all processes
    obj_func_glob=0.0;
    allreduce_cr_single(L1_error, obj_func_glob);
    L1_error       = obj_func_glob;
    obj_func_glob = 0.0;
    allreduce_cr_single_max(Linf_error, obj_func_glob);
    Linf_error     = obj_func_glob;///(ngrid_i*ngrid_j*ngrid_k);
}


void Grid::prepare_boundary_data_to_send(CUSTOMREAL* arr) {
    // store the pointers to the ghost layer's elements to be sent / to receive
    // node order should be always smaller to larger

    // j-k plane negative
    if (neighbors_id[0] != -1) {
        for (int i = 0; i < n_ghost_layers; i++) {
            for (int k = k_start_loc; k <= k_end_loc; k++) {
                #pragma omp simd
                for (int j = j_start_loc; j <= j_end_loc; j++) {
                    bin_s[i*loc_J_excl_ghost*loc_K_excl_ghost + (k-k_start_loc)*loc_J_excl_ghost + (j-j_start_loc)] = arr[I2V(i+n_ghost_layers,j,k)];
                }
            }
        }
    }
    // j-k plane positive
    if (neighbors_id[1] != -1) {
        for (int i = 0; i < n_ghost_layers; i++) {
            for (int k = k_start_loc; k <= k_end_loc; k++) {
                #pragma omp simd
                for (int j = j_start_loc; j <= j_end_loc; j++) {
                    bip_s[i*loc_J_excl_ghost*loc_K_excl_ghost + (k-k_start_loc)*loc_J_excl_ghost + (j-j_start_loc)] = arr[I2V(loc_I-n_ghost_layers*2+i,j,k)];
                }
            }
        }
    }
    // i-k plane negative
    if (neighbors_id[2] != -1) {
        for (int j = 0; j < n_ghost_layers; j++) {
            for (int k = k_start_loc; k <= k_end_loc; k++) {
                #pragma omp simd
                for (int i = i_start_loc; i <= i_end_loc; i++) {
                    bjn_s[j*loc_I_excl_ghost*loc_K_excl_ghost + (k-k_start_loc)*loc_I_excl_ghost + (i-i_start_loc)] = arr[I2V(i,j+n_ghost_layers,k)];
                }
            }
        }
    }
    // i-k plane positive
    if (neighbors_id[3] != -1) {
        for (int j = 0; j < n_ghost_layers; j++) {
            for (int k = k_start_loc; k <= k_end_loc; k++) {
                #pragma omp simd
                for (int i = i_start_loc; i <= i_end_loc; i++) {
                    bjp_s[j*loc_I_excl_ghost*loc_K_excl_ghost + (k-k_start_loc)*loc_I_excl_ghost + (i-i_start_loc)] = arr[I2V(i,loc_J-n_ghost_layers*2+j,k)];
                }
            }
        }
    }

    // i-j plane negative
    if (neighbors_id[4] != -1) {
        for (int k = 0; k < n_ghost_layers; k++) {
            for (int j = j_start_loc; j <= j_end_loc; j++) {
                #pragma omp simd
                for (int i = i_start_loc; i <= i_end_loc; i++) {
                    bkn_s[k*loc_I_excl_ghost*loc_J_excl_ghost + (j-j_start_loc)*loc_I_excl_ghost + (i-i_start_loc)] = arr[I2V(i,j,k+n_ghost_layers)];
                }
            }
        }
    }
    // i-j plane positive
    if (neighbors_id[5] != -1) {
        for (int k = 0; k < n_ghost_layers; k++) {
            for (int j = j_start_loc; j <= j_end_loc; j++) {
                #pragma omp simd
                for (int i = i_start_loc; i <= i_end_loc; i++) {
                    bkp_s[k*loc_I_excl_ghost*loc_J_excl_ghost + (j-j_start_loc)*loc_I_excl_ghost + (i-i_start_loc)] = arr[I2V(i,j,loc_K-n_ghost_layers*2+k)];
                }
            }
        }
    }

    //stdout_by_main("Preparation ghost layers done.");
}

void Grid::assign_received_data_to_ghost(CUSTOMREAL* arr) {
    // store the pointers to the ghost layer's elements to be sent / to receive
    // node order should be always smaller to larger

    // j-k plane negative
    if (neighbors_id[0] != -1) {
        for (int i = 0; i < n_ghost_layers; i++) {
            for (int k = k_start_loc; k <= k_end_loc; k++) {
                for (int j = j_start_loc; j <= j_end_loc; j++) {
                    arr[I2V(i,j,k)] = bin_r[i*loc_J_excl_ghost*loc_K_excl_ghost + (k-k_start_loc)*loc_J_excl_ghost + (j-j_start_loc)];
                }
            }
        }
    }
    // j-k plane positive
    if (neighbors_id[1] != -1) {
        for (int i = 0; i < n_ghost_layers; i++) {
            for (int k = k_start_loc; k <= k_end_loc; k++) {
                #pragma omp simd
                for (int j = j_start_loc; j <= j_end_loc; j++) {
                     arr[I2V(loc_I-n_ghost_layers+i,j,k)] = bip_r[i*loc_J_excl_ghost*loc_K_excl_ghost + (k-k_start_loc)*loc_J_excl_ghost + (j-j_start_loc)];
                }
            }
        }
    }
    // i-k plane negative
    if (neighbors_id[2] != -1) {
        for (int j = 0; j < n_ghost_layers; j++) {
            for (int k = k_start_loc; k <= k_end_loc; k++) {
                #pragma omp simd
                for (int i = i_start_loc; i <= i_end_loc; i++) {
                    arr[I2V(i,j,k)] = bjn_r[j*loc_I_excl_ghost*loc_K_excl_ghost + (k-k_start_loc)*loc_I_excl_ghost + (i-i_start_loc)];
                }
            }
        }
    }
    // i-k plane positive
    if (neighbors_id[3] != -1) {
        for (int j = 0; j < n_ghost_layers; j++) {
            for (int k = k_start_loc; k <= k_end_loc; k++) {
                #pragma omp simd
                for (int i = i_start_loc; i <= i_end_loc; i++) {
                    arr[I2V(i,loc_J-n_ghost_layers+j,k)] = bjp_r[j*loc_I_excl_ghost*loc_K_excl_ghost + (k-k_start_loc)*loc_I_excl_ghost + (i-i_start_loc)];
                }
            }
        }
    }

    // i-j plane negative
    if (neighbors_id[4] != -1) {
        for (int k = 0; k < n_ghost_layers; k++) {
            for (int j = j_start_loc; j <= j_end_loc; j++) {
                #pragma omp simd
                for (int i = i_start_loc; i <= i_end_loc; i++) {
                    arr[I2V(i,j,k)] = bkn_r[k*loc_I_excl_ghost*loc_J_excl_ghost + (j-j_start_loc)*loc_I_excl_ghost + (i-i_start_loc)];
                }
            }
        }
    }
    // i-j plane positive
    if (neighbors_id[5] != -1) {
        for (int k = 0; k < n_ghost_layers; k++) {
            for (int j = j_start_loc; j <= j_end_loc; j++) {
                #pragma omp simd
                for (int i = i_start_loc; i <= i_end_loc; i++) {
                    arr[I2V(i,j,loc_K-n_ghost_layers+k)] = bkp_r[k*loc_I_excl_ghost*loc_J_excl_ghost + (j-j_start_loc)*loc_I_excl_ghost + (i-i_start_loc)];
                }
            }
        }
    }

    //stdout_by_main("Preparation ghost layers done.");
}

void Grid::assign_received_data_to_ghost_av(CUSTOMREAL* arr) {
    // store the pointers to the ghost layer's elements to be sent / to receive
    // node order should be always smaller to larger

    // j-k plane negative
    if (neighbors_id[0] != -1) {
        for (int i = 0; i < n_ghost_layers; i++) {
            for (int k = k_start_loc; k <= k_end_loc; k++) {
                for (int j = j_start_loc; j <= j_end_loc; j++) {
                    arr[I2V(i,j,k)] += bin_r[i*loc_J_excl_ghost*loc_K_excl_ghost + (k-k_start_loc)*loc_J_excl_ghost + (j-j_start_loc)]/_2_CR;
                }
            }
        }
    }
    // j-k plane positive
    if (neighbors_id[1] != -1) {
        for (int i = 0; i < n_ghost_layers; i++) {
            for (int k = k_start_loc; k <= k_end_loc; k++) {
                for (int j = j_start_loc; j <= j_end_loc; j++) {
                     arr[I2V(loc_I-n_ghost_layers+i,j,k)] += bip_r[i*loc_J_excl_ghost*loc_K_excl_ghost + (k-k_start_loc)*loc_J_excl_ghost + (j-j_start_loc)]/_2_CR;
                }
            }
        }
    }
    // i-k plane negative
    if (neighbors_id[2] != -1) {
        for (int j = 0; j < n_ghost_layers; j++) {
            for (int k = k_start_loc; k <= k_end_loc; k++) {
                for (int i = i_start_loc; i <= i_end_loc; i++) {
                    arr[I2V(i,j,k)] += bjn_r[j*loc_I_excl_ghost*loc_K_excl_ghost + (k-k_start_loc)*loc_I_excl_ghost + (i-i_start_loc)]/_2_CR;
                }
            }
        }
    }
    // i-k plane positive
    if (neighbors_id[3] != -1) {
        for (int j = 0; j < n_ghost_layers; j++) {
            for (int k = k_start_loc; k <= k_end_loc; k++) {
                for (int i = i_start_loc; i <= i_end_loc; i++) {
                    arr[I2V(i,loc_J-n_ghost_layers+j,k)] += bjp_r[j*loc_I_excl_ghost*loc_K_excl_ghost + (k-k_start_loc)*loc_I_excl_ghost + (i-i_start_loc)]/_2_CR;
                }
            }
        }
    }

    // i-j plane negative
    if (neighbors_id[4] != -1) {
        for (int k = 0; k < n_ghost_layers; k++) {
            for (int j = j_start_loc; j <= j_end_loc; j++) {
                for (int i = i_start_loc; i <= i_end_loc; i++) {
                    arr[I2V(i,j,k)] += bkn_r[k*loc_I_excl_ghost*loc_J_excl_ghost + (j-j_start_loc)*loc_I_excl_ghost + (i-i_start_loc)]/_2_CR;
                }
            }
        }
    }
    // i-j plane positive
    if (neighbors_id[5] != -1) {
        for (int k = 0; k < n_ghost_layers; k++) {
            for (int j = j_start_loc; j <= j_end_loc; j++) {
                for (int i = i_start_loc; i <= i_end_loc; i++) {
                    arr[I2V(i,j,loc_K-n_ghost_layers+k)] += bkp_r[k*loc_I_excl_ghost*loc_J_excl_ghost + (j-j_start_loc)*loc_I_excl_ghost + (i-i_start_loc)]/_2_CR;
                }
            }
        }
    }

    //stdout_by_main("Preparation ghost layers done.");
}



void Grid::send_recev_boundary_data(CUSTOMREAL* arr){
    prepare_boundary_data_to_send(arr);

    // i-direction negative
    if (neighbors_id[0] != -1) {
        // send boundary layer to neighbor
        isend_cr(bin_s, n_grid_bound_i*n_ghost_layers, neighbors_id[0], mpi_reqs[0]);
        // receive boundary layer from neighbor
        irecv_cr(bin_r, n_grid_bound_i*n_ghost_layers, neighbors_id[0], mpi_reqs[0]);
    }
    // i-direction positive
    if (neighbors_id[1] != -1) {
        // send boundary layer to neighbor
        isend_cr(bip_s, n_grid_bound_i*n_ghost_layers, neighbors_id[1], mpi_reqs[1]);
        // receive boundary layer from neighbor
        irecv_cr(bip_r, n_grid_bound_i*n_ghost_layers, neighbors_id[1], mpi_reqs[1]);
    }
    // j-direction negative
    if (neighbors_id[2] != -1) {
        // send boundary layer to neighbor
        isend_cr(bjn_s, n_grid_bound_j*n_ghost_layers, neighbors_id[2], mpi_reqs[2]);
        // receive boundary layer from neighbor
        irecv_cr(bjn_r, n_grid_bound_j*n_ghost_layers, neighbors_id[2], mpi_reqs[2]);
    }
    // j-direction positive
    if (neighbors_id[3] != -1) {
        // send boundary layer to neighbor
        isend_cr(bjp_s, n_grid_bound_j*n_ghost_layers, neighbors_id[3], mpi_reqs[3]);
        // receive boundary layer from neighbor
        irecv_cr(bjp_r, n_grid_bound_j*n_ghost_layers, neighbors_id[3], mpi_reqs[3]);
    }
    // k-direction negative
    if (neighbors_id[4] != -1) {
        // send boundary layer to neighbor
        isend_cr(bkn_s, n_grid_bound_k*n_ghost_layers, neighbors_id[4], mpi_reqs[4]);
        // receive boundary layer from neighbor
        irecv_cr(bkn_r, n_grid_bound_k*n_ghost_layers, neighbors_id[4], mpi_reqs[4]);
    }
    // k-direction positive
    if (neighbors_id[5] != -1) {
        // send boundary layer to neighbor
        isend_cr(bkp_s, n_grid_bound_k*n_ghost_layers, neighbors_id[5], mpi_reqs[5]);
        // receive boundary layer from neighbor
        irecv_cr(bkp_r, n_grid_bound_k*n_ghost_layers, neighbors_id[5], mpi_reqs[5]);
    }

    // wait for finishing communication
    for (int i = 0; i < 6; i++) {
        if (neighbors_id[i] != -1) {
            wait_req(mpi_reqs[i]);
        }
    }

    // THIS SYNCHRONIZATION IS IMPORTANT.
    synchronize_all_inter();

    assign_received_data_to_ghost(arr);

}


void Grid::send_recev_boundary_data_av(CUSTOMREAL* arr){
    prepare_boundary_data_to_send(arr);

    // i-direction negative
    if (neighbors_id[0] != -1) {
        // send boundary layer to neighbor
        isend_cr(bin_s, n_grid_bound_i*n_ghost_layers, neighbors_id[0], mpi_reqs[0]);
        // receive boundary layer from neighbor
        irecv_cr(bin_r, n_grid_bound_i*n_ghost_layers, neighbors_id[0], mpi_reqs[0]);
    }
    // i-direction positive
    if (neighbors_id[1] != -1) {
        // send boundary layer to neighbor
        isend_cr(bip_s, n_grid_bound_i*n_ghost_layers, neighbors_id[1], mpi_reqs[1]);
        // receive boundary layer from neighbor
        irecv_cr(bip_r, n_grid_bound_i*n_ghost_layers, neighbors_id[1], mpi_reqs[1]);
    }
    // j-direction negative
    if (neighbors_id[2] != -1) {
        // send boundary layer to neighbor
        isend_cr(bjn_s, n_grid_bound_j*n_ghost_layers, neighbors_id[2], mpi_reqs[2]);
        // receive boundary layer from neighbor
        irecv_cr(bjn_r, n_grid_bound_j*n_ghost_layers, neighbors_id[2], mpi_reqs[2]);
    }
    // j-direction positive
    if (neighbors_id[3] != -1) {
        // send boundary layer to neighbor
        isend_cr(bjp_s, n_grid_bound_j*n_ghost_layers, neighbors_id[3], mpi_reqs[3]);
        // receive boundary layer from neighbor
        irecv_cr(bjp_r, n_grid_bound_j*n_ghost_layers, neighbors_id[3], mpi_reqs[3]);
    }
    // k-direction negative
    if (neighbors_id[4] != -1) {
        // send boundary layer to neighbor
        isend_cr(bkn_s, n_grid_bound_k*n_ghost_layers, neighbors_id[4], mpi_reqs[4]);
        // receive boundary layer from neighbor
        irecv_cr(bkn_r, n_grid_bound_k*n_ghost_layers, neighbors_id[4], mpi_reqs[4]);
    }
    // k-direction positive
    if (neighbors_id[5] != -1) {
        // send boundary layer to neighbor
        isend_cr(bkp_s, n_grid_bound_k*n_ghost_layers, neighbors_id[5], mpi_reqs[5]);
        // receive boundary layer from neighbor
        irecv_cr(bkp_r, n_grid_bound_k*n_ghost_layers, neighbors_id[5], mpi_reqs[5]);
    }

    // wait for finishing communication
    for (int i = 0; i < 6; i++) {
        if (neighbors_id[i] != -1) {
            wait_req(mpi_reqs[i]);
        }
    }

    assign_received_data_to_ghost_av(arr);

}


void Grid::prepare_boundary_data_to_send_kosumi(CUSTOMREAL* arr) {
    // ij nn
    if (neighbors_id_ij[0] != -1)
        for (int k = k_start_loc; k <= k_end_vis; k++)
            bij_nn_s[k-k_start_loc] = arr[I2V(i_start_loc,j_start_loc,k)];
    // ij np
    if (neighbors_id_ij[1] != -1)
        for (int k = k_start_loc; k <= k_end_vis; k++)
            bij_np_s[k-k_start_loc] = arr[I2V(i_start_loc,j_end_loc,k)];
    // ij pn
    if (neighbors_id_ij[2] != -1)
        for (int k = k_start_loc; k <= k_end_vis; k++)
            bij_pn_s[k-k_start_loc] = arr[I2V(i_end_loc,j_start_loc,k)];
    // ij pp
    if (neighbors_id_ij[3] != -1)
        for (int k = k_start_loc; k <= k_end_vis; k++)
            bij_pp_s[k-k_start_loc] = arr[I2V(i_end_loc,j_end_loc,k)];
    // jk nn
    if (neighbors_id_jk[0] != -1)
        for (int i = i_start_loc; i <= i_end_vis; i++)
            bjk_nn_s[i-i_start_loc] = arr[I2V(i,j_start_loc,k_start_loc)];
    // jk np
    if (neighbors_id_jk[1] != -1)
        for (int i = i_start_loc; i <= i_end_vis; i++)
            bjk_np_s[i-i_start_loc] = arr[I2V(i,j_start_loc,k_end_loc)];
    // jk pn
    if (neighbors_id_jk[2] != -1)
        for (int i = i_start_loc; i <= i_end_vis; i++)
            bjk_pn_s[i-i_start_loc] = arr[I2V(i,j_end_loc,k_start_loc)];
    // jk pp
    if (neighbors_id_jk[3] != -1)
        for (int i = i_start_loc; i <= i_end_vis; i++)
            bjk_pp_s[i-i_start_loc] = arr[I2V(i,j_end_loc,k_end_loc)];
    // ik nn
    if (neighbors_id_ik[0] != -1)
        for (int j = j_start_loc; j <= j_end_vis; j++)
            bik_nn_s[j-j_start_loc] = arr[I2V(i_start_loc,j,k_start_loc)];
    // ik np
    if (neighbors_id_ik[1] != -1)
        for (int j = j_start_loc; j <= j_end_vis; j++)
            bik_np_s[j-j_start_loc] = arr[I2V(i_start_loc,j,k_end_loc)];
    // ik pn
    if (neighbors_id_ik[2] != -1)
        for (int j = j_start_loc; j <= j_end_vis; j++)
            bik_pn_s[j-j_start_loc] = arr[I2V(i_end_loc,j,k_start_loc)];
    // ik pp
    if (neighbors_id_ik[3] != -1)
        for (int j = j_start_loc; j <= j_end_vis; j++)
            bik_pp_s[j-j_start_loc] = arr[I2V(i_end_loc,j,k_end_loc)];
    // ijk nnn
    if (neighbors_id_ijk[0] != -1)
        bijk_nnn_s[0] = arr[I2V(i_start_loc,j_start_loc,k_start_loc)];
    // ijk nnp
    if (neighbors_id_ijk[1] != -1)
        bijk_nnp_s[0] = arr[I2V(i_start_loc,j_start_loc,k_end_loc)];
    // ijk npn
    if (neighbors_id_ijk[2] != -1)
        bijk_npn_s[0] = arr[I2V(i_start_loc,j_end_loc,k_start_loc)];
    // ijk npp
    if (neighbors_id_ijk[3] != -1)
        bijk_npp_s[0] = arr[I2V(i_start_loc,j_end_loc,k_end_loc)];
    // ijk pnn
    if (neighbors_id_ijk[4] != -1)
        bijk_pnn_s[0] = arr[I2V(i_end_loc,j_start_loc,k_start_loc)];
    // ijk pnp
    if (neighbors_id_ijk[5] != -1)
        bijk_pnp_s[0] = arr[I2V(i_end_loc,j_start_loc,k_end_loc)];
    // ijk ppn
    if (neighbors_id_ijk[6] != -1)
        bijk_ppn_s[0] = arr[I2V(i_end_loc,j_end_loc,k_start_loc)];
    // ijk ppp
    if (neighbors_id_ijk[7] != -1)
        bijk_ppp_s[0] = arr[I2V(i_end_loc,j_end_loc,k_end_loc)];
}


void Grid::assign_received_data_to_ghost_kosumi(CUSTOMREAL* arr) {
    // ij nn
    if (neighbors_id_ij[0] != -1)
        for (int k = k_start_loc; k <= k_end_vis; k++)
            arr[I2V(i_start_loc-1,j_start_loc-1,k)] = bij_nn_r[k-k_start_loc];
    // ij np
    if (neighbors_id_ij[1] != -1)
        for (int k = k_start_loc; k <= k_end_vis; k++)
            arr[I2V(i_start_loc-1,j_end_vis,k)] = bij_np_r[k-k_start_loc];
    // ij pn
    if (neighbors_id_ij[2] != -1)
        for (int k = k_start_loc; k <= k_end_vis; k++)
            arr[I2V(i_end_vis,j_start_loc-1,k)] = bij_pn_r[k-k_start_loc];
    // ij pp
    if (neighbors_id_ij[3] != -1)
        for (int k = k_start_loc; k <= k_end_vis; k++)
            arr[I2V(i_end_vis,j_end_loc+1,k)] = bij_pp_r[k-k_start_loc];
    // jk nn
    if (neighbors_id_jk[0] != -1)
        for (int i = i_start_loc; i <= i_end_vis; i++)
            arr[I2V(i,j_start_loc-1,k_start_loc-1)] = bjk_nn_r[i-i_start_loc];
    // jk np
    if (neighbors_id_jk[1] != -1)
        for (int i = i_start_loc; i <= i_end_vis; i++)
            arr[I2V(i,j_start_loc-1,k_end_vis)] = bjk_np_r[i-i_start_loc];
    // jk pn
    if (neighbors_id_jk[2] != -1)
        for (int i = i_start_loc; i <= i_end_vis; i++)
            arr[I2V(i,j_end_vis,k_start_loc-1)] = bjk_pn_r[i-i_start_loc];
    // jk pp
    if (neighbors_id_jk[3] != -1)
        for (int i = i_start_loc; i <= i_end_vis; i++)
            arr[I2V(i,j_end_vis,k_end_vis)] = bjk_pp_r[i-i_start_loc];
    // ik nn
    if (neighbors_id_ik[0] != -1)
        for (int j = j_start_loc; j <= j_end_vis; j++)
            arr[I2V(i_start_loc-1,j,k_start_loc-1)] = bik_nn_r[j-j_start_loc];
    // ik np
    if (neighbors_id_ik[1] != -1)
        for (int j = j_start_loc; j <= j_end_vis; j++)
            arr[I2V(i_start_loc-1,j,k_end_vis)] = bik_np_r[j-j_start_loc];
    // ik pn
    if (neighbors_id_ik[2] != -1)
        for (int j = j_start_loc; j <= j_end_vis; j++)
            arr[I2V(i_end_vis,j,k_start_loc-1)] = bik_pn_r[j-j_start_loc];
    // ik pp
    if (neighbors_id_ik[3] != -1)
        for (int j = j_start_loc; j <= j_end_vis; j++)
            arr[I2V(i_end_vis,j,k_end_vis)] = bik_pp_r[j-j_start_loc];
    // ijk nnn
    if (neighbors_id_ijk[0] != -1)
        arr[I2V(i_start_loc-1,j_start_loc-1,k_start_loc-1)] = bijk_nnn_r[0];
    // ijk nnp
    if (neighbors_id_ijk[1] != -1)
        arr[I2V(i_start_loc-1,j_start_loc-1,k_end_vis)] = bijk_nnp_r[0];
    // ijk npn
    if (neighbors_id_ijk[2] != -1)
        arr[I2V(i_start_loc-1,j_end_vis,k_start_loc-1)] = bijk_npn_r[0];
    // ijk npp
    if (neighbors_id_ijk[3] != -1)
        arr[I2V(i_start_loc-1,j_end_vis,k_end_vis)] = bijk_npp_r[0];
    // ijk pnn
    if (neighbors_id_ijk[4] != -1)
        arr[I2V(i_end_vis,j_start_loc-1,k_start_loc-1)] = bijk_pnn_r[0];
    // ijk pnp
    if (neighbors_id_ijk[5] != -1)
        arr[I2V(i_end_vis,j_start_loc-1,k_end_vis)] = bijk_pnp_r[0];
    // ijk ppn
    if (neighbors_id_ijk[6] != -1)
        arr[I2V(i_end_vis,j_end_vis,k_start_loc-1)] = bijk_ppn_r[0];
    // ijk ppp
    if (neighbors_id_ijk[7] != -1)
        arr[I2V(i_end_vis,j_end_vis,k_end_vis)] = bijk_ppp_r[0];

}


void Grid::send_recev_boundary_data_kosumi(CUSTOMREAL* arr){

    prepare_boundary_data_to_send_kosumi(arr);

    // ij nn
    if (neighbors_id_ij[0] != -1){
        // send boundary layer to neighbor
        isend_cr(bij_nn_s, loc_K_vis, neighbors_id_ij[0], mpi_reqs_kosumi[0]);
        // receive boundary layer from neighbor
        irecv_cr(bij_nn_r, loc_K_vis, neighbors_id_ij[0], mpi_reqs_kosumi[0]);
    }
    // ij np
    if (neighbors_id_ij[1] != -1){
        // send boundary layer to neighbor
        isend_cr(bij_np_s, loc_K_vis, neighbors_id_ij[1], mpi_reqs_kosumi[1]);
        // receive boundary layer from neighbor
        irecv_cr(bij_np_r, loc_K_vis, neighbors_id_ij[1], mpi_reqs_kosumi[1]);
    }
    // ij pn
    if (neighbors_id_ij[2] != -1){
        // send boundary layer to neighbor
        isend_cr(bij_pn_s, loc_K_vis, neighbors_id_ij[2], mpi_reqs_kosumi[2]);
        // receive boundary layer from neighbor
        irecv_cr(bij_pn_r, loc_K_vis, neighbors_id_ij[2], mpi_reqs_kosumi[2]);
    }
    // ij pp
    if (neighbors_id_ij[3] != -1){
        // send boundary layer to neighbor
        isend_cr(bij_pp_s, loc_K_vis, neighbors_id_ij[3], mpi_reqs_kosumi[3]);
        // receive boundary layer from neighbor
        irecv_cr(bij_pp_r, loc_K_vis, neighbors_id_ij[3], mpi_reqs_kosumi[3]);
    }
    // jk nn
    if (neighbors_id_jk[0] != -1){
        // send boundary layer to neighbor
        isend_cr(bjk_nn_s, loc_I_vis, neighbors_id_jk[0], mpi_reqs_kosumi[4]);
        // receive boundary layer from neighbor
        irecv_cr(bjk_nn_r, loc_I_vis, neighbors_id_jk[0], mpi_reqs_kosumi[4]);
    }
    // jk np
    if (neighbors_id_jk[1] != -1){
        // send boundary layer to neighbor
        isend_cr(bjk_np_s, loc_I_vis, neighbors_id_jk[1], mpi_reqs_kosumi[5]);
        // receive boundary layer from neighbor
        irecv_cr(bjk_np_r, loc_I_vis, neighbors_id_jk[1], mpi_reqs_kosumi[5]);
    }
    // jk pn
    if (neighbors_id_jk[2] != -1){
        // send boundary layer to neighbor
        isend_cr(bjk_pn_s, loc_I_vis, neighbors_id_jk[2], mpi_reqs_kosumi[6]);
        // receive boundary layer from neighbor
        irecv_cr(bjk_pn_r, loc_I_vis, neighbors_id_jk[2], mpi_reqs_kosumi[6]);
    }
    // jk pp
    if (neighbors_id_jk[3] != -1){
        // send boundary layer to neighbor
        isend_cr(bjk_pp_s, loc_I_vis, neighbors_id_jk[3], mpi_reqs_kosumi[7]);
        // receive boundary layer from neighbor
        irecv_cr(bjk_pp_r, loc_I_vis, neighbors_id_jk[3], mpi_reqs_kosumi[7]);
    }
    // ik nn
    if (neighbors_id_ik[0] != -1){
        // send boundary layer to neighbor
        isend_cr(bik_nn_s, loc_J_vis, neighbors_id_ik[0], mpi_reqs_kosumi[8]);
        // receive boundary layer from neighbor
        irecv_cr(bik_nn_r, loc_J_vis, neighbors_id_ik[0], mpi_reqs_kosumi[8]);
    }
    // ik np
    if (neighbors_id_ik[1] != -1){
        // send boundary layer to neighbor
        isend_cr(bik_np_s, loc_J_vis, neighbors_id_ik[1], mpi_reqs_kosumi[9]);
        // receive boundary layer from neighbor
        irecv_cr(bik_np_r, loc_J_vis, neighbors_id_ik[1], mpi_reqs_kosumi[9]);
    }
    // ik pn
    if (neighbors_id_ik[2] != -1){
        // send boundary layer to neighbor
        isend_cr(bik_pn_s, loc_J_vis, neighbors_id_ik[2], mpi_reqs_kosumi[10]);
        // receive boundary layer from neighbor
        irecv_cr(bik_pn_r, loc_J_vis, neighbors_id_ik[2], mpi_reqs_kosumi[10]);
    }
    // ik pp
    if (neighbors_id_ik[3] != -1){
        // send boundary layer to neighbor
        isend_cr(bik_pp_s, loc_J_vis, neighbors_id_ik[3], mpi_reqs_kosumi[11]);
        // receive boundary layer from neighbor
        irecv_cr(bik_pp_r, loc_J_vis, neighbors_id_ik[3], mpi_reqs_kosumi[11]);
    }
    // ijk nnn
    if (neighbors_id_ijk[0] != -1){
        // send boundary layer to neighbor
        isend_cr(bijk_nnn_s, 1, neighbors_id_ijk[0], mpi_reqs_kosumi[12]);
        // receive boundary layer from neighbor
        irecv_cr(bijk_nnn_r, 1, neighbors_id_ijk[0], mpi_reqs_kosumi[12]);
    }
    // ijk nnp
    if (neighbors_id_ijk[1] != -1){
        // send boundary layer to neighbor
        isend_cr(bijk_nnp_s, 1, neighbors_id_ijk[1], mpi_reqs_kosumi[13]);
        // receive boundary layer from neighbor
        irecv_cr(bijk_nnp_r, 1, neighbors_id_ijk[1], mpi_reqs_kosumi[13]);
    }
    // ijk npn
    if (neighbors_id_ijk[2] != -1){
        // send boundary layer to neighbor
        isend_cr(bijk_npn_s, 1, neighbors_id_ijk[2], mpi_reqs_kosumi[14]);
        // receive boundary layer from neighbor
        irecv_cr(bijk_npn_r, 1, neighbors_id_ijk[2], mpi_reqs_kosumi[14]);
    }
    // ijk npp
    if (neighbors_id_ijk[3] != -1){
        // send boundary layer to neighbor
        isend_cr(bijk_npp_s, 1, neighbors_id_ijk[3], mpi_reqs_kosumi[15]);
        // receive boundary layer from neighbor
        irecv_cr(bijk_npp_r, 1, neighbors_id_ijk[3], mpi_reqs_kosumi[15]);
    }
    // ijk pnn
    if (neighbors_id_ijk[4] != -1){
        // send boundary layer to neighbor
        isend_cr(bijk_pnn_s, 1, neighbors_id_ijk[4], mpi_reqs_kosumi[16]);
        // receive boundary layer from neighbor
        irecv_cr(bijk_pnn_r, 1, neighbors_id_ijk[4], mpi_reqs_kosumi[16]);
    }
    // ijk pnp
    if (neighbors_id_ijk[5] != -1){
        // send boundary layer to neighbor
        isend_cr(bijk_pnp_s, 1, neighbors_id_ijk[5], mpi_reqs_kosumi[17]);
        // receive boundary layer from neighbor
        irecv_cr(bijk_pnp_r, 1, neighbors_id_ijk[5], mpi_reqs_kosumi[17]);
    }
    // ijk ppn
    if (neighbors_id_ijk[6] != -1){
        // send boundary layer to neighbor
        isend_cr(bijk_ppn_s, 1, neighbors_id_ijk[6], mpi_reqs_kosumi[18]);
        // receive boundary layer from neighbor
        irecv_cr(bijk_ppn_r, 1, neighbors_id_ijk[6], mpi_reqs_kosumi[18]);
    }
    // ijk ppp
    if (neighbors_id_ijk[7] != -1){
        // send boundary layer to neighbor
        isend_cr(bijk_ppp_s, 1, neighbors_id_ijk[7], mpi_reqs_kosumi[19]);
        // receive boundary layer from neighbor
        irecv_cr(bijk_ppp_r, 1, neighbors_id_ijk[7], mpi_reqs_kosumi[19]);
    }

    // wait for all communication to finish
    // ij nn
    if (neighbors_id_ij[0] != -1)
        wait_req(mpi_reqs_kosumi[0]);
    // ij np
    if (neighbors_id_ij[1] != -1)
        wait_req(mpi_reqs_kosumi[1]);
    // ij pn
    if (neighbors_id_ij[2] != -1)
        wait_req(mpi_reqs_kosumi[2]);
    // ij pp
    if (neighbors_id_ij[3] != -1)
        wait_req(mpi_reqs_kosumi[3]);
    // jk nn
    if (neighbors_id_jk[0] != -1)
        wait_req(mpi_reqs_kosumi[4]);
    // jk np
    if (neighbors_id_jk[1] != -1)
        wait_req(mpi_reqs_kosumi[5]);
    // jk pn
    if (neighbors_id_jk[2] != -1)
        wait_req(mpi_reqs_kosumi[6]);
    // jk pp
    if (neighbors_id_jk[3] != -1)
        wait_req(mpi_reqs_kosumi[7]);
    // ik nn
    if (neighbors_id_ik[0] != -1)
        wait_req(mpi_reqs_kosumi[8]);
    // ik np
    if (neighbors_id_ik[1] != -1)
        wait_req(mpi_reqs_kosumi[9]);
    // ik pn
    if (neighbors_id_ik[2] != -1)
        wait_req(mpi_reqs_kosumi[10]);
    // ik pp
    if (neighbors_id_ik[3] != -1)
        wait_req(mpi_reqs_kosumi[11]);
    // ijk nnn
    if (neighbors_id_ijk[0] != -1)
        wait_req(mpi_reqs_kosumi[12]);
    // ijk nnp
    if (neighbors_id_ijk[1] != -1)
        wait_req(mpi_reqs_kosumi[13]);
    // ijk npn
    if (neighbors_id_ijk[2] != -1)
        wait_req(mpi_reqs_kosumi[14]);
    // ijk npp
    if (neighbors_id_ijk[3] != -1)
        wait_req(mpi_reqs_kosumi[15]);
    // ijk pnn
    if (neighbors_id_ijk[4] != -1)
        wait_req(mpi_reqs_kosumi[16]);
    // ijk pnp
    if (neighbors_id_ijk[5] != -1)
        wait_req(mpi_reqs_kosumi[17]);
    // ijk ppn
    if (neighbors_id_ijk[6] != -1)
        wait_req(mpi_reqs_kosumi[18]);
    // ijk ppp
    if (neighbors_id_ijk[7] != -1)
        wait_req(mpi_reqs_kosumi[19]);

    assign_received_data_to_ghost_kosumi(arr);

}


void Grid::calc_T_plus_tau() {
    // calculate T_plus_tau
    for (int k_r = 0; k_r < loc_K; k_r++) {
        for (int j_lat = 0; j_lat < loc_J; j_lat++) {
            #pragma omp simd
            for (int i_lon = 0; i_lon < loc_I; i_lon++) {
                // T0*tau
                T_loc[I2V(i_lon,j_lat,k_r)] = T0v_loc[I2V(i_lon,j_lat,k_r)] * tau_loc[I2V(i_lon,j_lat,k_r)];

            }
        }
    }
}


void Grid::calc_residual() {
    for (int k_r = 0; k_r < loc_K; k_r++) {
        for (int j_lat = 0; j_lat < loc_J; j_lat++) {
            #pragma omp simd
            for (int i_lon = 0; i_lon < loc_I; i_lon++) {
                u_loc[I2V(i_lon,j_lat,k_r)] = u_loc[I2V(i_lon,j_lat,k_r)] - T_loc[I2V(i_lon,j_lat,k_r)];
            }
        }
    }
}

void Grid::write_inversion_grid_file(){

    std::ofstream ofs;

    inversion_grid_file_out = output_dir + "inversion_grid.txt";
    ofs.open(inversion_grid_file_out);

    if(subdom_main && id_subdomain == 0){       // main processor of subdomain && the first id of subdoumains
        for(int l = 0; l < n_inv_grids; l++){
            ofs << l << " " << ngrid_k_inv << " " << ngrid_j_inv << " " << ngrid_i_inv << std::endl;    // number of ivnersion grid
            // inversion grid of depth
            for(int k =0; k < ngrid_k_inv; k++)
                ofs << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ')
                    << radius2depth(r_loc_inv[I2V_INV_GRIDS_1DK(k,l)]) << " ";
            ofs << std::endl;
            for(int j =0; j < ngrid_j_inv; j++)
                ofs << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ')
                    << t_loc_inv[I2V_INV_GRIDS_1DJ(j,l)]*RAD2DEG << " ";
            ofs << std::endl;
            for(int i =0; i < ngrid_i_inv; i++)
                ofs << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ')
                    << p_loc_inv[I2V_INV_GRIDS_1DI(i,l)]*RAD2DEG << " ";
            ofs << std::endl;
        }
    }

}
