#include "grid_wrapper.cuh"

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
                std::vector<std::vector<CUSTOMREAL*>> const& vv_change){

    // store grid parameters
    grid_dv->loc_I_host = loc_I;
    grid_dv->loc_J_host = loc_J;
    grid_dv->loc_K_host = loc_K;

    // count node number
    grid_dv->n_nodes_total_host = loc_I*loc_J*loc_K;
    grid_dv->n_levels_host = ijk.size();
    int n_nodes_per_level[grid_dv->n_levels_host];
    for (int i = 0; i < grid_dv->n_levels_host; i++){
        n_nodes_per_level[i] = ijk.at(i).size();
    }

    // allocate memory on device
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->n_nodes_on_levels), grid_dv->n_levels_host),  0);

    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__k__0), grid_dv->n_nodes_total_host), 1);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_ip1j__k__0), grid_dv->n_nodes_total_host), 2);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_im1j__k__0), grid_dv->n_nodes_total_host), 3);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__jp1k__0), grid_dv->n_nodes_total_host), 4);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__jm1k__0), grid_dv->n_nodes_total_host), 5);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__kp10), grid_dv->n_nodes_total_host), 6);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__km10), grid_dv->n_nodes_total_host), 7);

    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__k__1), grid_dv->n_nodes_total_host), 8);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_ip1j__k__1), grid_dv->n_nodes_total_host), 9);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_im1j__k__1), grid_dv->n_nodes_total_host), 10);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__jp1k__1), grid_dv->n_nodes_total_host), 11);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__jm1k__1), grid_dv->n_nodes_total_host), 12);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__kp11), grid_dv->n_nodes_total_host), 13);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__km11), grid_dv->n_nodes_total_host), 14);

    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__k__2), grid_dv->n_nodes_total_host), 15);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_ip1j__k__2), grid_dv->n_nodes_total_host), 16);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_im1j__k__2), grid_dv->n_nodes_total_host), 17);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__jp1k__2), grid_dv->n_nodes_total_host), 18);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__jm1k__2), grid_dv->n_nodes_total_host), 19);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__kp12), grid_dv->n_nodes_total_host), 20);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__km12), grid_dv->n_nodes_total_host), 21);

    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__k__3), grid_dv->n_nodes_total_host), 22);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_ip1j__k__3), grid_dv->n_nodes_total_host), 23);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_im1j__k__3), grid_dv->n_nodes_total_host), 24);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__jp1k__3), grid_dv->n_nodes_total_host), 25);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__jm1k__3), grid_dv->n_nodes_total_host), 26);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__kp13), grid_dv->n_nodes_total_host), 27);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__km13), grid_dv->n_nodes_total_host), 28);

    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__k__4), grid_dv->n_nodes_total_host), 29);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_ip1j__k__4), grid_dv->n_nodes_total_host), 30);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_im1j__k__4), grid_dv->n_nodes_total_host), 31);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__jp1k__4), grid_dv->n_nodes_total_host), 32);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__jm1k__4), grid_dv->n_nodes_total_host), 33);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__kp14), grid_dv->n_nodes_total_host), 34);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__km14), grid_dv->n_nodes_total_host), 35);

    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__k__5), grid_dv->n_nodes_total_host), 36);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_ip1j__k__5), grid_dv->n_nodes_total_host), 37);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_im1j__k__5), grid_dv->n_nodes_total_host), 38);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__jp1k__5), grid_dv->n_nodes_total_host), 39);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__jm1k__5), grid_dv->n_nodes_total_host), 40);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__kp15), grid_dv->n_nodes_total_host), 41);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__km15), grid_dv->n_nodes_total_host), 42);

    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__k__6), grid_dv->n_nodes_total_host), 43);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_ip1j__k__6), grid_dv->n_nodes_total_host), 44);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_im1j__k__6), grid_dv->n_nodes_total_host), 45);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__jp1k__6), grid_dv->n_nodes_total_host), 46);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__jm1k__6), grid_dv->n_nodes_total_host), 47);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__kp16), grid_dv->n_nodes_total_host), 48);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__km16), grid_dv->n_nodes_total_host), 49);

    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__k__7), grid_dv->n_nodes_total_host), 50);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_ip1j__k__7), grid_dv->n_nodes_total_host), 51);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_im1j__k__7), grid_dv->n_nodes_total_host), 52);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__jp1k__7), grid_dv->n_nodes_total_host), 53);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__jm1k__7), grid_dv->n_nodes_total_host), 54);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__kp17), grid_dv->n_nodes_total_host), 55);
    print_CUDA_error_if_any(allocate_memory_on_device_i((void**) &(grid_dv->vv_i__j__km17), grid_dv->n_nodes_total_host), 56);

    // transfer data from host to device


}

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
                std::vector<std::vector<CUSTOMREAL*>> const& vv_change){



}

