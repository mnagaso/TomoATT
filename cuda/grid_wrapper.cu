#include "grid_wrapper.cuh"

void cuda_initialize_grid_1st(std::vector< std::vector<int> >& ijk, Grid_on_device* grid_dv, int const& loc_I, int const& loc_J, int const& loc_K,
                CUSTOMREAL const& dp, CUSTOMREAL const& dt, CUSTOMREAL const& dr, \
                std::vector<std::vector<int*>>        & vv_i__j__k__, \
                std::vector<std::vector<int*>>        & vv_ip1j__k__, \
                std::vector<std::vector<int*>>        & vv_im1j__k__, \
                std::vector<std::vector<int*>>        & vv_i__jp1k__, \
                std::vector<std::vector<int*>>        & vv_i__jm1k__, \
                std::vector<std::vector<int*>>        & vv_i__j__kp1, \
                std::vector<std::vector<int*>>        & vv_i__j__km1, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_fac_a, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_fac_b, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_fac_c, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_fac_f, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_T0v, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_T0r, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_T0t, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_T0p, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_fun, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_change){

    // store grid parameters
    grid_dv->loc_I_host = loc_I;
    grid_dv->loc_J_host = loc_J;
    grid_dv->loc_K_host = loc_K;
    grid_dv->dr_host = dr;
    grid_dv->dt_host = dt;
    grid_dv->dp_host = dp;

    // count node number
    grid_dv->n_nodes_total_host = 0;
    grid_dv->n_levels_host = ijk.size();
    // allocate grid_dv->n_nodes_on_levels_host
    grid_dv->n_nodes_on_levels_host = new int[grid_dv->n_levels_host];

    for (int i=0; i<grid_dv->n_levels_host; i++){
        grid_dv->n_nodes_on_levels_host[i] = ijk[i].size();
        grid_dv->n_nodes_total_host += grid_dv->n_nodes_on_levels_host[i];
    }

    // allocate memory on device
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_i(grid_dv->n_nodes_on_levels, grid_dv->n_nodes_on_levels_host, grid_dv->n_levels_host), 11);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__k___0, vv_i__j__k__.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 2);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip1j__k___0, vv_ip1j__k__.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 2);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im1j__k___0, vv_im1j__k__.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 3);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp1k___0, vv_i__jp1k__.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 4);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm1k___0, vv_i__jm1k__.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 5);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp1_0, vv_i__j__kp1.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 6);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km1_0, vv_i__j__km1.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 7);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__k___1, vv_i__j__k__.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 8);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip1j__k___1, vv_ip1j__k__.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 9);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im1j__k___1, vv_im1j__k__.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 10);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp1k___1, vv_i__jp1k__.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 11);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm1k___1, vv_i__jm1k__.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 12);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp1_1, vv_i__j__kp1.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 13);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km1_1, vv_i__j__km1.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 14);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__k___2, vv_i__j__k__.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 15);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip1j__k___2, vv_ip1j__k__.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 16);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im1j__k___2, vv_im1j__k__.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 17);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp1k___2, vv_i__jp1k__.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 18);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm1k___2, vv_i__jm1k__.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 19);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp1_2, vv_i__j__kp1.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 20);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km1_2, vv_i__j__km1.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 21);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__k___3, vv_i__j__k__.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 22);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip1j__k___3, vv_ip1j__k__.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 23);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im1j__k___3, vv_im1j__k__.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 24);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp1k___3, vv_i__jp1k__.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 25);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm1k___3, vv_i__jm1k__.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 26);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp1_3, vv_i__j__kp1.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 27);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km1_3, vv_i__j__km1.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 28);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__k___4, vv_i__j__k__.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 29);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip1j__k___4, vv_ip1j__k__.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 30);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im1j__k___4, vv_im1j__k__.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 31);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp1k___4, vv_i__jp1k__.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 32);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm1k___4, vv_i__jm1k__.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 33);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp1_4, vv_i__j__kp1.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 34);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km1_4, vv_i__j__km1.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 35);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__k___5, vv_i__j__k__.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 36);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip1j__k___5, vv_ip1j__k__.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 37);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im1j__k___5, vv_im1j__k__.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 38);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp1k___5, vv_i__jp1k__.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 39);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm1k___5, vv_i__jm1k__.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 40);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp1_5, vv_i__j__kp1.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 41);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km1_5, vv_i__j__km1.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 42);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__k___6, vv_i__j__k__.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 43);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip1j__k___6, vv_ip1j__k__.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 44);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im1j__k___6, vv_im1j__k__.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 45);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp1k___6, vv_i__jp1k__.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 46);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm1k___6, vv_i__jm1k__.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 47);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp1_6, vv_i__j__kp1.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 48);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km1_6, vv_i__j__km1.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 49);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__k___7, vv_i__j__k__.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 50);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip1j__k___7, vv_ip1j__k__.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 51);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im1j__k___7, vv_im1j__k__.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 52);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp1k___7, vv_i__jp1k__.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 53);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm1k___7, vv_i__jm1k__.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 54);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp1_7, vv_i__j__kp1.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 55);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km1_7, vv_i__j__km1.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 56);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_a_0,  vv_fac_a.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 57);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_b_0,  vv_fac_b.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 58);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_c_0,  vv_fac_c.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 59);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_f_0,  vv_fac_f.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 60);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0v_0,    vv_T0v.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 61);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0r_0,    vv_T0r.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 62);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0t_0,    vv_T0t.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 63);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0p_0,    vv_T0p.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 64);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_fun_0,    vv_fun.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 65);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(grid_dv->vv_change_0, vv_change.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 66);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_a_1,  vv_fac_a.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 67);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_b_1,  vv_fac_b.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 68);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_c_1,  vv_fac_c.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 69);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_f_1,  vv_fac_f.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 70);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0v_1,    vv_T0v.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 71);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0r_1,    vv_T0r.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 72);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0t_1,    vv_T0t.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 73);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0p_1,    vv_T0p.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 74);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_fun_1,    vv_fun.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 75);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(grid_dv->vv_change_1, vv_change.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 76);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_a_2,  vv_fac_a.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 77);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_b_2,  vv_fac_b.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 78);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_c_2,  vv_fac_c.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 79);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_f_2,  vv_fac_f.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 80);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0v_2,    vv_T0v.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 81);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0r_2,    vv_T0r.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 82);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0t_2,    vv_T0t.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 83);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0p_2,    vv_T0p.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 84);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_fun_2,    vv_fun.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 85);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(grid_dv->vv_change_2, vv_change.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 86);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_a_3,  vv_fac_a.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 77);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_b_3,  vv_fac_b.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 78);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_c_3,  vv_fac_c.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 79);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_f_3,  vv_fac_f.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 80);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0v_3,    vv_T0v.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 81);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0r_3,    vv_T0r.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 82);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0t_3,    vv_T0t.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 83);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0p_3,    vv_T0p.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 84);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_fun_3,    vv_fun.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 85);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(grid_dv->vv_change_3, vv_change.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 86);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_a_4,  vv_fac_a.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 77);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_b_4,  vv_fac_b.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 78);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_c_4,  vv_fac_c.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 79);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_f_4,  vv_fac_f.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 80);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0v_4,    vv_T0v.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 81);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0r_4,    vv_T0r.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 82);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0t_4,    vv_T0t.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 83);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0p_4,    vv_T0p.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 84);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_fun_4,    vv_fun.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 85);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(grid_dv->vv_change_4, vv_change.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 86);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_a_5,  vv_fac_a.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 77);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_b_5,  vv_fac_b.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 78);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_c_5,  vv_fac_c.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 79);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_f_5,  vv_fac_f.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 80);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0v_5,    vv_T0v.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 81);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0r_5,    vv_T0r.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 82);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0t_5,    vv_T0t.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 83);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0p_5,    vv_T0p.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 84);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_fun_5,    vv_fun.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 85);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(grid_dv->vv_change_5, vv_change.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 86);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_a_6,  vv_fac_a.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 77);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_b_6,  vv_fac_b.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 78);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_c_6,  vv_fac_c.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 79);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_f_6,  vv_fac_f.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 80);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0v_6,    vv_T0v.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 81);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0r_6,    vv_T0r.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 82);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0t_6,    vv_T0t.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 83);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0p_6,    vv_T0p.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 84);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_fun_6,    vv_fun.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 85);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(grid_dv->vv_change_6, vv_change.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 86);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_a_7,  vv_fac_a.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 77);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_b_7,  vv_fac_b.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 78);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_c_7,  vv_fac_c.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 79);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_f_7,  vv_fac_f.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 80);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0v_7,    vv_T0v.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 81);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0r_7,    vv_T0r.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 82);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0t_7,    vv_T0t.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 83);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0p_7,    vv_T0p.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 84);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_fun_7,    vv_fun.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 85);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(grid_dv->vv_change_7, vv_change.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 86);

    // allocate tau
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**)&(grid_dv->tau), grid_dv->n_nodes_total_host), 87);


}

void cuda_initialize_grid_3rd(std::vector< std::vector<int> >& ijk, Grid_on_device* grid_dv, int const& loc_I, int const& loc_J, int const& loc_K,
                CUSTOMREAL const& dp, CUSTOMREAL const& dt, CUSTOMREAL const& dr, \
                std::vector<std::vector<int*>>        & vv_i__j__k__, \
                std::vector<std::vector<int*>>        & vv_ip1j__k__, \
                std::vector<std::vector<int*>>        & vv_im1j__k__, \
                std::vector<std::vector<int*>>        & vv_i__jp1k__, \
                std::vector<std::vector<int*>>        & vv_i__jm1k__, \
                std::vector<std::vector<int*>>        & vv_i__j__kp1, \
                std::vector<std::vector<int*>>        & vv_i__j__km1, \
                std::vector<std::vector<int*>>        & vv_ip2j__k__, \
                std::vector<std::vector<int*>>        & vv_im2j__k__, \
                std::vector<std::vector<int*>>        & vv_i__jp2k__, \
                std::vector<std::vector<int*>>        & vv_i__jm2k__, \
                std::vector<std::vector<int*>>        & vv_i__j__kp2, \
                std::vector<std::vector<int*>>        & vv_i__j__km2, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_fac_a, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_fac_b, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_fac_c, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_fac_f, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_T0v, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_T0r, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_T0t, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_T0p, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_fun, \
                std::vector<std::vector<CUSTOMREAL*>> & vv_change){

    grid_dv->if_3rd_order = true;

    // store grid parameters
    grid_dv->loc_I_host = loc_I;
    grid_dv->loc_J_host = loc_J;
    grid_dv->loc_K_host = loc_K;

    // count node number
    grid_dv->n_nodes_total_host = 0;
    grid_dv->n_levels_host = ijk.size();
    grid_dv->n_nodes_on_levels_host = new int[grid_dv->n_levels_host];

    for (int i = 0; i < grid_dv->n_levels_host; i++){
        grid_dv->n_nodes_on_levels_host[i] = ijk.at(i).size();
        grid_dv->n_nodes_total_host += grid_dv->n_nodes_on_levels_host[i];
    }

    // allocate memory on device
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_i(grid_dv->n_nodes_on_levels, grid_dv->n_nodes_on_levels_host, grid_dv->n_levels_host), 11);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__k___0, vv_i__j__k__.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 1001);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip1j__k___0, vv_ip1j__k__.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 2);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im1j__k___0, vv_im1j__k__.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 3);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp1k___0, vv_i__jp1k__.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 4);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm1k___0, vv_i__jm1k__.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 5);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp1_0, vv_i__j__kp1.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 6);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km1_0, vv_i__j__km1.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 7);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip2j__k___0, vv_ip2j__k__.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 8);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im2j__k___0, vv_im2j__k__.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 9);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp2k___0, vv_i__jp2k__.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 10);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm2k___0, vv_i__jm2k__.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 11);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp2_0, vv_i__j__kp2.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 12);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km2_0, vv_i__j__km2.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 13);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__k___1, vv_i__j__k__.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 8);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip1j__k___1, vv_ip1j__k__.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 9);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im1j__k___1, vv_im1j__k__.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 10);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp1k___1, vv_i__jp1k__.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 11);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm1k___1, vv_i__jm1k__.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 12);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp1_1, vv_i__j__kp1.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 13);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km1_1, vv_i__j__km1.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 14);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip2j__k___1, vv_ip2j__k__.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 8);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im2j__k___1, vv_im2j__k__.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 9);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp2k___1, vv_i__jp2k__.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 10);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm2k___1, vv_i__jm2k__.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 11);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp2_1, vv_i__j__kp2.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 12);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km2_1, vv_i__j__km2.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 13);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__k___2, vv_i__j__k__.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 15);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip1j__k___2, vv_ip1j__k__.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 16);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im1j__k___2, vv_im1j__k__.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 17);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp1k___2, vv_i__jp1k__.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 18);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm1k___2, vv_i__jm1k__.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 19);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp1_2, vv_i__j__kp1.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 20);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km1_2, vv_i__j__km1.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 21);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip2j__k___2, vv_ip2j__k__.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 8);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im2j__k___2, vv_im2j__k__.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 9);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp2k___2, vv_i__jp2k__.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 10);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm2k___2, vv_i__jm2k__.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 11);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp2_2, vv_i__j__kp2.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 12);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km2_2, vv_i__j__km2.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 13);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__k___3, vv_i__j__k__.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 22);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip1j__k___3, vv_ip1j__k__.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 23);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im1j__k___3, vv_im1j__k__.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 24);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp1k___3, vv_i__jp1k__.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 25);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm1k___3, vv_i__jm1k__.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 26);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp1_3, vv_i__j__kp1.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 27);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km1_3, vv_i__j__km1.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 28);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip2j__k___3, vv_ip2j__k__.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 8);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im2j__k___3, vv_im2j__k__.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 9);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp2k___3, vv_i__jp2k__.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 10);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm2k___3, vv_i__jm2k__.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 11);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp2_3, vv_i__j__kp2.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 12);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km2_3, vv_i__j__km2.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 13);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__k___4, vv_i__j__k__.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 29);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip1j__k___4, vv_ip1j__k__.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 30);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im1j__k___4, vv_im1j__k__.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 31);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp1k___4, vv_i__jp1k__.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 32);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm1k___4, vv_i__jm1k__.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 33);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp1_4, vv_i__j__kp1.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 34);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km1_4, vv_i__j__km1.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 35);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip2j__k___4, vv_ip2j__k__.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 8);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im2j__k___4, vv_im2j__k__.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 9);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp2k___4, vv_i__jp2k__.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 10);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm2k___4, vv_i__jm2k__.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 11);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp2_4, vv_i__j__kp2.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 12);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km2_4, vv_i__j__km2.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 13);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__k___5, vv_i__j__k__.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 36);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip1j__k___5, vv_ip1j__k__.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 37);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im1j__k___5, vv_im1j__k__.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 38);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp1k___5, vv_i__jp1k__.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 39);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm1k___5, vv_i__jm1k__.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 40);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp1_5, vv_i__j__kp1.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 41);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km1_5, vv_i__j__km1.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 42);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip2j__k___5, vv_ip2j__k__.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 8);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im2j__k___5, vv_im2j__k__.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 9);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp2k___5, vv_i__jp2k__.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 10);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm2k___5, vv_i__jm2k__.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 11);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp2_5, vv_i__j__kp2.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 12);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km2_5, vv_i__j__km2.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 13);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__k___6, vv_i__j__k__.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 43);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip1j__k___6, vv_ip1j__k__.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 44);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im1j__k___6, vv_im1j__k__.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 45);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp1k___6, vv_i__jp1k__.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 46);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm1k___6, vv_i__jm1k__.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 47);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp1_6, vv_i__j__kp1.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 48);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km1_6, vv_i__j__km1.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 49);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip2j__k___6, vv_ip2j__k__.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 8);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im2j__k___6, vv_im2j__k__.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 9);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp2k___6, vv_i__jp2k__.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 10);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm2k___6, vv_i__jm2k__.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 11);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp2_6, vv_i__j__kp2.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 12);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km2_6, vv_i__j__km2.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 13);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__k___7, vv_i__j__k__.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 50);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip1j__k___7, vv_ip1j__k__.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 51);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im1j__k___7, vv_im1j__k__.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 52);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp1k___7, vv_i__jp1k__.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 53);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm1k___7, vv_i__jm1k__.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 54);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp1_7, vv_i__j__kp1.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 55);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km1_7, vv_i__j__km1.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 56);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_ip2j__k___7, vv_ip2j__k__.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 8);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_im2j__k___7, vv_im2j__k__.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 9);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jp2k___7, vv_i__jp2k__.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 10);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__jm2k___7, vv_i__jm2k__.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 11);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__kp2_7, vv_i__j__kp2.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 12);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_i(grid_dv->vv_i__j__km2_7, vv_i__j__km2.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 13);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_a_0,  vv_fac_a.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 57);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_b_0,  vv_fac_b.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 58);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_c_0,  vv_fac_c.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 59);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_f_0,  vv_fac_f.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 60);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0v_0,    vv_T0v.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 61);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0r_0,    vv_T0r.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 62);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0t_0,    vv_T0t.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 63);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0p_0,    vv_T0p.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 64);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_fun_0,    vv_fun.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 65);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(grid_dv->vv_change_0, vv_change.at(0), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 66);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_a_1,  vv_fac_a.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 67);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_b_1,  vv_fac_b.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 68);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_c_1,  vv_fac_c.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 69);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_f_1,  vv_fac_f.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 70);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0v_1,    vv_T0v.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 71);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0r_1,    vv_T0r.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 72);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0t_1,    vv_T0t.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 73);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0p_1,    vv_T0p.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 74);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_fun_1,    vv_fun.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 75);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(grid_dv->vv_change_1, vv_change.at(1), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 76);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_a_2,  vv_fac_a.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 77);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_b_2,  vv_fac_b.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 78);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_c_2,  vv_fac_c.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 79);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_f_2,  vv_fac_f.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 80);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0v_2,    vv_T0v.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 81);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0r_2,    vv_T0r.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 82);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0t_2,    vv_T0t.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 83);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0p_2,    vv_T0p.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 84);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_fun_2,    vv_fun.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 85);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(grid_dv->vv_change_2, vv_change.at(2), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 86);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_a_3,  vv_fac_a.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 77);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_b_3,  vv_fac_b.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 78);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_c_3,  vv_fac_c.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 79);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_f_3,  vv_fac_f.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 80);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0v_3,    vv_T0v.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 81);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0r_3,    vv_T0r.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 82);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0t_3,    vv_T0t.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 83);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0p_3,    vv_T0p.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 84);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_fun_3,    vv_fun.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 85);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(grid_dv->vv_change_3, vv_change.at(3), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 86);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_a_4,  vv_fac_a.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 77);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_b_4,  vv_fac_b.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 78);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_c_4,  vv_fac_c.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 79);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_f_4,  vv_fac_f.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 80);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0v_4,    vv_T0v.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 81);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0r_4,    vv_T0r.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 82);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0t_4,    vv_T0t.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 83);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0p_4,    vv_T0p.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 84);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_fun_4,    vv_fun.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 85);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(grid_dv->vv_change_4, vv_change.at(4), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 86);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_a_5,  vv_fac_a.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 77);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_b_5,  vv_fac_b.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 78);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_c_5,  vv_fac_c.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 79);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_f_5,  vv_fac_f.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 80);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0v_5,    vv_T0v.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 81);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0r_5,    vv_T0r.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 82);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0t_5,    vv_T0t.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 83);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0p_5,    vv_T0p.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 84);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_fun_5,    vv_fun.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 85);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(grid_dv->vv_change_5, vv_change.at(5), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 86);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_a_6,  vv_fac_a.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 77);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_b_6,  vv_fac_b.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 78);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_c_6,  vv_fac_c.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 79);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_f_6,  vv_fac_f.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 80);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0v_6,    vv_T0v.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 81);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0r_6,    vv_T0r.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 82);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0t_6,    vv_T0t.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 83);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0p_6,    vv_T0p.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 84);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_fun_6,    vv_fun.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 85);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(grid_dv->vv_change_6, vv_change.at(6), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 86);

    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_a_7,  vv_fac_a.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 77);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_b_7,  vv_fac_b.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 78);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_c_7,  vv_fac_c.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 79);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv( grid_dv->vv_fac_f_7,  vv_fac_f.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 80);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0v_7,    vv_T0v.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 81);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0r_7,    vv_T0r.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 82);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0t_7,    vv_T0t.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 83);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_T0p_7,    vv_T0p.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 84);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(   grid_dv->vv_fun_7,    vv_fun.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 85);
    print_CUDA_error_if_any(allocate_and_copy_host_to_device_flattened_cv(grid_dv->vv_change_7, vv_change.at(7), grid_dv->n_nodes_total_host, grid_dv->n_nodes_on_levels_host), 86);

    // allocate tau
    print_CUDA_error_if_any(allocate_memory_on_device_cv((void**)&(grid_dv->tau), grid_dv->n_nodes_total_host), 87);

}


void cuda_finalize_grid(Grid_on_device* grid_dv){
    // deallocate memory on device
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->n_nodes_on_levels), 10000);
    delete [] grid_dv->n_nodes_on_levels_host;
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__k___0), 1);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_ip1j__k___0), 2);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_im1j__k___0), 3);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jp1k___0), 4);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jm1k___0), 5);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__kp1_0), 6);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__km1_0), 7);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__k___1), 8);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_ip1j__k___1), 9);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_im1j__k___1), 10);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jp1k___1), 11);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jm1k___1), 12);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__kp1_1), 13);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__km1_1), 14);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__k___2), 15);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_ip1j__k___2), 16);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_im1j__k___2), 17);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jp1k___2), 18);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jm1k___2), 19);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__kp1_2), 20);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__km1_2), 21);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__k___3), 22);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_ip1j__k___3), 23);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_im1j__k___3), 24);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jp1k___3), 25);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jm1k___3), 26);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__kp1_3), 27);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__km1_3), 28);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__k___4), 29);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_ip1j__k___4), 30);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_im1j__k___4), 31);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jp1k___4), 32);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jm1k___4), 33);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__kp1_4), 34);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__km1_4), 35);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__k___5), 36);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_ip1j__k___5), 37);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_im1j__k___5), 38);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jp1k___5), 39);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jm1k___5), 40);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__kp1_5), 41);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__km1_5), 42);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__k___6), 43);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_ip1j__k___6), 44);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_im1j__k___6), 45);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jp1k___6), 46);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jm1k___6), 47);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__kp1_6), 48);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__km1_6), 49);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__k___7), 50);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_ip1j__k___7), 51);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_im1j__k___7), 52);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jp1k___7), 53);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jm1k___7), 54);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__kp1_7), 55);
    print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__km1_7), 56);

    if(grid_dv->if_3rd_order){
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_ip2j__k___0), 10008);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_im2j__k___0), 10009);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jp2k___0), 10010);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jm2k___0), 10011);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__kp2_0), 10012);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__km2_0), 10013);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_ip2j__k___1), 10008);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_im2j__k___1), 10009);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jp2k___1), 10010);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jm2k___1), 10011);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__kp2_1), 10012);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__km2_1), 10013);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_ip2j__k___2), 10008);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_im2j__k___2), 10009);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jp2k___2), 10010);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jm2k___2), 10011);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__kp2_2), 10012);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__km2_2), 10013);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_ip2j__k___3), 10008);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_im2j__k___3), 10009);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jp2k___3), 10010);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jm2k___3), 10011);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__kp2_3), 10012);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__km2_3), 10013);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_ip2j__k___4), 10008);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_im2j__k___4), 10009);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jp2k___4), 10010);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jm2k___4), 10011);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__kp2_4), 10012);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__km2_4), 10013);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_ip2j__k___5), 10008);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_im2j__k___5), 10009);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jp2k___5), 10010);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jm2k___5), 10011);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__kp2_5), 10012);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__km2_5), 10013);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_ip2j__k___6), 10008);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_im2j__k___6), 10009);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jp2k___6), 10010);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jm2k___6), 10011);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__kp2_6), 10012);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__km2_6), 10013);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_ip2j__k___7), 10008);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_im2j__k___7), 10009);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jp2k___7), 10010);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__jm2k___7), 10011);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__kp2_7), 10012);
        print_CUDA_error_if_any(deallocate_memory_on_device_i(grid_dv->vv_i__j__km2_7), 10013);
    }

    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_a_0), 10057);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_b_0), 10058);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_c_0), 10059);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_f_0), 10060);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0v_0), 10061);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0r_0), 10062);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0t_0), 10063);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0p_0), 10064);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_fun_0), 10065);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_dv->vv_change_0), 10066);

    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_a_1), 10067);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_b_1), 10068);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_c_1), 10069);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_f_1), 10070);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0v_1), 10071);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0r_1), 10072);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0t_1), 10073);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0p_1), 10074);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_fun_1), 10075);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_dv->vv_change_1), 10076);

    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_a_2), 10077);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_b_2), 10078);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_c_2), 10079);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_f_2), 10080);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0v_2), 10081);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0r_2), 10082);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0t_2), 10083);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0p_2), 10084);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_fun_2), 10085);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_dv->vv_change_2), 10086);

    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_a_3), 10077);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_b_3), 10078);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_c_3), 10079);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_f_3), 10080);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0v_3), 10081);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0r_3), 10082);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0t_3), 10083);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0p_3), 10084);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_fun_3), 10085);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_dv->vv_change_3), 10086);

    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_a_4), 10077);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_b_4), 10078);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_c_4), 10079);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_f_4), 10080);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0v_4), 10081);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0r_4), 10082);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0t_4), 10083);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0p_4), 10084);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_fun_4), 10085);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_dv->vv_change_4), 10086);

    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_a_5), 10077);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_b_5), 10078);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_c_5), 10079);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_f_5), 10080);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0v_5), 10081);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0r_5), 10082);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0t_5), 10083);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0p_5), 10084);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_fun_5), 10085);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_dv->vv_change_5), 10086);

    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_a_6), 10077);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_b_6), 10078);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_c_6), 10079);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_f_6), 10080);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0v_6), 10081);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0r_6), 10082);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0t_6), 10083);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0p_6), 10084);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_fun_6), 10085);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_dv->vv_change_6), 10086);

    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_a_7), 10077);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_b_7), 10078);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_c_7), 10079);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv( grid_dv->vv_fac_f_7), 10080);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0v_7), 10081);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0r_7), 10082);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0t_7), 10083);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_T0p_7), 10084);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(   grid_dv->vv_fun_7), 10085);
    print_CUDA_error_if_any(deallocate_memory_on_device_cv(grid_dv->vv_change_7), 10086);


}


// copy tau from host to device
void cuda_copy_tau_to_device(Grid_on_device* grid_dv, CUSTOMREAL* tau_h){
    print_CUDA_error_if_any(copy_host_to_device_cv(grid_dv->tau, tau_h, grid_dv->n_nodes_total_host), 10087);
}


// copy tau from device to host
void cuda_copy_tau_to_host(Grid_on_device* grid_dv, CUSTOMREAL* tau_h){
    print_CUDA_error_if_any(copy_device_to_host_cv(tau_h, grid_dv->tau, grid_dv->n_nodes_total_host), 10088);
}

