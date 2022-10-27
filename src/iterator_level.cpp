#include "iterator_level.h"

#ifdef USE_AVX
#include "vectorized_sweep.h"
#endif

Iterator_level::Iterator_level(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                : Iterator(IP, grid, src, io, first_init, is_teleseismic_in, is_second_run_in) {
    // do nothing
}


void Iterator_level::do_sweep_adj(int iswp, Grid& grid, InputParams& IP){

    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;
    int n_levels = ijk_for_this_subproc.size();

    for (int i_level = 0; i_level < n_levels; i_level++) {
        size_t n_nodes = ijk_for_this_subproc[i_level].size();

        #pragma omp simd
        for (size_t i_node = 0; i_node < n_nodes; i_node++) {

            V2I(ijk_for_this_subproc[i_level][i_node], iip, jjt, kkr);

            if (r_dirc < 0) kkr = nr-kkr; //kk-1;
            else            kkr = kkr-1;  //nr-kk;
            if (t_dirc < 0) jjt = nt-jjt; //jj-1;
            else            jjt = jjt-1;  //nt-jj;
            if (p_dirc < 0) iip = np-iip; //ii-1;
            else            iip = iip-1;  //np-ii;

            //
            // calculate stencils
            //
            calculate_stencil_adj(grid, iip, jjt, kkr);
        } // end ijk

        // mpi synchronization
        synchronize_all_sub();

    } // end loop i_level

}


Iterator_level_tele::Iterator_level_tele(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                : Iterator(IP, grid, src, io, first_init, is_teleseismic_in, is_second_run_in) {
    // do nothing
}


void Iterator_level_tele::do_sweep_adj(int iswp, Grid& grid, InputParams& IP){
    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;
    int n_levels = ijk_for_this_subproc.size();

    for (int i_level = 0; i_level < n_levels; i_level++) {
        size_t n_nodes = ijk_for_this_subproc[i_level].size();

        #pragma omp simd
        for (size_t i_node = 0; i_node < n_nodes; i_node++) {

            V2I(ijk_for_this_subproc[i_level][i_node], iip, jjt, kkr);

            if (r_dirc < 0) kkr = nr-kkr; //kk-1;
            else            kkr = kkr-1;  //nr-kk;
            if (t_dirc < 0) jjt = nt-jjt; //jj-1;
            else            jjt = jjt-1;  //nt-jj;
            if (p_dirc < 0) iip = np-iip; //ii-1;
            else            iip = iip-1;  //np-ii;

            //
            // calculate stencils
            //
            if (iip != 0    && jjt != 0    && kkr != 0 \
             && iip != np-1 && jjt != nt-1 && kkr != nr-1) {
                // calculate stencils
                calculate_stencil_adj(grid, iip, jjt, kkr);
            } else {
                calculate_boundary_nodes_tele_adj(grid, iip, jjt, kkr);
            }
        } // end ijk

        // mpi synchronization
        synchronize_all_sub();

    } // end loop i_level

}


Iterator_level_1st_order::Iterator_level_1st_order(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                         : Iterator_level(IP, grid, src, io, first_init, is_teleseismic_in, is_second_run_in) {
    // initialization is done in the base class
}


void Iterator_level_1st_order::do_sweep(int iswp, Grid& grid, InputParams& IP){
    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;
    int n_levels = ijk_for_this_subproc.size();

    for (int i_level = 0; i_level < n_levels; i_level++) {
        size_t n_nodes = ijk_for_this_subproc[i_level].size();

        #pragma omp simd
        for (size_t i_node = 0; i_node < n_nodes; i_node++) {

            V2I(ijk_for_this_subproc[i_level][i_node], iip, jjt, kkr);

            if (r_dirc < 0) kkr = nr-kkr; //kk-1;
            else            kkr = kkr-1;  //nr-kk;
            if (t_dirc < 0) jjt = nt-jjt; //jj-1;
            else            jjt = jjt-1;  //nt-jj;
            if (p_dirc < 0) iip = np-iip; //ii-1;
            else            iip = iip-1;  //np-ii;

            //
            // calculate stencils
            //
            if (grid.is_changed[I2V(iip, jjt, kkr)]) {
                calculate_stencil_1st_order(grid, iip, jjt, kkr);
            } // is_changed == true
        } // end ijk

        // mpi synchronization
        synchronize_all_sub();

    } // end loop i_level

    // update boundary
    if (subdom_main) {
        calculate_boundary_nodes(grid);
    }
}


Iterator_level_3rd_order::Iterator_level_3rd_order(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                         : Iterator_level(IP, grid, src, io, first_init, is_teleseismic_in, is_second_run_in) {
    // initialization is done in the base class
}


void Iterator_level_3rd_order::do_sweep(int iswp, Grid& grid, InputParams& IP){


#ifndef USE_AVX

    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;
    int n_levels = ijk_for_this_subproc.size();

    for (int i_level = 0; i_level < n_levels; i_level++) {
        size_t n_nodes = ijk_for_this_subproc[i_level].size();

        #pragma omp simd
        for (size_t i_node = 0; i_node < n_nodes; i_node++) {

            V2I(ijk_for_this_subproc[i_level][i_node], iip, jjt, kkr);

            if (r_dirc < 0) kkr = nr-kkr; //kk-1;
            else            kkr = kkr-1;  //nr-kk;
            if (t_dirc < 0) jjt = nt-jjt; //jj-1;
            else            jjt = jjt-1;  //nt-jj;
            if (p_dirc < 0) iip = np-iip; //ii-1;
            else            iip = iip-1;  //np-ii;

            //
            // calculate stencils
            //
            if (grid.is_changed[I2V(iip, jjt, kkr)]) {
                calculate_stencil_3rd_order(grid, iip, jjt, kkr);
            } // is_changed == true
        } // end ijk

        // mpi synchronization
        synchronize_all_sub();

    } // end loop i_level

#else // __AVX__

    // store stencil coefs
    __m256d v_pp1;
    __m256d v_pp2;
    __m256d v_pt1;
    __m256d v_pt2;
    __m256d v_pr1;
    __m256d v_pr2;

    // measure time for only loop
    //auto start = std::chrono::high_resolution_clock::now();

    int n_levels = ijk_for_this_subproc.size();
    for (int i_level = 0; i_level < n_levels; i_level++) {
        int n_nodes = ijk_for_this_subproc.at(i_level).size();
        //std::cout << "n_nodes = " << n_nodes << std::endl;

        int num_iter = n_nodes / NSIMD + (n_nodes % NSIMD == 0 ? 0 : 1);

        // make alias to preloaded data
        __m256d* v_iip    = (__m256d*) vv_iip.at(iswp).at(i_level);
        __m256d* v_jjt    = (__m256d*) vv_jjt.at(iswp).at(i_level);
        __m256d* v_kkr    = (__m256d*) vv_kkr.at(iswp).at(i_level);

        __m256d* v_fac_a  = (__m256d*) vv_fac_a.at(iswp).at(i_level);
        __m256d* v_fac_b  = (__m256d*) vv_fac_b.at(iswp).at(i_level);
        __m256d* v_fac_c  = (__m256d*) vv_fac_c.at(iswp).at(i_level);
        __m256d* v_fac_f  = (__m256d*) vv_fac_f.at(iswp).at(i_level);
        __m256d* v_T0v    = (__m256d*) vv_T0v.at(iswp).at(i_level);
        __m256d* v_T0r    = (__m256d*) vv_T0r.at(iswp).at(i_level);
        __m256d* v_T0t    = (__m256d*) vv_T0t.at(iswp).at(i_level);
        __m256d* v_T0p    = (__m256d*) vv_T0p.at(iswp).at(i_level);
        __m256d* v_fun    = (__m256d*) vv_fun.at(iswp).at(i_level);
        __m256d* v_change = (__m256d*) vv_change.at(iswp).at(i_level);

        // alias for dumped index
        int* dump_icc = vv_icc.at(iswp).at(i_level);
        int* dump_jcc = vv_jcc.at(iswp).at(i_level);
        int* dump_kcc = vv_kcc.at(iswp).at(i_level);
        int* dump_ip1 = vv_ip1.at(iswp).at(i_level);
        int* dump_jp1 = vv_jp1.at(iswp).at(i_level);
        int* dump_kp1 = vv_kp1.at(iswp).at(i_level);
        int* dump_im1 = vv_im1.at(iswp).at(i_level);
        int* dump_jm1 = vv_jm1.at(iswp).at(i_level);
        int* dump_km1 = vv_km1.at(iswp).at(i_level);
        int* dump_ip2 = vv_ip2.at(iswp).at(i_level);
        int* dump_jp2 = vv_jp2.at(iswp).at(i_level);
        int* dump_kp2 = vv_kp2.at(iswp).at(i_level);
        int* dump_im2 = vv_im2.at(iswp).at(i_level);
        int* dump_jm2 = vv_jm2.at(iswp).at(i_level);
        int* dump_km2 = vv_km2.at(iswp).at(i_level);

        // load data of all nodes in one level on temporal aligned array
        for (int i_vec = 0; i_vec < num_iter; i_vec++) {

            __m256d v_c__    = load_mem_gen_to_m256d(grid.tau_loc,  &dump_icc[i_vec], &dump_jcc[i_vec], &dump_kcc[i_vec]);
            __m256d v_p__    = load_mem_gen_to_m256d(grid.tau_loc,  &dump_ip1[i_vec], &dump_jcc[i_vec], &dump_kcc[i_vec]);
            __m256d v_m__    = load_mem_gen_to_m256d(grid.tau_loc,  &dump_im1[i_vec], &dump_jcc[i_vec], &dump_kcc[i_vec]);
            __m256d v__p_    = load_mem_gen_to_m256d(grid.tau_loc,  &dump_icc[i_vec], &dump_jp1[i_vec], &dump_kcc[i_vec]);
            __m256d v__m_    = load_mem_gen_to_m256d(grid.tau_loc,  &dump_icc[i_vec], &dump_jm1[i_vec], &dump_kcc[i_vec]);
            __m256d v___p    = load_mem_gen_to_m256d(grid.tau_loc,  &dump_icc[i_vec], &dump_jcc[i_vec], &dump_kp1[i_vec]);
            __m256d v___m    = load_mem_gen_to_m256d(grid.tau_loc,  &dump_icc[i_vec], &dump_jcc[i_vec], &dump_km1[i_vec]);
            __m256d v_pp____ = load_mem_gen_to_m256d(grid.tau_loc,  &dump_ip2[i_vec], &dump_jcc[i_vec], &dump_kcc[i_vec]);
            __m256d v_mm____ = load_mem_gen_to_m256d(grid.tau_loc,  &dump_im2[i_vec], &dump_jcc[i_vec], &dump_kcc[i_vec]);
            __m256d v___pp__ = load_mem_gen_to_m256d(grid.tau_loc,  &dump_icc[i_vec], &dump_jp2[i_vec], &dump_kcc[i_vec]);
            __m256d v___mm__ = load_mem_gen_to_m256d(grid.tau_loc,  &dump_icc[i_vec], &dump_jm2[i_vec], &dump_kcc[i_vec]);
            __m256d v_____pp = load_mem_gen_to_m256d(grid.tau_loc,  &dump_icc[i_vec], &dump_jcc[i_vec], &dump_kp2[i_vec]);
            __m256d v_____mm = load_mem_gen_to_m256d(grid.tau_loc,  &dump_icc[i_vec], &dump_jcc[i_vec], &dump_km2[i_vec]);

            // loop over all nodes in one level
            fake_stencil_3rd_pre_simd(v_iip[i_vec], v_jjt[i_vec], v_kkr[i_vec], v_c__, v_p__, v_m__, v__p_, v__m_, v___p, v___m, \
                                      v_pp____, v_mm____, v___pp__, v___mm__, v_____pp, v_____mm, \
                                      v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                      dp, dt, dr, loc_I, loc_J, loc_J);

            //// calculate updated value on c
            fake_stencil_3rd_apre_simd(v_c__, v_fac_a[i_vec], v_fac_b[i_vec], v_fac_c[i_vec], v_fac_f[i_vec], \
                                       v_T0v[i_vec], v_T0p[i_vec]  , v_T0t[i_vec]  , v_T0r[i_vec]  , v_fun[i_vec]  , v_change[i_vec], \
                                       v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                       dp, dt, dr);

            // store v_c__ to dump_c__
            _mm256_store_pd(dump_c__, v_c__);

            for (int i = 0; i < NSIMD; i++) {
                int tmp_ijk = I2V(dump_icc[i_vec*NSIMD+i], dump_jcc[i_vec*NSIMD+i], dump_kcc[i_vec*NSIMD+i]);
                grid.tau_loc[tmp_ijk] = dump_c__[i];
            }



        } // end of i_vec loop

        // mpi synchronization
        synchronize_all_sub();

    } // end of i_level loop


#endif // USE_AVX


    // update boundary
    if (subdom_main) {
        calculate_boundary_nodes(grid);
    }

}


Iterator_level_1st_order_tele::Iterator_level_1st_order_tele(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                         : Iterator_level_tele(IP, grid, src, io, first_init, is_teleseismic_in, is_second_run_in) {
    // initialization is done in the base class
}


void Iterator_level_1st_order_tele::do_sweep(int iswp, Grid& grid, InputParams& IP){
    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;
    int n_levels = ijk_for_this_subproc.size();

    for (int i_level = 0; i_level < n_levels; i_level++) {
        size_t n_nodes = ijk_for_this_subproc[i_level].size();

        #pragma omp simd
        for (size_t i_node = 0; i_node < n_nodes; i_node++) {

            V2I(ijk_for_this_subproc[i_level][i_node], iip, jjt, kkr);

            if (r_dirc < 0) kkr = nr-kkr; //kk-1;
            else            kkr = kkr-1;  //nr-kk;
            if (t_dirc < 0) jjt = nt-jjt; //jj-1;
            else            jjt = jjt-1;  //nt-jj;
            if (p_dirc < 0) iip = np-iip; //ii-1;
            else            iip = iip-1;  //np-ii;

            //
            // calculate stencils
            //
            if (iip != 0 && iip != np-1 && jjt != 0 && jjt != nt-1 && kkr != 0 && kkr != nr-1) {
                // calculate stencils
                calculate_stencil_1st_order_tele(grid, iip, jjt, kkr);
            } else {
                // update boundary
                calculate_boundary_nodes_tele(grid, iip, jjt, kkr);
            }
        } // end ijk

        // mpi synchronization
        synchronize_all_sub();

    } // end loop i_level
}


Iterator_level_3rd_order_tele::Iterator_level_3rd_order_tele(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                         : Iterator_level_tele(IP, grid, src, io, first_init, is_teleseismic_in, is_second_run_in) {
    // initialization is done in the base class
}


void Iterator_level_3rd_order_tele::do_sweep(int iswp, Grid& grid, InputParams& IP){
    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;
    int n_levels = ijk_for_this_subproc.size();

    for (int i_level = 0; i_level < n_levels; i_level++) {
        size_t n_nodes = ijk_for_this_subproc[i_level].size();

        #pragma omp simd
        for (size_t i_node = 0; i_node < n_nodes; i_node++) {

            V2I(ijk_for_this_subproc[i_level][i_node], iip, jjt, kkr);

            if (r_dirc < 0) kkr = nr-kkr; //kk-1;
            else            kkr = kkr-1;  //nr-kk;
            if (t_dirc < 0) jjt = nt-jjt; //jj-1;
            else            jjt = jjt-1;  //nt-jj;
            if (p_dirc < 0) iip = np-iip; //ii-1;
            else            iip = iip-1;  //np-ii;

            //
            // calculate stencils
            //
            if (iip != 0 && iip != np-1 && jjt != 0 && jjt != nt-1 && kkr != 0 && kkr != nr-1) {
                // calculate stencils
                calculate_stencil_3rd_order_tele(grid, iip, jjt, kkr);
            } else {
                // update boundary
                calculate_boundary_nodes_tele(grid, iip, jjt, kkr);
            }
        } // end ijk

        // mpi synchronization
        synchronize_all_sub();

    } // end loop i_level

}