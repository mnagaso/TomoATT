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

    // print length of ijk_for_this_subproc_optim
    //std::cout << "ijk_for_this_subproc_optim.size() = " << ijk_for_this_subproc_optim.size() << std::endl;
    //for (auto& ijk : ijk_for_this_subproc_optim) {
    //    std::cout << "ijk.size() = " << ijk.size() << std::endl;
    //    for (auto& ijk2 : ijk) {
    //        std::cout << "ijk2.size() = " << ijk2.size() << std::endl;
    //    }
    //}

    int iip, jjt, kkr;

    // measure time for only loop
    //auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<std::vector<int>>> ijk_for_this_subproc_optim_tmp ;

    ijk_for_this_subproc_optim_tmp = ijk_for_this_subproc_swps.at(iswp);
    int n_levels = ijk_for_this_subproc_optim_tmp.size();
    for (int i_level = 0; i_level < n_levels; i_level++) {
        int n_nodes = ijk_for_this_subproc_optim_tmp[i_level].size();
        //std::cout << "n_nodes = " << n_nodes << std::endl;

        int num_iter = n_nodes / NSIMD + (n_nodes % NSIMD == 0 ? 0 : 1);

        // load data of all nodes in one level on temporal aligned array
        for (int i = 0; i < num_iter; i++) {

            for (int ii = 0; ii < NSIMD; ii++){
                int i_node = i*NSIMD + ii;

                //if (i_node >= n_nodes)
                //    break;

                if (i_node < n_nodes){
//                    iip = ijk_for_this_subproc_optim_tmp[i_level][i_node][0];
//                    jjt = ijk_for_this_subproc_optim_tmp[i_level][i_node][1];
//                    kkr = ijk_for_this_subproc_optim_tmp[i_level][i_node][2];
                    iip = ijk_for_this_subproc_optim_tmp.at(i_level).at(i_node).at(0);
                    jjt = ijk_for_this_subproc_optim_tmp.at(i_level).at(i_node).at(1);
                    kkr = ijk_for_this_subproc_optim_tmp.at(i_level).at(i_node).at(2);
                } else {
                    // for the case that n_nodes is not multiple of NSIMD
                    iip = 0;
                    jjt = 0;
                    kkr = 0;
                }

                //std::cout << "load_stencil_data, iip = " << iip << ", jjt = " << jjt << ", kkr = " << kkr << ", I2V(iip, jjt, kkr) = " << I2V(iip, jjt, kkr) << ", i_node " << i_node << std::endl;
                dump_iip[ii] = (CUSTOMREAL)iip;
                dump_jjt[ii] = (CUSTOMREAL)jjt;
                dump_kkr[ii] = (CUSTOMREAL)kkr;

                dump_icc[ii] = iip;
                dump_jcc[ii] = jjt;
                dump_kcc[ii] = kkr;
                dump_ip1[ii] = (iip+1);
                dump_im1[ii] = (iip-1);
                dump_jp1[ii] = (jjt+1);
                dump_jm1[ii] = (jjt-1);
                dump_kp1[ii] = (kkr+1);
                dump_km1[ii] = (kkr-1);
                dump_ip2[ii] = (iip+2);
                dump_im2[ii] = (iip-2);
                dump_jp2[ii] = (jjt+2);
                dump_jm2[ii] = (jjt-2);
                dump_kp2[ii] = (kkr+2);
                dump_km2[ii] = (kkr-2);

//                // load all data at the same time
//                //load_stencil_data(grid.tau_loc, grid.fac_a_loc, grid.fac_b_loc, grid.fac_c_loc, grid.fac_f_loc, \
//                //                  grid.T0v_loc, grid.T0r_loc, grid.T0t_loc, grid.T0p_loc, grid.fun_loc, grid.is_changed, \
//                //                  iip, jjt, kkr, ii, \
//                //                  dump_c__, dump_fac_a, dump_fac_b, dump_fac_c, dump_fac_f, \
//                //                  dump_T0v, dump_T0r, dump_T0t, dump_T0p, dump_fun, dump_change, \
//                //                  dump_p__, dump_m__, dump__p_, dump__m_, dump___p, dump___m, dump_pp____, dump_mm____, dump___pp__, dump___mm__, dump_____pp, dump_____mm, \
//                //                  loc_I, loc_J, loc_K);
//
//                // load each data separately
//                //load_stencil_tau(
//                //    grid.tau_loc, \
//                //    iip, jjt, kkr, ii, \
//                //    dump_c__, dump_p__, dump_m__, dump__p_, dump__m_, dump___p, dump___m, dump_pp____, dump_mm____, dump___pp__, dump___mm__, dump_____pp, dump_____mm, \
//                //    loc_I, loc_J, loc_K);
//                int ip1 = iip+1, im1 = iip-1, jp1 = jjt+1, jm1 = jjt-1, kp1 = kkr+1, km1 = kkr-1;
//                int ip2 = iip+2, im2 = iip-2, jp2 = jjt+2, jm2 = jjt-2, kp2 = kkr+2, km2 = kkr-2;
//                load_mem_gen(grid.tau_loc, dump_c__   , iip, jjt, kkr, ii  );
//                load_mem_gen(grid.tau_loc, dump_p__   , ip1, jjt, kkr, ii);
//                load_mem_gen(grid.tau_loc, dump_m__   , im1, jjt, kkr, ii);
//                load_mem_gen(grid.tau_loc, dump__p_   , iip, jp1, kkr, ii);
//                load_mem_gen(grid.tau_loc, dump__m_   , iip, jm1, kkr, ii);
//                load_mem_gen(grid.tau_loc, dump___p   , iip, jjt, kp1, ii);
//                load_mem_gen(grid.tau_loc, dump___m   , iip, jjt, km1, ii);
//                load_mem_gen(grid.tau_loc, dump_pp____, ip2, jjt, kkr, ii);
//                load_mem_gen(grid.tau_loc, dump_mm____, im2, jjt, kkr, ii);
//                load_mem_gen(grid.tau_loc, dump___pp__, iip, jp2, kkr, ii);
//                load_mem_gen(grid.tau_loc, dump___mm__, iip, jm2, kkr, ii);
//                load_mem_gen(grid.tau_loc, dump_____pp, iip, jjt, kp2, ii);
//                load_mem_gen(grid.tau_loc, dump_____mm, iip, jjt, km2, ii);
//
//                load_mem_gen(grid.fac_a_loc, dump_fac_a, iip, jjt, kkr, ii);
//                load_mem_gen(grid.fac_b_loc, dump_fac_b, iip, jjt, kkr, ii);
//                load_mem_gen(grid.fac_c_loc, dump_fac_c, iip, jjt, kkr, ii);
//                load_mem_gen(grid.fac_f_loc, dump_fac_f, iip, jjt, kkr, ii);
//                load_mem_gen(grid.T0v_loc, dump_T0v, iip, jjt, kkr, ii);
//                load_mem_gen(grid.T0r_loc, dump_T0r, iip, jjt, kkr, ii);
//                load_mem_gen(grid.T0t_loc, dump_T0t, iip, jjt, kkr, ii);
//                load_mem_gen(grid.T0p_loc, dump_T0p, iip, jjt, kkr, ii);
//                load_mem_gen(grid.fun_loc, dump_fun, iip, jjt, kkr, ii);
//                load_mem_bool(grid.is_changed, dump_change, iip, jjt, kkr, ii);

            }

//            // alias for temporal aligned array
            __m256d* v_iip = (__m256d*)dump_iip;
            __m256d* v_jjt = (__m256d*)dump_jjt;
            __m256d* v_kkr = (__m256d*)dump_kkr;
//
//            __m256d* v_c__    = (__m256d*)dump_c__;
//            __m256d* v_p__    = (__m256d*)dump_p__;
//            __m256d* v_m__    = (__m256d*)dump_m__;
//            __m256d* v__p_    = (__m256d*)dump__p_;
//            __m256d* v__m_    = (__m256d*)dump__m_;
//            __m256d* v___p    = (__m256d*)dump___p;
//            __m256d* v___m    = (__m256d*)dump___m;
//            __m256d* v_pp____ = (__m256d*)dump_pp____;
//            __m256d* v_mm____ = (__m256d*)dump_mm____;
//            __m256d* v___pp__ = (__m256d*)dump___pp__;
//            __m256d* v___mm__ = (__m256d*)dump___mm__;
//            __m256d* v_____pp = (__m256d*)dump_____pp;
//            __m256d* v_____mm = (__m256d*)dump_____mm;
//
//            __m256d* v_fac_a  = (__m256d*)dump_fac_a;
//            __m256d* v_fac_b  = (__m256d*)dump_fac_b;
//            __m256d* v_fac_c  = (__m256d*)dump_fac_c;
//            __m256d* v_fac_f  = (__m256d*)dump_fac_f;
//            __m256d* v_T0v    = (__m256d*)dump_T0v;
//            __m256d* v_T0r    = (__m256d*)dump_T0r;
//            __m256d* v_T0t    = (__m256d*)dump_T0t;
//            __m256d* v_T0p    = (__m256d*)dump_T0p;
//            __m256d* v_fun    = (__m256d*)dump_fun;
//            __m256d* v_change = (__m256d*)dump_change;
//
//            // loop over all nodes in one level
//            int i_vec=0;
//            fake_stencil_3rd_pre_simd(v_iip[i_vec], v_jjt[i_vec], v_kkr[i_vec], v_c__[i_vec], v_p__[i_vec], v_m__[i_vec], v__p_[i_vec], v__m_[i_vec], v___p[i_vec], v___m[i_vec], \
//                                      v_pp____[i_vec], v_mm____[i_vec], v___pp__[i_vec], v___mm__[i_vec], v_____pp[i_vec], v_____mm[i_vec], \
//                                      v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
//                                      dp, dt, dr, loc_I, loc_J, loc_J);
//
//            //if (i_level == 198 && i_vec == 780)
//            //    std::cout << "nan will be happened here." << std::endl;
//
//            //// calculate updated value on c
//            fake_stencil_3rd_apre_simd(v_c__[i_vec], v_fac_a[i_vec], v_fac_b[i_vec], v_fac_c[i_vec], v_fac_f[i_vec], \
//                                       v_T0v[i_vec], v_T0p[i_vec], v_T0t[i_vec], v_T0r[i_vec], v_fun[i_vec], v_change[i_vec], \
//                                       v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
//                                       dp, dt, dr);


            __m256d v_c__    = load_mem_gen_to_m256d(grid.tau_loc,  dump_icc, dump_jcc, dump_kcc);
            __m256d v_p__    = load_mem_gen_to_m256d(grid.tau_loc,  dump_ip1, dump_jcc, dump_kcc);
            __m256d v_m__    = load_mem_gen_to_m256d(grid.tau_loc,  dump_im1, dump_jcc, dump_kcc);
            __m256d v__p_    = load_mem_gen_to_m256d(grid.tau_loc,  dump_icc, dump_jp1, dump_kcc);
            __m256d v__m_    = load_mem_gen_to_m256d(grid.tau_loc,  dump_icc, dump_jm1, dump_kcc);
            __m256d v___p    = load_mem_gen_to_m256d(grid.tau_loc,  dump_icc, dump_jcc, dump_kp1);
            __m256d v___m    = load_mem_gen_to_m256d(grid.tau_loc,  dump_icc, dump_jcc, dump_km1);
            __m256d v_pp____ = load_mem_gen_to_m256d(grid.tau_loc,  dump_ip2, dump_jcc, dump_kcc);
            __m256d v_mm____ = load_mem_gen_to_m256d(grid.tau_loc,  dump_im2, dump_jcc, dump_kcc);
            __m256d v___pp__ = load_mem_gen_to_m256d(grid.tau_loc,  dump_icc, dump_jp2, dump_kcc);
            __m256d v___mm__ = load_mem_gen_to_m256d(grid.tau_loc,  dump_icc, dump_jm2, dump_kcc);
            __m256d v_____pp = load_mem_gen_to_m256d(grid.tau_loc,  dump_icc, dump_jcc, dump_kp2);
            __m256d v_____mm = load_mem_gen_to_m256d(grid.tau_loc,  dump_icc, dump_jcc, dump_km2);

            // those parameters can be pre-loaded only once before the loop
            __m256d v_fac_a  = load_mem_gen_to_m256d(grid.fac_a_loc,  dump_icc, dump_jcc, dump_kcc);
            __m256d v_fac_b  = load_mem_gen_to_m256d(grid.fac_b_loc,  dump_icc, dump_jcc, dump_kcc);
            __m256d v_fac_c  = load_mem_gen_to_m256d(grid.fac_c_loc,  dump_icc, dump_jcc, dump_kcc);
            __m256d v_fac_f  = load_mem_gen_to_m256d(grid.fac_f_loc,  dump_icc, dump_jcc, dump_kcc);
            __m256d v_T0v    = load_mem_gen_to_m256d(grid.T0v_loc,    dump_icc, dump_jcc, dump_kcc);
            __m256d v_T0r    = load_mem_gen_to_m256d(grid.T0r_loc,    dump_icc, dump_jcc, dump_kcc);
            __m256d v_T0t    = load_mem_gen_to_m256d(grid.T0t_loc,    dump_icc, dump_jcc, dump_kcc);
            __m256d v_T0p    = load_mem_gen_to_m256d(grid.T0p_loc,    dump_icc, dump_jcc, dump_kcc);
            __m256d v_fun    = load_mem_gen_to_m256d(grid.fun_loc,    dump_icc, dump_jcc, dump_kcc);
            __m256d v_change = load_mem_bool_to_m256d(grid.is_changed,  dump_icc, dump_jcc, dump_kcc);

            // loop over all nodes in one level
            fake_stencil_3rd_pre_simd(v_iip[0], v_jjt[0], v_kkr[0], v_c__, v_p__, v_m__, v__p_, v__m_, v___p, v___m, \
                                      v_pp____, v_mm____, v___pp__, v___mm__, v_____pp, v_____mm, \
                                      v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                      dp, dt, dr, loc_I, loc_J, loc_J);

            //if (i_level == 198 && i_vec == 780)
            //    std::cout << "nan will be happened here." << std::endl;

            //// calculate updated value on c
            fake_stencil_3rd_apre_simd(v_c__, v_fac_a, v_fac_b, v_fac_c, v_fac_f, \
                                       v_T0v, v_T0p  , v_T0t  , v_T0r  , v_fun  , v_change, \
                                       v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                       dp, dt, dr);

            // store v_c__ to dump_c__
            _mm256_store_pd(dump_c__, v_c__);

            for (int i = 0; i < NSIMD; i++) {
                int tmp_ijk = I2V(dump_iip[i], dump_jjt[i], dump_kkr[i]);
                grid.tau_loc[tmp_ijk] = dump_c__[i];
            }



        } // end of i_vec loop

        // mpi synchronization
        synchronize_all_sub();

    } // end of i_level loop

//    // count nan in grid.tau_loc
//    int count_nan = 0;
//    int nnodes = loc_I*loc_J*loc_K;
//    for (int i = 0; i < nnodes; i++) {
//        if (std::isnan(grid.tau_loc[i])) {
//            count_nan++;
//        }
//    }
//
//    std::cout << "count_nan = " << count_nan << std::endl;


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