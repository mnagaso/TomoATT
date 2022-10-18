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

    // if max_length % NSIMD != then we need to pad it

    __m256i v_iip;
    __m256i v_jjt;
    __m256i v_kkr;

    __m256d v_c__;
    __m256d v_p__;
    __m256d v_m__;
    __m256d v__p_;
    __m256d v__m_;
    __m256d v___p;
    __m256d v___m;
    __m256d v_pp____;
    __m256d v_mm____;
    __m256d v___pp__;
    __m256d v___mm__;
    __m256d v_____pp;
    __m256d v_____mm;

    __m256d v_fac_a;
    __m256d v_fac_b;
    __m256d v_fac_c;
    __m256d v_fac_f;
    __m256d v_T0v;
    __m256d v_T0r;
    __m256d v_T0t;
    __m256d v_T0p;
    __m256d v_fun;

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
        //for (int i_node = 0; i_node < n_nodes; i_node++) {
        for (int i = 0; i < num_iter; i++) {
            for (int ii = 0; ii < NSIMD; ii++){
                int i_node = i*NSIMD + ii;

                if (i_node >= n_nodes)
                    break;

                if (i_node < n_nodes){
                    iip = ijk_for_this_subproc_optim_tmp[i_level][i_node][0];
                    jjt = ijk_for_this_subproc_optim_tmp[i_level][i_node][1];
                    kkr = ijk_for_this_subproc_optim_tmp[i_level][i_node][2];
                } else {
                    iip = 0;
                    jjt = 0;
                    kkr = 0;
                }

                //std::cout << "load_stencil_data, iip = " << iip << ", jjt = " << jjt << ", kkr = " << kkr << ", I2V(iip, jjt, kkr) = " << I2V(iip, jjt, kkr) << ", i_node " << i_node << std::endl;
                dump_iip[i_node] = iip;
                dump_jjt[i_node] = jjt;
                dump_kkr[i_node] = kkr;


                load_stencil_data(grid.tau_loc, grid.fac_a_loc, grid.fac_b_loc, grid.fac_c_loc, grid.fac_f_loc, \
                                    grid.T0v_loc, grid.T0r_loc, grid.T0t_loc, grid.T0p_loc, grid.fun_loc, \
                                    iip, jjt, kkr, i_node, \
                                    dump_c__, dump_fac_a, dump_fac_b, dump_fac_c, dump_fac_f, \
                                    dump_T0v, dump_T0r, dump_T0t, dump_T0p, dump_fun, \
                                    dump_p__, dump_m__, dump__p_, dump__m_, dump___p, dump___m, dump_pp____, dump_mm____, dump___pp__, dump___mm__, dump_____pp, dump_____mm, \
                                    loc_I, loc_J, loc_K);

            }
        }

        for (int i_node = 0; i_node < num_iter; i_node++) {
            // load stencil data

            v_iip = _mm256_load_si256((__m256i*)&dump_iip[i_node*NSIMD]);
            v_jjt = _mm256_load_si256((__m256i*)&dump_jjt[i_node*NSIMD]);
            v_kkr = _mm256_load_si256((__m256i*)&dump_kkr[i_node*NSIMD]);

            v_c__ = _mm256_load_pd(&dump_c__[i_node*NSIMD]);
            v_p__ = _mm256_load_pd(&dump_p__[i_node*NSIMD]);
            v_m__ = _mm256_load_pd(&dump_m__[i_node*NSIMD]);
            v__p_ = _mm256_load_pd(&dump__p_[i_node*NSIMD]);
            v__m_ = _mm256_load_pd(&dump__m_[i_node*NSIMD]);
            v___p = _mm256_load_pd(&dump___p[i_node*NSIMD]);
            v___m = _mm256_load_pd(&dump___m[i_node*NSIMD]);
            v_pp____ = _mm256_load_pd(&dump_pp____[i_node*NSIMD]);
            v_mm____ = _mm256_load_pd(&dump_mm____[i_node*NSIMD]);
            v___pp__ = _mm256_load_pd(&dump___pp__[i_node*NSIMD]);
            v___mm__ = _mm256_load_pd(&dump___mm__[i_node*NSIMD]);
            v_____pp = _mm256_load_pd(&dump_____pp[i_node*NSIMD]);
            v_____mm = _mm256_load_pd(&dump_____mm[i_node*NSIMD]);

            v_fac_a = _mm256_load_pd(&dump_fac_a[i_node*NSIMD]);
            v_fac_b = _mm256_load_pd(&dump_fac_b[i_node*NSIMD]);
            v_fac_c = _mm256_load_pd(&dump_fac_c[i_node*NSIMD]);
            v_fac_f = _mm256_load_pd(&dump_fac_f[i_node*NSIMD]);
            v_T0v = _mm256_load_pd(&dump_T0v[i_node*NSIMD]);
            v_T0r = _mm256_load_pd(&dump_T0r[i_node*NSIMD]);
            v_T0t = _mm256_load_pd(&dump_T0t[i_node*NSIMD]);
            v_T0p = _mm256_load_pd(&dump_T0p[i_node*NSIMD]);
            v_fun = _mm256_load_pd(&dump_fun[i_node*NSIMD]);


            // store stencil coefs
            v_pp1 = _mm256_setzero_pd();
            v_pp2 = _mm256_setzero_pd();
            v_pt1 = _mm256_setzero_pd();
            v_pt2 = _mm256_setzero_pd();
            v_pr1 = _mm256_setzero_pd();
            v_pr2 = _mm256_setzero_pd();

            // calculate stencil coefs
            fake_stencil_3rd_pre_simd(v_iip, v_jjt, v_kkr, v_c__, v_p__, v_m__, v__p_, v__m_, v___p, v___m, \
                                      v_pp____, v_mm____, v___pp__, v___mm__, v_____pp, v_____mm, \
                                      v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                      dp, dt, dr, loc_I, loc_J, loc_J);

            //// calculate updated value on c
            fake_stencil_3rd_apre_simd(v_c__, v_fac_a, v_fac_b, v_fac_c, v_fac_f, \
                                       v_T0v, v_T0r, v_T0t, v_T0p, v_fun, \
                                       v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                       dp, dt, dr);

            // unload v_c__ to dump_c__
            _mm256_store_pd(&dump_c__[i_node*NSIMD], v_c__);

        }


        for (int i = 0; i < n_nodes; i++) {

            int tmp_ijk = I2V(dump_iip[i], dump_jjt[i], dump_kkr[i]);
            grid.tau_loc[tmp_ijk] = dump_c__[i];
        }


        // below is the original code without intrinsics

        //for (int i_node = 0; i_node < n_nodes; i_node++) {
        //    int tmp_ijk = ijk_for_this_subproc_optim_tmp.at(i_level).at(i_node);
        //    V2I(tmp_ijk, iip, jjt, kkr);
        //    fake_stencil_3rd_pre(c, iip, jjt, kkr, dump_pp1, dump_pp2, dump_pt1, dump_pt2, dump_pr1, dump_pr2, i_node);
        //}

        //for (int i_node = 0; i_node < n_nodes; i_node++) {
        //    int tmp_ijk = ijk_for_this_subproc_optim_tmp.at(i_level).at(i_node);
        //    V2I(tmp_ijk, iip, jjt, kkr);
        //    fake_stencil_3rd_apre(a,b,c, iip, jjt, kkr, dump_pp1, dump_pp2, dump_pt1, dump_pt2, dump_pr1, dump_pr2, i_node);
        //}
    }


#endif // __AVX__

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