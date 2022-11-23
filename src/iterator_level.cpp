#include "iterator_level.h"

#ifdef USE_SIMD
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

    if(!use_gpu){

#if !defined USE_SIMD

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

#elif USE_AVX512 || USE_AVX

        //
        __mTd v_DP_inv      = _mmT_set1_pd(1.0/dp);
        __mTd v_DT_inv      = _mmT_set1_pd(1.0/dt);
        __mTd v_DR_inv      = _mmT_set1_pd(1.0/dr);
        __mTd v_DP_inv_half = _mmT_set1_pd(1.0/dp*0.5);
        __mTd v_DT_inv_half = _mmT_set1_pd(1.0/dt*0.5);
        __mTd v_DR_inv_half = _mmT_set1_pd(1.0/dr*0.5);

        // store stencil coefs
        __mTd v_pp1;
        __mTd v_pp2;
        __mTd v_pt1;
        __mTd v_pt2;
        __mTd v_pr1;
        __mTd v_pr2;

        // measure time for only loop
        //auto start = std::chrono::high_resolution_clock::now();

        int n_levels = ijk_for_this_subproc.size();
        for (int i_level = 0; i_level < n_levels; i_level++) {
            int n_nodes = ijk_for_this_subproc.at(i_level).size();
            //std::cout << "n_nodes = " << n_nodes << std::endl;

            int num_iter = n_nodes / NSIMD + (n_nodes % NSIMD == 0 ? 0 : 1);

            // make alias to preloaded data
            __mTd* v_iip    = (__mTd*) vv_iip.at(iswp).at(i_level);
            __mTd* v_jjt    = (__mTd*) vv_jjt.at(iswp).at(i_level);
            __mTd* v_kkr    = (__mTd*) vv_kkr.at(iswp).at(i_level);

            __mTd* v_fac_a  = (__mTd*) vv_fac_a.at(iswp).at(i_level);
            __mTd* v_fac_b  = (__mTd*) vv_fac_b.at(iswp).at(i_level);
            __mTd* v_fac_c  = (__mTd*) vv_fac_c.at(iswp).at(i_level);
            __mTd* v_fac_f  = (__mTd*) vv_fac_f.at(iswp).at(i_level);
            __mTd* v_T0v    = (__mTd*) vv_T0v.at(iswp).at(i_level);
            __mTd* v_T0r    = (__mTd*) vv_T0r.at(iswp).at(i_level);
            __mTd* v_T0t    = (__mTd*) vv_T0t.at(iswp).at(i_level);
            __mTd* v_T0p    = (__mTd*) vv_T0p.at(iswp).at(i_level);
            __mTd* v_fun    = (__mTd*) vv_fun.at(iswp).at(i_level);
            __mTd* v_change = (__mTd*) vv_change.at(iswp).at(i_level);

            // alias for dumped index
            int* dump_ijk   = vv_i__j__k__.at(iswp).at(i_level);
            int* dump_ip1jk = vv_ip1j__k__.at(iswp).at(i_level);
            int* dump_im1jk = vv_im1j__k__.at(iswp).at(i_level);
            int* dump_ijp1k = vv_i__jp1k__.at(iswp).at(i_level);
            int* dump_ijm1k = vv_i__jm1k__.at(iswp).at(i_level);
            int* dump_ijkp1 = vv_i__j__kp1.at(iswp).at(i_level);
            int* dump_ijkm1 = vv_i__j__km1.at(iswp).at(i_level);

            // load data of all nodes in one level on temporal aligned array
            for (int _i_vec = 0; _i_vec < num_iter; _i_vec++) {

                int i_vec = _i_vec * NSIMD;
                __mTd v_c__    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijk[i_vec]);
                __mTd v_p__    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ip1jk[i_vec]);
                __mTd v_m__    = load_mem_gen_to_mTd(grid.tau_loc, &dump_im1jk[i_vec]);
                __mTd v__p_    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijp1k[i_vec]);
                __mTd v__m_    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijm1k[i_vec]);
                __mTd v___p    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijkp1[i_vec]);
                __mTd v___m    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijkm1[i_vec]);

                // loop over all nodes in one level
                vect_stencil_1st_pre_simd(v_iip[_i_vec], v_jjt[_i_vec], v_kkr[_i_vec], \
                                          v_c__, \
                                          v_p__,    v_m__,    v__p_,    v__m_,    v___p,    v___m, \
                                          v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                          v_DP_inv, v_DT_inv, v_DR_inv, \
                                          v_DP_inv_half, v_DT_inv_half, v_DR_inv_half, \
                                          loc_I, loc_J, loc_K);

                // calculate updated value on c
                vect_stencil_1st_3rd_apre_simd(v_c__, v_fac_a[_i_vec], v_fac_b[_i_vec], v_fac_c[_i_vec], v_fac_f[_i_vec], \
                                               v_T0v[_i_vec], v_T0p[_i_vec]  , v_T0t[_i_vec]  , v_T0r[_i_vec]  , v_fun[_i_vec]  , v_change[_i_vec], \
                                               v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                               v_DP_inv, v_DT_inv, v_DR_inv);

                // store v_c__ to dump_c__
                _mmT_store_pd(dump_c__, v_c__);


                for (int i = 0; i < NSIMD; i++) {
                    if(i_vec+i>=n_nodes) break;

                    grid.tau_loc[dump_ijk[i_vec+i]] = dump_c__[i];
                }



            } // end of i_vec loop

            // mpi synchronization
            synchronize_all_sub();

        } // end of i_level loop

#elif USE_ARM_SVE

        svbool_t pg;
        //
        __mTd v_DP_inv      = svdup_f64(1.0/dp);
        __mTd v_DT_inv      = svdup_f64(1.0/dt);
        __mTd v_DR_inv      = svdup_f64(1.0/dr);
        __mTd v_DP_inv_half = svdup_f64(1.0/dp*0.5);
        __mTd v_DT_inv_half = svdup_f64(1.0/dt*0.5);
        __mTd v_DR_inv_half = svdup_f64(1.0/dr*0.5);

        // store stencil coefs
        __mTd v_pp1;
        __mTd v_pp2;
        __mTd v_pt1;
        __mTd v_pt2;
        __mTd v_pr1;
        __mTd v_pr2;

        __mTd v_c__   ;
        __mTd v_p__   ;
        __mTd v_m__   ;
        __mTd v__p_   ;
        __mTd v__m_   ;
        __mTd v___p   ;
        __mTd v___m   ;

        __mTd v_iip_   ;
        __mTd v_jjt_   ;
        __mTd v_kkr_   ;
        __mTd v_fac_a_ ;
        __mTd v_fac_b_ ;
        __mTd v_fac_c_ ;
        __mTd v_fac_f_ ;
        __mTd v_T0v_   ;
        __mTd v_T0r_   ;
        __mTd v_T0t_   ;
        __mTd v_T0p_   ;
        __mTd v_fun_   ;
        __mTd v_change_;


        // measure time for only loop
        //auto start = std::chrono::high_resolution_clock::now();

        int n_levels = ijk_for_this_subproc.size();
        for (int i_level = 0; i_level < n_levels; i_level++) {
            int n_nodes = ijk_for_this_subproc.at(i_level).size();
            //std::cout << "n_nodes = " << n_nodes << std::endl;

            int num_iter = n_nodes / NSIMD + (n_nodes % NSIMD == 0 ? 0 : 1);

            // make alias to preloaded data
            CUSTOMREAL* v_iip    = vv_iip.at(iswp).at(i_level);
            CUSTOMREAL* v_jjt    = vv_jjt.at(iswp).at(i_level);
            CUSTOMREAL* v_kkr    = vv_kkr.at(iswp).at(i_level);

            CUSTOMREAL* v_fac_a  = vv_fac_a.at(iswp).at(i_level);
            CUSTOMREAL* v_fac_b  = vv_fac_b.at(iswp).at(i_level);
            CUSTOMREAL* v_fac_c  = vv_fac_c.at(iswp).at(i_level);
            CUSTOMREAL* v_fac_f  = vv_fac_f.at(iswp).at(i_level);
            CUSTOMREAL* v_T0v    = vv_T0v.at(iswp).at(i_level);
            CUSTOMREAL* v_T0r    = vv_T0r.at(iswp).at(i_level);
            CUSTOMREAL* v_T0t    = vv_T0t.at(iswp).at(i_level);
            CUSTOMREAL* v_T0p    = vv_T0p.at(iswp).at(i_level);
            CUSTOMREAL* v_fun    = vv_fun.at(iswp).at(i_level);
            CUSTOMREAL* v_change = vv_change.at(iswp).at(i_level);

            // alias for dumped index
            uint64_t* dump_ijk   = vv_i__j__k__.at(iswp).at(i_level);
            uint64_t* dump_ip1jk = vv_ip1j__k__.at(iswp).at(i_level);
            uint64_t* dump_im1jk = vv_im1j__k__.at(iswp).at(i_level);
            uint64_t* dump_ijp1k = vv_i__jp1k__.at(iswp).at(i_level);
            uint64_t* dump_ijm1k = vv_i__jm1k__.at(iswp).at(i_level);
            uint64_t* dump_ijkp1 = vv_i__j__kp1.at(iswp).at(i_level);
            uint64_t* dump_ijkm1 = vv_i__j__km1.at(iswp).at(i_level);

            // load data of all nodes in one level on temporal aligned array
            for (int _i_vec = 0; _i_vec < num_iter; _i_vec++) {
                int i_vec = _i_vec * NSIMD;

                pg = svwhilelt_b64(i_vec, n_nodes);

                v_c__    = load_mem_gen_to_mTd(pg, grid.tau_loc,  &dump_ijk[i_vec]);
                v_p__    = load_mem_gen_to_mTd(pg, grid.tau_loc,  &dump_ip1jk[i_vec]);
                v_m__    = load_mem_gen_to_mTd(pg, grid.tau_loc,  &dump_im1jk[i_vec]);
                v__p_    = load_mem_gen_to_mTd(pg, grid.tau_loc,  &dump_ijp1k[i_vec]);
                v__m_    = load_mem_gen_to_mTd(pg, grid.tau_loc,  &dump_ijm1k[i_vec]);
                v___p    = load_mem_gen_to_mTd(pg, grid.tau_loc,  &dump_ijkp1[i_vec]);
                v___m    = load_mem_gen_to_mTd(pg, grid.tau_loc,  &dump_ijkm1[i_vec]);

                // load v_iip, v_jjt, v_kkr
                v_iip_   = svld1_vnum_f64(pg, v_iip   , _i_vec);
                v_jjt_   = svld1_vnum_f64(pg, v_jjt   , _i_vec);
                v_kkr_   = svld1_vnum_f64(pg, v_kkr   , _i_vec);
                v_fac_a_ = svld1_vnum_f64(pg, v_fac_a , _i_vec);
                v_fac_b_ = svld1_vnum_f64(pg, v_fac_b , _i_vec);
                v_fac_c_ = svld1_vnum_f64(pg, v_fac_c , _i_vec);
                v_fac_f_ = svld1_vnum_f64(pg, v_fac_f , _i_vec);
                v_T0v_   = svld1_vnum_f64(pg, v_T0v   , _i_vec);
                v_T0r_   = svld1_vnum_f64(pg, v_T0r   , _i_vec);
                v_T0t_   = svld1_vnum_f64(pg, v_T0t   , _i_vec);
                v_T0p_   = svld1_vnum_f64(pg, v_T0p   , _i_vec);
                v_fun_   = svld1_vnum_f64(pg, v_fun   , _i_vec);
                v_change_= svld1_vnum_f64(pg, v_change, _i_vec);

                // loop over all nodes in one level
                vect_stencil_1st_pre_simd(pg, v_iip_, v_jjt_, v_kkr_, \
                                          v_c__, \
                                          v_p__,    v_m__,    v__p_,    v__m_,    v___p,    v___m, \
                                          v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                          v_DP_inv, v_DT_inv, v_DR_inv, \
                                          v_DP_inv_half, v_DT_inv_half, v_DR_inv_half, \
                                          loc_I, loc_J, loc_K);

                // calculate updated value on c
                vect_stencil_1st_3rd_apre_simd(pg, v_c__, v_fac_a_, v_fac_b_, v_fac_c_, v_fac_f_, \
                                               v_T0v_, v_T0p_, v_T0t_, v_T0r_, v_fun_, v_change_, \
                                               v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                               v_DP_inv, v_DT_inv, v_DR_inv);

                // store v_c__ to dump_c__
                svst1_scatter_u64index_f64(pg, grid.tau_loc, svld1_u64(pg,&dump_ijk[i_vec]), v_c__);


            } // end of i_vec loop

            // mpi synchronization
            synchronize_all_sub();

        } // end of i_level loop

#endif // ifndef USE_SIMD

    } // end of if !use_gpu
    else { // if use_gpu

#if defined USE_CUDA

        // copy tau to device
        cuda_copy_tau_to_device(gpu_grid, grid.tau_loc);

        // run iteration
        cuda_run_iteration_forward(gpu_grid, iswp);

        // copy tau to host
        cuda_copy_tau_to_host(gpu_grid, grid.tau_loc);

#else // !defiend USE_CUDA
        // exit code
        std::cout << "Error: USE_CUDA is not defined" << std::endl;
        exit(1);
#endif

    } // end of if use_gpu

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

    if(!use_gpu) {

#if !defined USE_SIMD

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

#elif USE_AVX512 || USE_AVX

        //
        __mTd v_DP_inv      = _mmT_set1_pd(1.0/dp);
        __mTd v_DT_inv      = _mmT_set1_pd(1.0/dt);
        __mTd v_DR_inv      = _mmT_set1_pd(1.0/dr);
        __mTd v_DP_inv_half = _mmT_set1_pd(1.0/dp*0.5);
        __mTd v_DT_inv_half = _mmT_set1_pd(1.0/dt*0.5);
        __mTd v_DR_inv_half = _mmT_set1_pd(1.0/dr*0.5);

        // store stencil coefs
        __mTd v_pp1;
        __mTd v_pp2;
        __mTd v_pt1;
        __mTd v_pt2;
        __mTd v_pr1;
        __mTd v_pr2;

        // measure time for only loop
        //auto start = std::chrono::high_resolution_clock::now();

        int n_levels = ijk_for_this_subproc.size();
        for (int i_level = 0; i_level < n_levels; i_level++) {
            int n_nodes = ijk_for_this_subproc.at(i_level).size();
            //std::cout << "n_nodes = " << n_nodes << std::endl;

            int num_iter = n_nodes / NSIMD + (n_nodes % NSIMD == 0 ? 0 : 1);

            // make alias to preloaded data
            __mTd* v_iip    = (__mTd*) vv_iip.at(iswp).at(i_level);
            __mTd* v_jjt    = (__mTd*) vv_jjt.at(iswp).at(i_level);
            __mTd* v_kkr    = (__mTd*) vv_kkr.at(iswp).at(i_level);

            __mTd* v_fac_a  = (__mTd*) vv_fac_a.at(iswp).at(i_level);
            __mTd* v_fac_b  = (__mTd*) vv_fac_b.at(iswp).at(i_level);
            __mTd* v_fac_c  = (__mTd*) vv_fac_c.at(iswp).at(i_level);
            __mTd* v_fac_f  = (__mTd*) vv_fac_f.at(iswp).at(i_level);
            __mTd* v_T0v    = (__mTd*) vv_T0v.at(iswp).at(i_level);
            __mTd* v_T0r    = (__mTd*) vv_T0r.at(iswp).at(i_level);
            __mTd* v_T0t    = (__mTd*) vv_T0t.at(iswp).at(i_level);
            __mTd* v_T0p    = (__mTd*) vv_T0p.at(iswp).at(i_level);
            __mTd* v_fun    = (__mTd*) vv_fun.at(iswp).at(i_level);
            __mTd* v_change = (__mTd*) vv_change.at(iswp).at(i_level);

            // alias for dumped index
            int* dump_ijk   = vv_i__j__k__.at(iswp).at(i_level);
            int* dump_ip1jk = vv_ip1j__k__.at(iswp).at(i_level);
            int* dump_im1jk = vv_im1j__k__.at(iswp).at(i_level);
            int* dump_ijp1k = vv_i__jp1k__.at(iswp).at(i_level);
            int* dump_ijm1k = vv_i__jm1k__.at(iswp).at(i_level);
            int* dump_ijkp1 = vv_i__j__kp1.at(iswp).at(i_level);
            int* dump_ijkm1 = vv_i__j__km1.at(iswp).at(i_level);
            int* dump_ip2jk = vv_ip2j__k__.at(iswp).at(i_level);
            int* dump_im2jk = vv_im2j__k__.at(iswp).at(i_level);
            int* dump_ijp2k = vv_i__jp2k__.at(iswp).at(i_level);
            int* dump_ijm2k = vv_i__jm2k__.at(iswp).at(i_level);
            int* dump_ijkp2 = vv_i__j__kp2.at(iswp).at(i_level);
            int* dump_ijkm2 = vv_i__j__km2.at(iswp).at(i_level);

            // load data of all nodes in one level on temporal aligned array
            for (int _i_vec = 0; _i_vec < num_iter; _i_vec++) {

                int i_vec = _i_vec * NSIMD;
                __mTd v_c__    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijk  [i_vec]);
                __mTd v_p__    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ip1jk[i_vec]);
                __mTd v_m__    = load_mem_gen_to_mTd(grid.tau_loc, &dump_im1jk[i_vec]);
                __mTd v__p_    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijp1k[i_vec]);
                __mTd v__m_    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijm1k[i_vec]);
                __mTd v___p    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijkp1[i_vec]);
                __mTd v___m    = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijkm1[i_vec]);
                __mTd v_pp____ = load_mem_gen_to_mTd(grid.tau_loc, &dump_ip2jk[i_vec]);
                __mTd v_mm____ = load_mem_gen_to_mTd(grid.tau_loc, &dump_im2jk[i_vec]);
                __mTd v___pp__ = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijp2k[i_vec]);
                __mTd v___mm__ = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijm2k[i_vec]);
                __mTd v_____pp = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijkp2[i_vec]);
                __mTd v_____mm = load_mem_gen_to_mTd(grid.tau_loc, &dump_ijkm2[i_vec]);

                // loop over all nodes in one level
                vect_stencil_3rd_pre_simd(v_iip[_i_vec], v_jjt[_i_vec], v_kkr[_i_vec], \
                                          v_c__, \
                                          v_p__,    v_m__,    v__p_,    v__m_,    v___p,    v___m, \
                                          v_pp____, v_mm____, v___pp__, v___mm__, v_____pp, v_____mm, \
                                          v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                          v_DP_inv, v_DT_inv, v_DR_inv, \
                                          v_DP_inv_half, v_DT_inv_half, v_DR_inv_half, \
                                          loc_I, loc_J, loc_K);

                //// calculate updated value on c
                vect_stencil_1st_3rd_apre_simd(v_c__, v_fac_a[_i_vec], v_fac_b[_i_vec], v_fac_c[_i_vec], v_fac_f[_i_vec], \
                                               v_T0v[_i_vec], v_T0p[_i_vec]  , v_T0t[_i_vec]  , v_T0r[_i_vec]  , v_fun[_i_vec]  , v_change[_i_vec], \
                                               v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                               v_DP_inv, v_DT_inv, v_DR_inv);

                // store v_c__ to dump_c__
                _mmT_store_pd(dump_c__, v_c__);


                for (int i = 0; i < NSIMD; i++) {
                    if(i_vec+i>=n_nodes) break;

                    grid.tau_loc[dump_ijk[i_vec+i]] = dump_c__[i];
                }



            } // end of i_vec loop

            // mpi synchronization
            synchronize_all_sub();

        } // end of i_level loop

#elif USE_ARM_SVE
        //
        svbool_t pg;

        __mTd v_DP_inv      = svdup_f64(1.0/dp);
        __mTd v_DT_inv      = svdup_f64(1.0/dt);
        __mTd v_DR_inv      = svdup_f64(1.0/dr);
        __mTd v_DP_inv_half = svdup_f64(1.0/dp*0.5);
        __mTd v_DT_inv_half = svdup_f64(1.0/dt*0.5);
        __mTd v_DR_inv_half = svdup_f64(1.0/dr*0.5);

        // store stencil coefs
        __mTd v_pp1;
        __mTd v_pp2;
        __mTd v_pt1;
        __mTd v_pt2;
        __mTd v_pr1;
        __mTd v_pr2;

        __mTd v_c__   ;
        __mTd v_p__   ;
        __mTd v_m__   ;
        __mTd v__p_   ;
        __mTd v__m_   ;
        __mTd v___p   ;
        __mTd v___m   ;
        __mTd v_pp____;
        __mTd v_mm____;
        __mTd v___pp__;
        __mTd v___mm__;
        __mTd v_____pp;
        __mTd v_____mm;

        __mTd v_iip_   ;
        __mTd v_jjt_   ;
        __mTd v_kkr_   ;
        __mTd v_fac_a_ ;
        __mTd v_fac_b_ ;
        __mTd v_fac_c_ ;
        __mTd v_fac_f_ ;
        __mTd v_T0v_   ;
        __mTd v_T0r_   ;
        __mTd v_T0t_   ;
        __mTd v_T0p_   ;
        __mTd v_fun_   ;
        __mTd v_change_;


        // measure time for only loop
        //auto start = std::chrono::high_resolution_clock::now();

        int n_levels = ijk_for_this_subproc.size();
        for (int i_level = 0; i_level < n_levels; i_level++) {
            int n_nodes = ijk_for_this_subproc.at(i_level).size();
            //std::cout << "n_nodes = " << n_nodes << std::endl;

            int num_iter = n_nodes / NSIMD + (n_nodes % NSIMD == 0 ? 0 : 1);

            // make alias to preloaded data
            CUSTOMREAL* v_iip    = vv_iip.at(iswp).at(i_level);
            CUSTOMREAL* v_jjt    = vv_jjt.at(iswp).at(i_level);
            CUSTOMREAL* v_kkr    = vv_kkr.at(iswp).at(i_level);

            CUSTOMREAL* v_fac_a  = vv_fac_a.at(iswp).at(i_level);
            CUSTOMREAL* v_fac_b  = vv_fac_b.at(iswp).at(i_level);
            CUSTOMREAL* v_fac_c  = vv_fac_c.at(iswp).at(i_level);
            CUSTOMREAL* v_fac_f  = vv_fac_f.at(iswp).at(i_level);
            CUSTOMREAL* v_T0v    = vv_T0v.at(iswp).at(i_level);
            CUSTOMREAL* v_T0r    = vv_T0r.at(iswp).at(i_level);
            CUSTOMREAL* v_T0t    = vv_T0t.at(iswp).at(i_level);
            CUSTOMREAL* v_T0p    = vv_T0p.at(iswp).at(i_level);
            CUSTOMREAL* v_fun    = vv_fun.at(iswp).at(i_level);
            CUSTOMREAL* v_change = vv_change.at(iswp).at(i_level);

            // alias for dumped index
            uint64_t* dump_ijk   = vv_i__j__k__.at(iswp).at(i_level);
            uint64_t* dump_ip1jk = vv_ip1j__k__.at(iswp).at(i_level);
            uint64_t* dump_im1jk = vv_im1j__k__.at(iswp).at(i_level);
            uint64_t* dump_ijp1k = vv_i__jp1k__.at(iswp).at(i_level);
            uint64_t* dump_ijm1k = vv_i__jm1k__.at(iswp).at(i_level);
            uint64_t* dump_ijkp1 = vv_i__j__kp1.at(iswp).at(i_level);
            uint64_t* dump_ijkm1 = vv_i__j__km1.at(iswp).at(i_level);
            uint64_t* dump_ip2jk = vv_ip2j__k__.at(iswp).at(i_level);
            uint64_t* dump_im2jk = vv_im2j__k__.at(iswp).at(i_level);
            uint64_t* dump_ijp2k = vv_i__jp2k__.at(iswp).at(i_level);
            uint64_t* dump_ijm2k = vv_i__jm2k__.at(iswp).at(i_level);
            uint64_t* dump_ijkp2 = vv_i__j__kp2.at(iswp).at(i_level);
            uint64_t* dump_ijkm2 = vv_i__j__km2.at(iswp).at(i_level);

            // load data of all nodes in one level on temporal aligned array
            for (int _i_vec = 0; _i_vec < num_iter; _i_vec++) {
                int i_vec = _i_vec * NSIMD;

                pg = svwhilelt_b64(i_vec, n_nodes);

                v_c__    = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ijk  [i_vec]);
                v_p__    = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ip1jk[i_vec]);
                v_m__    = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_im1jk[i_vec]);
                v__p_    = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ijp1k[i_vec]);
                v__m_    = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ijm1k[i_vec]);
                v___p    = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ijkp1[i_vec]);
                v___m    = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ijkm1[i_vec]);
                v_pp____ = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ip2jk[i_vec]);
                v_mm____ = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_im2jk[i_vec]);
                v___pp__ = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ijp2k[i_vec]);
                v___mm__ = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ijm2k[i_vec]);
                v_____pp = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ijkp2[i_vec]);
                v_____mm = load_mem_gen_to_mTd(pg, grid.tau_loc, &dump_ijkm2[i_vec]);

                // load v_iip, v_jjt, v_kkr
                v_iip_   = svld1_vnum_f64(pg, v_iip   , _i_vec);
                v_jjt_   = svld1_vnum_f64(pg, v_jjt   , _i_vec);
                v_kkr_   = svld1_vnum_f64(pg, v_kkr   , _i_vec);
                v_fac_a_ = svld1_vnum_f64(pg, v_fac_a , _i_vec);
                v_fac_b_ = svld1_vnum_f64(pg, v_fac_b , _i_vec);
                v_fac_c_ = svld1_vnum_f64(pg, v_fac_c , _i_vec);
                v_fac_f_ = svld1_vnum_f64(pg, v_fac_f , _i_vec);
                v_T0v_   = svld1_vnum_f64(pg, v_T0v   , _i_vec);
                v_T0r_   = svld1_vnum_f64(pg, v_T0r   , _i_vec);
                v_T0t_   = svld1_vnum_f64(pg, v_T0t   , _i_vec);
                v_T0p_   = svld1_vnum_f64(pg, v_T0p   , _i_vec);
                v_fun_   = svld1_vnum_f64(pg, v_fun   , _i_vec);
                v_change_= svld1_vnum_f64(pg, v_change, _i_vec);

                // loop over all nodes in one level
                vect_stencil_3rd_pre_simd(pg, v_iip_, v_jjt_, v_kkr_, \
                                          v_c__, \
                                          v_p__,    v_m__,    v__p_,    v__m_,    v___p,    v___m, \
                                          v_pp____, v_mm____, v___pp__, v___mm__, v_____pp, v_____mm, \
                                          v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                          v_DP_inv, v_DT_inv, v_DR_inv, \
                                          v_DP_inv_half, v_DT_inv_half, v_DR_inv_half, \
                                          loc_I, loc_J, loc_K);

                //// calculate updated value on c
                vect_stencil_1st_3rd_apre_simd(pg, v_c__, v_fac_a_, v_fac_b_, v_fac_c_, v_fac_f_, \
                                               v_T0v_, v_T0p_, v_T0t_, v_T0r_, v_fun_, v_change_, \
                                               v_pp1, v_pp2, v_pt1, v_pt2, v_pr1, v_pr2, \
                                               v_DP_inv, v_DT_inv, v_DR_inv);

                // store v_c__ to dump_c__
                svst1_scatter_u64index_f64(pg, grid.tau_loc, svld1_u64(pg,&dump_ijk[i_vec]), v_c__);
            } // end of i_vec loop

            // mpi synchronization
            synchronize_all_sub();

        } // end of i_level loop

#endif // ifndef USE_SIMD

    } // end of if !use_gpu
    else { // if use_gpu

#if defined USE_CUDA

        // copy tau to device
        cuda_copy_tau_to_device(gpu_grid, grid.tau_loc);

        // run iteration
        cuda_run_iteration_forward(gpu_grid, iswp);

        // copy tau to host
        cuda_copy_tau_to_host(gpu_grid, grid.tau_loc);

#else // !defiend USE_CUDA
        // exit code
        std::cout << "Error: USE_CUDA is not defined" << std::endl;
        exit(1);
#endif

    } // end of if use_gpu

    // update boundary
    if (subdom_main) {
        calculate_boundary_nodes(grid);
    }

}

// Iterator_level_1st_order_upwind::Iterator_level_1st_order_upwind(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in)
//                                 : Iterator_level(IP, grid, src, io, first_init, is_teleseismic_in) {
//     // initialization is done in the base class
// }

// void Iterator_level_1st_order_upwind::do_sweep(int iswp, Grid& grid, InputParams& IP){
//     // set sweep direction
//     set_sweep_direction(iswp);

//     int iip, jjt, kkr;

//     for (int i_level = st_level; i_level <= ed_level; i_level++) {
//         for (auto& ijk : ijk_for_this_subproc[i_level-st_level]) {
//             int kk = ijk.at(2);
//             int jj = ijk.at(1);
//             int ii = ijk.at(0);

//             if (r_dirc < 0) kkr = nr-kk; //kk-1;
//             else            kkr = kk-1;  //nr-kk;
//             if (t_dirc < 0) jjt = nt-jj; //jj-1;
//             else            jjt = jj-1;  //nt-jj;
//             if (p_dirc < 0) iip = np-ii; //ii-1;
//             else            iip = ii-1;  //np-ii;

//             //
//             // calculate stencils
//             //
//             if (grid.is_changed[I2V(iip, jjt, kkr)]) {
//                 calculate_stencil_3rd_order(grid, iip, jjt, kkr);
//             } // is_changed == true
//         } // end ijk

//         // mpi synchronization
//         synchronize_all_sub();

//     } // end loop i_level

//     // update boundary
//     if (subdom_main) {
//         calculate_boundary_nodes(grid);
//     }
// }

Iterator_level_1st_order_tele::Iterator_level_1st_order_tele(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                         : Iterator_level_tele(IP, grid, src, io, first_init, is_teleseismic_in, is_second_run_in) {
    // initialization is done in the base class
}


void Iterator_level_1st_order_tele::do_sweep(int iswp, Grid& grid, InputParams& IP){

    if(!use_gpu) {

#if !defined USE_SIMD

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

#elif USE_AVX512 || USE_AVX

#elif USE_ARM_SVE

#endif // ifndef USE_SIMD

    } // end of if !use_gpu
    else { // if use_gpu

#if defined USE_CUDA

        //// copy tau to device
        //cuda_copy_tau_to_device(gpu_grid, grid.tau_loc);

        //// run iteration
        //cuda_run_iteration_forward_tele(gpu_grid, iswp);

        //// copy tau to host
        //cuda_copy_tau_to_host(gpu_grid, grid.tau_loc);

#else // !defiend USE_CUDA
        // exit code
        std::cout << "Error: USE_CUDA is not defined" << std::endl;
        exit(1);
#endif

    } // end of if use_gpu


}


Iterator_level_3rd_order_tele::Iterator_level_3rd_order_tele(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                         : Iterator_level_tele(IP, grid, src, io, first_init, is_teleseismic_in, is_second_run_in) {
    // initialization is done in the base class
}


void Iterator_level_3rd_order_tele::do_sweep(int iswp, Grid& grid, InputParams& IP){

    if(!use_gpu) {

#if !defined USE_SIMD

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

#elif USE_AVX512 || USE_AVX

#elif USE_ARM_SVE

#endif // ifndef USE_SIMD

    } // end of if !use_gpu
    else { // if use_gpu

#if defined USE_CUDA

        //// copy tau to device
        //cuda_copy_tau_to_device(gpu_grid, grid.tau_loc);

        //// run iteration
        //cuda_run_iteration_forward_tele(gpu_grid, iswp);

        //// copy tau to host
        //cuda_copy_tau_to_host(gpu_grid, grid.tau_loc);

#else // !defiend USE_CUDA
        // exit code
        std::cout << "Error: USE_CUDA is not defined" << std::endl;
        exit(1);
#endif

    } // end of if use_gpu


}