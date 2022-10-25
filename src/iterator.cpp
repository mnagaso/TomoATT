#include "iterator.h"

Iterator::Iterator(InputParams& IP, Grid& grid, Source& src, IO_utils& io, \
                   bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
         : is_second_run(is_second_run_in) {
    if(n_subprocs > 1) {

        // share necessary values between subprocs
        np = loc_I;
        nt = loc_J;
        nr = loc_K;
        dr = grid.dr;
        dt = grid.dt;
        dp = grid.dp;

        if (first_init) {
            broadcast_i_single_sub(nr,0);
            broadcast_i_single_sub(nt,0);
            broadcast_i_single_sub(np,0);
            broadcast_cr_single_sub(dr,0);
            broadcast_cr_single_sub(dt,0);
            broadcast_cr_single_sub(dp,0);

            // initialize n_elms for subprocs
            if (!subdom_main){
                loc_I = np;
                loc_J = nt;
                loc_K = nr;
                grid.dr = dr;
                grid.dt = dt;
                grid.dp = dp;
            }
        }

        // check if teleseismic source
        is_teleseismic = is_teleseismic_in;
        broadcast_bool_single_sub(is_teleseismic,0);

        // check initialized values
        if (if_verbose){
            std::cout << "nr: " << nr << " nt: " << nt << " np: " << np << std::endl;
            std::cout << "dr: " << dr << " dt: " << dt << " dp: " << dp << std::endl;
            std::cout << "loc_I: " << loc_I << " loc_J: " << loc_J << " loc_K: " << loc_K << std::endl;
        }


    } else {
        np = loc_I;
        nt = loc_J;
        nr = loc_K;
        dr = grid.dr;
        dt = grid.dt;
        dp = grid.dp;
        is_teleseismic = is_teleseismic_in;
    }

    // set initial and end indices of level set
    st_level = 6;
    ed_level = nr+nt+np-3;

    // initialize factors etc.
    initialize_arrays(IP, grid, src);

}


Iterator::~Iterator() {

#ifdef USE_AVX
    free(dump_c__);// center of stencil
    free(dump_p__);
    free(dump_m__);
    free(dump__p_);
    free(dump__m_);
    free(dump___p);
    free(dump___m);
    // second orders
    free(dump_pp____);
    free(dump_mm____);
    free(dump___pp__);
    free(dump___mm__);
    free(dump_____pp);
    free(dump_____mm);

    // center of fac_a fac_b fac_c fac_f T0v T0r T0t T0p fun
    free(dump_fac_a );
    free(dump_fac_b );
    free(dump_fac_c );
    free(dump_fac_f );
    free(dump_T0v   );
    free(dump_T0r   );
    free(dump_T0t   );
    free(dump_T0p   );
    free(dump_fun   );
    free(dump_change);
    free(dump_iip);
    free(dump_jjt);
    free(dump_kkr);
#endif

}


void Iterator::initialize_arrays(InputParams& IP, Grid& grid, Source& src) {
    if(if_verbose && myrank == 0) std::cout << "(re) initializing arrays" << std::endl;

    if (!is_second_run) { // field initialization has already been done in the second run
        if (subdom_main) {
            if (!is_teleseismic) {
                // set initial a b c and calculate a0 b0 c0 f0
                grid.setup_factors(src);
                // calculate T0 T0r T0t T0p and initialize tau
                grid.initialize_fields(src);
            } else {
                // copy T_loc arrival time on domain's boundaries.
                grid.initialize_fields_teleseismic(src, IP.get_src_point(id_sim_src));
            }
        }
    }

    // assign processes for each sweeping level
    if (IP.get_sweep_type() == SWEEP_TYPE_LEVEL) assign_processes_for_levels();

}


void Iterator::run_iteration_forward(InputParams& IP, Grid& grid, IO_utils& io, bool& first_init) {

    if(if_verbose) stdout_by_main("--- start iteration forward. ---");

    // start timer
    std::string iter_str = "iteration_forward";
    Timer timer_iter(iter_str);

    if(!use_gpu){
        // in cpu mode
        iter_count = 0; cur_diff_L1 = HUGE_VAL; cur_diff_Linf = HUGE_VAL;

        // calculate the differcence from the true solution
        if (if_test && subdom_main) {
            if (!is_teleseismic)
                grid.calc_L1_and_Linf_error(ini_err_L1, ini_err_Linf);
            else
                grid.calc_L1_and_Linf_diff_tele(cur_diff_L1, cur_diff_Linf);

            if (myrank==0)
                std::cout << "initial err values L1, inf: " << ini_err_L1 << ", " << ini_err_Linf << std::endl;
        }

        if (subdom_main) {
            grid.calc_L1_and_Linf_diff(cur_diff_L1, cur_diff_Linf);
            if (myrank==0 && if_verbose)
                std::cout << "initial diff values L1, inf: " << cur_diff_L1 << ", " << cur_diff_Linf << std::endl;
        }

        // start iteration
        while (true) {

            // store tau for comparison
            if (subdom_main){
                if(!is_teleseismic)
                    grid.tau2tau_old();
                else
                    grid.T2tau_old();
            }

            // do sweeping for all direction
            for (int iswp = nswp-1; iswp > -1; iswp--) {
                do_sweep(iswp, grid, IP);

#ifdef FREQ_SYNC_GHOST
                // synchronize ghost cells everytime after sweeping of each direction
                if (subdom_main){
                    if (!is_teleseismic)
                        grid.send_recev_boundary_data(grid.tau_loc);
                    else
                        grid.send_recev_boundary_data(grid.T_loc);
                }
#endif
            }

#ifndef FREQ_SYNC_GHOST
            // synchronize ghost cells everytime after sweeping of all directions
            // as the same method with Detrixhe2016
            if (subdom_main){
                if (!is_teleseismic)
                    grid.send_recev_boundary_data(grid.tau_loc);
                else
                    grid.send_recev_boundary_data(grid.T_loc);
            }
#endif

            // calculate the objective function
            // if converged, break the loop
            if (subdom_main) {
                if (!is_teleseismic)
                    grid.calc_L1_and_Linf_diff(cur_diff_L1, cur_diff_Linf);
                else
                    grid.calc_L1_and_Linf_diff_tele(cur_diff_L1, cur_diff_Linf);

                if(if_test) {
                    grid.calc_L1_and_Linf_error(cur_err_L1, cur_err_Linf);
                    //std::cout << "Theoretical difference: " << cur_err_L1 << ' ' << cur_err_Linf << std::endl;
                }
            }

            // broadcast the diff values
            broadcast_cr_single_sub(cur_diff_L1, 0);
            broadcast_cr_single_sub(cur_diff_Linf, 0);
            if (if_test) {
                broadcast_cr_single_sub(cur_err_L1, 0);
                broadcast_cr_single_sub(cur_err_Linf, 0);
            }

            //if (iter_count==0)
            //    std::cout << "id_sim, sub_rank, cur_diff_L1, cur_diff_Linf: " << id_sim << ", " << sub_rank << ", " << cur_diff_L1 << ", " << cur_diff_Linf << std::endl;

            // debug store temporal T fields
            //io.write_tmp_tau_h5(grid, iter_count);

            //if (cur_diff_L1 < IP.get_conv_tol() && cur_diff_Linf < IP.get_conv_tol()) { // MNMN: let us use only L1 because Linf stop decreasing when using numbers of subdomains.
            if (cur_diff_L1 < IP.get_conv_tol()) {
                //stdout_by_main("--- iteration converged. ---");
                goto iter_end;
            } else if (IP.get_max_iter() <= iter_count) {
                stdout_by_main("--- iteration reached to the maximum number of iterations. ---");
                goto iter_end;
            } else {
                if(myrank==0 && if_verbose)
                    std::cout << "iteration " << iter_count << ": " << cur_diff_L1 << ", " << cur_diff_Linf << ", " << timer_iter.get_t_delta() << "\n";
                iter_count++;
            }
        }

    iter_end:
        if (myrank==0){
            if (if_verbose){
                std::cout << "Converged at iteration " << iter_count << ": " << cur_diff_L1 << ", " << cur_diff_Linf << std::endl;
            }
            if (if_test)
                std::cout << "errors at iteration " << iter_count << ": " << cur_err_L1 << ", " << cur_err_Linf << std::endl;
        }

        // calculate T
        // teleseismic case will update T_loc, so T = T0*tau is not necessary.
        if (subdom_main && !is_teleseismic) grid.calc_T_plus_tau();

    } else {

        // GPU mode
    #ifdef USE_CUDA
        if (!is_teleseismic){
            // transfer field to GPU grid for iteration (only need to call for the first source of first inversion step)
            if(first_init)
                grid.initialize_gpu_grid(ijk_for_this_subproc);
            else
                grid.reinitialize_gpu_grid(); // update only fac_abcf

            // run iteration
            cuda_run_iterate_forward(grid.gpu_grid, IP.get_conv_tol(), IP.get_max_iter(), IP.get_stencil_order());
            // copy result on device to host
            cuda_copy_T_loc_tau_loc_to_host(grid.T_loc, grid.tau_loc, grid.gpu_grid);
        } else {
            std::cout << "ERROR: GPU mode does not support teleseismic source." << std::endl;
            exit(1);
        }
    #endif

    }

    // check the time for iteration
    if (inter_sub_rank==0 && subdom_main) {
        timer_iter.stop_timer();

        // output time in file
        std::ofstream ofs;
        ofs.open("time.txt", std::ios::app);
        ofs << "Converged at iteration " << iter_count << ", L1 " << cur_diff_L1 << ", Linf " << cur_diff_Linf << ", total time[s] " << timer_iter.get_t() << std::endl;

    }


}


void Iterator::run_iteration_adjoint(InputParams& IP, Grid& grid, IO_utils& io) {
    if(if_verbose) stdout_by_main("--- start iteration adjoint. ---");
    iter_count = 0;
    CUSTOMREAL cur_diff_L1_dummy = HUGE_VAL;
               cur_diff_Linf = HUGE_VAL;

    // initialize delta and Tadj_loc (here we use the array tau_old instead of delta for memory saving)
    if (subdom_main)
        init_delta_and_Tadj(grid, IP);
    if(if_verbose) std::cout << "checker point 1, myrank: " << myrank << ", id_sim: " << id_sim << ", id_subdomain: " << id_subdomain
               << ", subdom_main:" << subdom_main << ", world_rank: " << world_rank << std::endl;
    // fix the boundary of the adjoint field
    if (subdom_main){
        fix_boundary_Tadj(grid);
    }

    // start timer
    std::string iter_str = "iteration_adjoint";
    Timer timer_iter(iter_str);
    if(if_verbose) stdout_by_main("--- adjoint fast sweeping ... ---");

    // start iteration
    while (true) {
        // do sweeping for all direction
        for (int iswp = 0; iswp < nswp; iswp++) {
            do_sweep_adj(iswp, grid, IP);

#ifdef FREQ_SYNC_GHOST
            // copy the values of communication nodes and ghost nodes
            if (subdom_main)
                grid.send_recev_boundary_data(grid.tau_loc);
#endif
        }

#ifndef FREQ_SYNC_GHOST
        // copy the values of communication nodes and ghost nodes
        if (subdom_main)
            grid.send_recev_boundary_data(grid.tau_loc);
#endif

        // calculate the objective function
        // if converged, break the loop
        if (subdom_main) {
            grid.calc_L1_and_Linf_diff_adj(cur_diff_L1_dummy, cur_diff_Linf);
        }

        // broadcast the diff values
        //broadcast_cr_single_sub(cur_diff_L1, 0);
        broadcast_cr_single_sub(cur_diff_Linf, 0);

        // debug store temporal T fields
        //io.write_tmp_tau_h5(grid, iter_count);

        // store tau -> Tadj
        if (subdom_main)
            grid.update_Tadj();

        if (cur_diff_Linf < IP.get_conv_tol()) {
            //stdout_by_main("--- adjoint iteration converged. ---");
            goto iter_end;
        } else if (IP.get_max_iter() <= iter_count) {
            stdout_by_main("--- adjoint iteration reached the maximum number of iterations. ---");
            goto iter_end;
        } else {
            if(myrank==0 && if_verbose)
                std::cout << "iteration adj. " << iter_count << ": " << cur_diff_Linf << std::endl;
            iter_count++;
        }


    }

iter_end:
    if (myrank==0){
        if (if_verbose)
            std::cout << "Converged at adjoint iteration " << iter_count << ": "  << cur_diff_Linf << std::endl;
        if (if_test)
            std::cout << "errors at iteration " << iter_count << ": " << cur_err_Linf << std::endl;
    }

    // check the time for iteration
    if (inter_sub_rank==0 && subdom_main) timer_iter.stop_timer();


}


void Iterator::init_delta_and_Tadj(Grid& grid, InputParams& IP) {
    if(if_verbose) std::cout << "initializing delta and Tadj" << std::endl;

    for (int k = 0; k < nr; k++) {
        for (int j = 0; j < nt; j++) {
            for (int i = 0; i < np; i++) {
                grid.tau_old_loc[I2V(i,j,k)] = _0_CR; // use tau_old_loc for delta
                grid.tau_loc[I2V(i,j,k)]     = _0_CR; // use tau_loc for Tadj_loc (later copy to Tadj_loc)
                grid.Tadj_loc[I2V(i,j,k)]    = 9999999.9;
            }
        }
    }

    // set the contributions of stations from the each station
    std::vector<SrcRec>& receivers = IP.get_rec_points(id_sim_src); // get receivers

    int DEBUG_REC_COUNT = 0;

    // loop all receivers
    for (auto& rec: receivers) {
        CUSTOMREAL delta_lon = grid.get_delta_lon();
        CUSTOMREAL delta_lat = grid.get_delta_lat();
        CUSTOMREAL delta_r   = grid.get_delta_r();

        if(!rec.is_rec_pair){        // absolute travel time
            // get positions
            CUSTOMREAL rec_lon = rec.lon*DEG2RAD;
            CUSTOMREAL rec_lat = rec.lat*DEG2RAD;
            CUSTOMREAL rec_r = depth2radius(rec.dep);

            // check if the receiver is in this subdomain
            if (grid.get_lon_min_loc() <= rec_lon && rec_lon <= grid.get_lon_max_loc()  && \
                grid.get_lat_min_loc() <= rec_lat && rec_lat <= grid.get_lat_max_loc()  && \
                grid.get_r_min_loc()   <= rec_r   && rec_r   <= grid.get_r_max_loc()   ) {

                DEBUG_REC_COUNT++;

                // descretize receiver position (LOCAL ID)
                int i_rec_loc =  std::floor((rec_lon - grid.get_lon_min_loc()) / delta_lon);
                int j_rec_loc =  std::floor((rec_lat - grid.get_lat_min_loc()) / delta_lat);
                int k_rec_loc =  std::floor((rec_r   - grid.get_r_min_loc())   / delta_r);

                // discretized receiver position
                CUSTOMREAL dis_rec_lon = grid.p_loc_1d[i_rec_loc];
                CUSTOMREAL dis_rec_lat = grid.t_loc_1d[j_rec_loc];
                CUSTOMREAL dis_rec_r   = grid.r_loc_1d[k_rec_loc];

                // relative position errors
                CUSTOMREAL e_lon = std::min(_1_CR,(rec_lon - dis_rec_lon)/delta_lon);
                CUSTOMREAL e_lat = std::min(_1_CR,(rec_lat - dis_rec_lat)/delta_lat);
                CUSTOMREAL e_r   = std::min(_1_CR,(rec_r   - dis_rec_r)  /delta_r);

                // set delta values
                grid.tau_old_loc[I2V(i_rec_loc,j_rec_loc,k_rec_loc)]       += rec.t_adj*(1.0-e_lon)*(1.0-e_lat)*(1.0-e_r)/(delta_lon*delta_lat*delta_r*my_square(rec_r)*std::cos(rec_lat));
                grid.tau_old_loc[I2V(i_rec_loc,j_rec_loc,k_rec_loc+1)]     += rec.t_adj*(1.0-e_lon)*(1.0-e_lat)*     e_r /(delta_lon*delta_lat*delta_r*my_square(rec_r)*std::cos(rec_lat));
                grid.tau_old_loc[I2V(i_rec_loc,j_rec_loc+1,k_rec_loc)]     += rec.t_adj*(1.0-e_lon)*     e_lat* (1.0-e_r)/(delta_lon*delta_lat*delta_r*my_square(rec_r)*std::cos(rec_lat));
                grid.tau_old_loc[I2V(i_rec_loc,j_rec_loc+1,k_rec_loc+1)]   += rec.t_adj*(1.0-e_lon)*     e_lat*      e_r /(delta_lon*delta_lat*delta_r*my_square(rec_r)*std::cos(rec_lat));
                grid.tau_old_loc[I2V(i_rec_loc+1,j_rec_loc,k_rec_loc)]     += rec.t_adj*     e_lon *(1.0-e_lat)*(1.0-e_r)/(delta_lon*delta_lat*delta_r*my_square(rec_r)*std::cos(rec_lat));
                grid.tau_old_loc[I2V(i_rec_loc+1,j_rec_loc,k_rec_loc+1)]   += rec.t_adj*     e_lon *(1.0-e_lat)*     e_r /(delta_lon*delta_lat*delta_r*my_square(rec_r)*std::cos(rec_lat));
                grid.tau_old_loc[I2V(i_rec_loc+1,j_rec_loc+1,k_rec_loc)]   += rec.t_adj*     e_lon *     e_lat *(1.0-e_r)/(delta_lon*delta_lat*delta_r*my_square(rec_r)*std::cos(rec_lat));
                grid.tau_old_loc[I2V(i_rec_loc+1,j_rec_loc+1,k_rec_loc+1)] += rec.t_adj*     e_lon *     e_lat *     e_r /(delta_lon*delta_lat*delta_r*my_square(rec_r)*std::cos(rec_lat));
            }
        } else {        // differential travel time

            for(int index_rec = 0; index_rec < 2; index_rec++){
                // get positions
                CUSTOMREAL rec_lon = rec.lon_pair[index_rec]*DEG2RAD;
                CUSTOMREAL rec_lat = rec.lat_pair[index_rec]*DEG2RAD;
                CUSTOMREAL rec_r = depth2radius(rec.dep_pair[index_rec]);
                CUSTOMREAL rec_adj = rec.ddt_adj_pair[index_rec];

                // check if the receiver is in this subdomain
                if (grid.get_lon_min_loc() <= rec_lon && rec_lon <= grid.get_lon_max_loc()  && \
                    grid.get_lat_min_loc() <= rec_lat && rec_lat <= grid.get_lat_max_loc()  && \
                    grid.get_r_min_loc()   <= rec_r   && rec_r   <= grid.get_r_max_loc()   ) {

                    DEBUG_REC_COUNT++;

                    // descretize receiver position (LOCAL ID)
                    int i_rec_loc =  std::floor((rec_lon - grid.get_lon_min_loc()) / delta_lon);
                    int j_rec_loc =  std::floor((rec_lat - grid.get_lat_min_loc()) / delta_lat);
                    int k_rec_loc =  std::floor((rec_r   - grid.get_r_min_loc())   / delta_r);

                    // discretized receiver position
                    CUSTOMREAL dis_rec_lon = grid.p_loc_1d[i_rec_loc];
                    CUSTOMREAL dis_rec_lat = grid.t_loc_1d[j_rec_loc];
                    CUSTOMREAL dis_rec_r   = grid.r_loc_1d[k_rec_loc];

                    // relative position errors
                    CUSTOMREAL e_lon = std::min(_1_CR,(rec_lon - dis_rec_lon)/delta_lon);
                    CUSTOMREAL e_lat = std::min(_1_CR,(rec_lat - dis_rec_lat)/delta_lat);
                    CUSTOMREAL e_r   = std::min(_1_CR,(rec_r   - dis_rec_r)  /delta_r);

                    // set delta values
                    grid.tau_old_loc[I2V(i_rec_loc,j_rec_loc,k_rec_loc)]       += rec_adj*(1.0-e_lon)*(1.0-e_lat)*(1.0-e_r)/(delta_lon*delta_lat*delta_r*my_square(rec_r)*std::cos(rec_lat));
                    grid.tau_old_loc[I2V(i_rec_loc,j_rec_loc,k_rec_loc+1)]     += rec_adj*(1.0-e_lon)*(1.0-e_lat)*     e_r /(delta_lon*delta_lat*delta_r*my_square(rec_r)*std::cos(rec_lat));
                    grid.tau_old_loc[I2V(i_rec_loc,j_rec_loc+1,k_rec_loc)]     += rec_adj*(1.0-e_lon)*     e_lat* (1.0-e_r)/(delta_lon*delta_lat*delta_r*my_square(rec_r)*std::cos(rec_lat));
                    grid.tau_old_loc[I2V(i_rec_loc,j_rec_loc+1,k_rec_loc+1)]   += rec_adj*(1.0-e_lon)*     e_lat*      e_r /(delta_lon*delta_lat*delta_r*my_square(rec_r)*std::cos(rec_lat));
                    grid.tau_old_loc[I2V(i_rec_loc+1,j_rec_loc,k_rec_loc)]     += rec_adj*     e_lon *(1.0-e_lat)*(1.0-e_r)/(delta_lon*delta_lat*delta_r*my_square(rec_r)*std::cos(rec_lat));
                    grid.tau_old_loc[I2V(i_rec_loc+1,j_rec_loc,k_rec_loc+1)]   += rec_adj*     e_lon *(1.0-e_lat)*     e_r /(delta_lon*delta_lat*delta_r*my_square(rec_r)*std::cos(rec_lat));
                    grid.tau_old_loc[I2V(i_rec_loc+1,j_rec_loc+1,k_rec_loc)]   += rec_adj*     e_lon *     e_lat *(1.0-e_r)/(delta_lon*delta_lat*delta_r*my_square(rec_r)*std::cos(rec_lat));
                    grid.tau_old_loc[I2V(i_rec_loc+1,j_rec_loc+1,k_rec_loc+1)] += rec_adj*     e_lon *     e_lat *     e_r /(delta_lon*delta_lat*delta_r*my_square(rec_r)*std::cos(rec_lat));

                }
            }

        }

    }

    // communicate tau_old_loc to all processors
    grid.send_recev_boundary_data(grid.tau_old_loc);

    if(if_verbose) std::cout << "DEBUG_REC_COUNT: " << DEBUG_REC_COUNT << std::endl;

}


void Iterator::fix_boundary_Tadj(Grid& grid) {

    if (!is_teleseismic){
        // r, theta boundary
        if (grid.i_first())
            for (int ir = 0; ir < nr; ir++)
                for (int it = 0; it < nt; it++)
                    grid.tau_loc[I2V(0,it,ir)]    = _0_CR;
        if (grid.i_last())
            for (int ir = 0; ir < nr; ir++)
                for (int it = 0; it < nt; it++)
                    grid.tau_loc[I2V(np-1,it,ir)] = _0_CR;
        // r, phi boundary
        if (grid.j_first())
            for (int ir = 0; ir < nr; ir++)
                for (int ip = 0; ip < np; ip++)
                    grid.tau_loc[I2V(ip,0,ir)]    = _0_CR;
        if (grid.j_last())
            for (int ir = 0; ir < nr; ir++)
                for (int ip = 0; ip < np; ip++)
                    grid.tau_loc[I2V(ip,nt-1,ir)] = _0_CR;
        // theta, phi boundary
        if (grid.k_first())
            for (int it = 0; it < nt; it++)
                for (int ip = 0; ip < np; ip++)
                    grid.tau_loc[I2V(ip,it,0)]    = _0_CR;
        if (grid.k_last())
            for (int it = 0; it < nt; it++)
                for (int ip = 0; ip < np; ip++)
                    grid.tau_loc[I2V(ip,it,nr-1)] = _0_CR;
    } else {
        for (int ir = 0; ir < nr; ir++)
            for (int it = 0; it < nt; it++)
                for (int ip = 0; ip < np; ip++)
                    calculate_boundary_nodes_tele_adj(grid,ir,it,ip);
    }

}


// assign intra-node processes for each sweeping level
void Iterator::assign_processes_for_levels() {
    // allocate memory for process range
    std::vector<int> n_nodes_of_levels;

    // teleseismic case need to iterate outermost layer
    int st_id, st_id2;
    if (!is_teleseismic){
        st_id  = 2;
        st_id2 = 1;
    } else {
        st_id  = 1;
        st_id2 = 0;
    }

    for (int level = st_level; level <= ed_level; level++) {
        int kleft  = std::max(st_id, level-np-nt+2);
        int kright = std::min(level-4, nr-1)   ;

        int count_n_nodes = 0;

        for (int kk = kleft; kk <= kright; kk++) {
            int jleft  = std::max(st_id, level-kk-np+st_id2);
            int jright = std::min(level-kk-2, nt-1);

            int n_jlines = jright - jleft + 1;
            count_n_nodes += n_jlines;
        } // end loop kk

        n_nodes_of_levels.push_back(count_n_nodes); // store the number of nodes of each level
    } // end loop level

    // find max in n_nodes_of_levels
    max_n_nodes_plane = 0;
    for (size_t i = 0; i < n_nodes_of_levels.size(); i++)
        if (n_nodes_of_levels[i] > max_n_nodes_plane)
            max_n_nodes_plane = n_nodes_of_levels[i];


    int n_grids_this_subproc = 0; // count the number of grids calculated by this subproc

    // assign nodes on each level plane to processes
    for (int level = st_level; level <= ed_level; level++) {
        int kleft  = std::max(st_id, level-np-nt+2);
        int kright = std::min(level-4, nr-1)   ;

        std::vector<int> asigned_nodes_on_this_level;
        int grid_count = 0;

        // n grids calculated by each subproc
        int n_grids_each = static_cast<int>(n_nodes_of_levels[level-st_level] / n_subprocs);
        int n_grids_by_this = n_grids_each;
        int i_grid_start = n_grids_each*sub_rank;

        // add modulo for last sub_rank
        if (sub_rank == n_subprocs-1)
            n_grids_by_this += static_cast<int>(n_nodes_of_levels[level-st_level] % n_subprocs);


        for (int kk = kleft; kk <= kright; kk++) {
            int jleft  = std::max(st_id, level-kk-np+st_id2);
            int jright = std::min(level-kk-2, nt-1);

            for (int jj = jleft; jj <= jright; jj++) {
                int ii = level - kk - jj;

                // check if this node should be assigned to this process
                //if (grid_count%n_subprocs == sub_rank) {
                if (grid_count >= i_grid_start && grid_count < i_grid_start+n_grids_by_this) {
                    int tmp_ijk = I2V(ii,jj,kk);
                    asigned_nodes_on_this_level.push_back(tmp_ijk);
                    n_grids_this_subproc++;
                }

                grid_count++;

            } // end loop jj
        } // end loop kk

        // store the node ids of each level
        ijk_for_this_subproc.push_back(asigned_nodes_on_this_level);
    } // end loop level

    if(if_verbose)
        std::cout << "n total grids calculated by sub_rank " << sub_rank << ": " << n_grids_this_subproc << std::endl;

#ifdef USE_AVX
    // prepare node stack for all sweeping directions
    int iip, jjt, kkr;
    for (int iswp=0; iswp < 8; iswp++){
        set_sweep_direction(iswp);

        std::vector<std::vector<std::vector<int>>> ijk_for_this_subproc_optim_tmp;
        int n_levels = ijk_for_this_subproc.size();
        for (int i_level = 0; i_level < n_levels; i_level++) {
            int n_nodes = ijk_for_this_subproc[i_level].size();
            std::vector<std::vector<int>> ijk_for_this_subproc_optim_tmp_tmp;

            for (int i_node = 0; i_node < n_nodes; i_node++) {
                int tmp_ijk = ijk_for_this_subproc[i_level][i_node];
                V2I(tmp_ijk, iip, jjt, kkr);
                if (r_dirc < 0) kkr = loc_K-kkr; //kk-1;
                else            kkr = kkr-1;  //nr-kk;
                if (t_dirc < 0) jjt = loc_J-jjt; //jj-1;
                else            jjt = jjt-1;  //nt-jj;
                if (p_dirc < 0) iip = loc_I-iip; //ii-1;
                else            iip = iip-1;  //np-ii;
                //tmp_ijk = I2V(iip, jjt, kkr);
                std::vector<int> tmp_ijk_vec = {iip, jjt, kkr};

                ijk_for_this_subproc_optim_tmp_tmp.push_back(tmp_ijk_vec);
            }
            ijk_for_this_subproc_optim_tmp.push_back(ijk_for_this_subproc_optim_tmp_tmp);
        }

        ijk_for_this_subproc_swps.push_back(ijk_for_this_subproc_optim_tmp);
    }

    int NSIMD = 4;
    int pad_length = 0;
    int max_length = max_n_nodes_plane;
    if (max_length % NSIMD != 0)
        pad_length = NSIMD - max_length % NSIMD;
//    int dump_length = (max_length + pad_length);
    int dump_length = NSIMD;


    // stencil dumps
    // first orders
    const int ALIGN     = 32;//sizeof(CUSTOMREAL); // double: 64, float: 32
    //const int ALIGN_int = 32;//sizeof(int);        // int: 32

    dump_c__ = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));// center of stencil
    dump_p__ = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump_m__ = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump__p_ = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump__m_ = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump___p = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump___m = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    // second orders
    dump_pp____ = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump_mm____ = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump___pp__ = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump___mm__ = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump_____pp = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump_____mm = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));

    // center of fac_a fac_b fac_c fac_f T0v T0r T0t T0p fun
    dump_fac_a  = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump_fac_b  = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump_fac_c  = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump_fac_f  = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump_T0v    = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump_T0r    = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump_T0t    = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump_T0p    = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump_fun    = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump_change = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));

//    dump_iip = (int*) aligned_alloc(ALIGN_int, dump_length*sizeof(int));
//    dump_jjt = (int*) aligned_alloc(ALIGN_int, dump_length*sizeof(int));
//    dump_kkr = (int*) aligned_alloc(ALIGN_int, dump_length*sizeof(int));
    dump_iip = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump_jjt = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));
    dump_kkr = (CUSTOMREAL*) aligned_alloc(ALIGN, dump_length*sizeof(CUSTOMREAL));

    dump_icc = (int*) aligned_alloc(ALIGN, dump_length*sizeof(int));
    dump_jcc = (int*) aligned_alloc(ALIGN, dump_length*sizeof(int));
    dump_kcc = (int*) aligned_alloc(ALIGN, dump_length*sizeof(int));

    dump_ip1 = (int*) aligned_alloc(ALIGN, dump_length*sizeof(int));
    dump_jp1 = (int*) aligned_alloc(ALIGN, dump_length*sizeof(int));
    dump_kp1 = (int*) aligned_alloc(ALIGN, dump_length*sizeof(int));

    dump_im1 = (int*) aligned_alloc(ALIGN, dump_length*sizeof(int));
    dump_jm1 = (int*) aligned_alloc(ALIGN, dump_length*sizeof(int));
    dump_km1 = (int*) aligned_alloc(ALIGN, dump_length*sizeof(int));

    dump_ip2 = (int*) aligned_alloc(ALIGN, dump_length*sizeof(int));
    dump_jp2 = (int*) aligned_alloc(ALIGN, dump_length*sizeof(int));
    dump_kp2 = (int*) aligned_alloc(ALIGN, dump_length*sizeof(int));

    dump_im2 = (int*) aligned_alloc(ALIGN, dump_length*sizeof(int));
    dump_jm2 = (int*) aligned_alloc(ALIGN, dump_length*sizeof(int));
    dump_km2 = (int*) aligned_alloc(ALIGN, dump_length*sizeof(int));



#endif


}


void Iterator::calculate_stencil_1st_order(Grid& grid, int& iip, int& jjt, int&kkr){
    sigr = SWEEPING_COEFF*std::sqrt(grid.fac_a_loc[I2V(iip, jjt, kkr)])*grid.T0v_loc[I2V(iip, jjt, kkr)];
    sigt = SWEEPING_COEFF*std::sqrt(grid.fac_b_loc[I2V(iip, jjt, kkr)])*grid.T0v_loc[I2V(iip, jjt, kkr)];
    sigp = SWEEPING_COEFF*std::sqrt(grid.fac_c_loc[I2V(iip, jjt, kkr)])*grid.T0v_loc[I2V(iip, jjt, kkr)];
    coe  = _1_CR/((sigr/dr)+(sigt/dt)+(sigp/dp));

    pp1 = (grid.tau_loc[I2V(iip  , jjt  , kkr  )] - grid.tau_loc[I2V(iip-1, jjt  , kkr  )])/dp;
    pp2 = (grid.tau_loc[I2V(iip+1, jjt  , kkr  )] - grid.tau_loc[I2V(iip  , jjt  , kkr  )])/dp;

    pt1 = (grid.tau_loc[I2V(iip  , jjt  , kkr  )] - grid.tau_loc[I2V(iip  , jjt-1, kkr  )])/dt;
    pt2 = (grid.tau_loc[I2V(iip  , jjt+1, kkr  )] - grid.tau_loc[I2V(iip  , jjt  , kkr  )])/dt;

    pr1 = (grid.tau_loc[I2V(iip  , jjt  , kkr  )] - grid.tau_loc[I2V(iip  , jjt  , kkr-1)])/dr;
    pr2 = (grid.tau_loc[I2V(iip  , jjt  , kkr+1)] - grid.tau_loc[I2V(iip  , jjt  , kkr  )])/dr;

    // LF Hamiltonian
    Htau = calc_LF_Hamiltonian(grid, pp1, pp2, pt1, pt2, pr1, pr2, iip, jjt, kkr);

    grid.tau_loc[I2V(iip, jjt, kkr)] += coe*(grid.fun_loc[I2V(iip, jjt, kkr)] - Htau) \
                                      + coe*(sigr*(pr2-pr1)/_2_CR + sigt*(pt2-pt1)/_2_CR + sigp*(pp2-pp1)/_2_CR);

}


void Iterator::calculate_stencil_3rd_order(Grid& grid, int& iip, int& jjt, int&kkr){
    sigr = SWEEPING_COEFF*std::sqrt(grid.fac_a_loc[I2V(iip, jjt, kkr)])*grid.T0v_loc[I2V(iip, jjt, kkr)];
    sigt = SWEEPING_COEFF*std::sqrt(grid.fac_b_loc[I2V(iip, jjt, kkr)])*grid.T0v_loc[I2V(iip, jjt, kkr)];
    sigp = SWEEPING_COEFF*std::sqrt(grid.fac_c_loc[I2V(iip, jjt, kkr)])*grid.T0v_loc[I2V(iip, jjt, kkr)];
    coe  = _1_CR/((sigr/dr)+(sigt/dt)+(sigp/dp));


    // direction p
    if (iip == 1) {
        pp1 = (grid.tau_loc[I2V(iip,jjt,kkr)] \
             - grid.tau_loc[I2V(iip-1,jjt,kkr)]) / dp;

        wp2 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.tau_loc[I2V(iip,  jjt,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip+1,jjt,kkr)] \
                                                           +      grid.tau_loc[I2V(iip+2,jjt,kkr)]) ) \

                                         / (eps + my_square(      grid.tau_loc[I2V(iip-1,jjt,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,  jjt,kkr)] \
                                                           +      grid.tau_loc[I2V(iip+1,jjt,kkr)]) )) );

        pp2 = (_1_CR - wp2) * (         grid.tau_loc[I2V(iip+1,jjt,kkr)] \
                              -         grid.tau_loc[I2V(iip-1,jjt,kkr)]) / _2_CR / dp \
                   + wp2  * ( - _3_CR * grid.tau_loc[I2V(iip,jjt,kkr)] \
                              + _4_CR * grid.tau_loc[I2V(iip+1,jjt,kkr)] \
                              -         grid.tau_loc[I2V(iip+2,jjt,kkr)] ) / _2_CR / dp;

    } else if (iip == np-2) {
        wp1 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.tau_loc[I2V(iip,  jjt,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip-1,jjt,kkr)] \
                                                           +      grid.tau_loc[I2V(iip-2,jjt,kkr)]) ) \

                                         / (eps + my_square(      grid.tau_loc[I2V(iip+1,jjt,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,  jjt,kkr)] \
                                                           +      grid.tau_loc[I2V(iip-1,jjt,kkr)]) )) );

        pp1 = (_1_CR - wp1) * (           grid.tau_loc[I2V(iip+1,jjt,kkr)] \
                                -         grid.tau_loc[I2V(iip-1,jjt,kkr)]) / _2_CR / dp \
                     + wp1  * ( + _3_CR * grid.tau_loc[I2V(iip,jjt,kkr)] \
                                - _4_CR * grid.tau_loc[I2V(iip-1,jjt,kkr)] \
                                +         grid.tau_loc[I2V(iip-2,jjt,kkr)] ) / _2_CR / dp;

        pp2 = (grid.tau_loc[I2V(iip+1,jjt,kkr)] - grid.tau_loc[I2V(iip,jjt,kkr)]) / dp;

    } else {
        wp1 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.tau_loc[I2V(iip,  jjt,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip-1,jjt,kkr)] \
                                                           +      grid.tau_loc[I2V(iip-2,jjt,kkr)]) ) \

                                         / (eps + my_square(      grid.tau_loc[I2V(iip+1,jjt,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,  jjt,kkr)] \
                                                           +      grid.tau_loc[I2V(iip-1,jjt,kkr)]) )) );

        pp1 = (_1_CR - wp1) * (         grid.tau_loc[I2V(iip+1,jjt,kkr)] \
                              -         grid.tau_loc[I2V(iip-1,jjt,kkr)]) / _2_CR / dp \
                   + wp1  * ( + _3_CR * grid.tau_loc[I2V(iip,jjt,kkr)] \
                              - _4_CR * grid.tau_loc[I2V(iip-1,jjt,kkr)] \
                              +         grid.tau_loc[I2V(iip-2,jjt,kkr)] ) / _2_CR / dp;

        wp2 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.tau_loc[I2V(iip,  jjt,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip+1,jjt,kkr)] \
                                                           +      grid.tau_loc[I2V(iip+2,jjt,kkr)]) ) \

                                         / (eps + my_square(      grid.tau_loc[I2V(iip-1,jjt,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,  jjt,kkr)] \
                                                           +      grid.tau_loc[I2V(iip+1,jjt,kkr)]) )) );

        pp2 = (_1_CR - wp2) * (           grid.tau_loc[I2V(iip+1,jjt,kkr)] \
                                -         grid.tau_loc[I2V(iip-1,jjt,kkr)]) / _2_CR / dp \
                     + wp2  * ( - _3_CR * grid.tau_loc[I2V(iip  ,jjt,kkr)] \
                                + _4_CR * grid.tau_loc[I2V(iip+1,jjt,kkr)] \
                                -         grid.tau_loc[I2V(iip+2,jjt,kkr)] ) / _2_CR / dp;

    }

    // direction t
    if (jjt == 1) {
        pt1 = (grid.tau_loc[I2V(iip,jjt  ,kkr)] \
             - grid.tau_loc[I2V(iip,jjt-1,kkr)]) / dt;

        wt2 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.tau_loc[I2V(iip,  jjt,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,jjt+1,kkr)] \
                                                           +      grid.tau_loc[I2V(iip,jjt+2,kkr)]) ) \

                                         / (eps + my_square(      grid.tau_loc[I2V(iip,jjt-1,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,  jjt,kkr)] \
                                                           +      grid.tau_loc[I2V(iip,jjt+1,kkr)]) )) );

        pt2 = (_1_CR - wt2) * (         grid.tau_loc[I2V(iip,jjt+1,kkr)] \
                              -         grid.tau_loc[I2V(iip,jjt-1,kkr)]) / _2_CR / dt \
                   + wt2  * ( - _3_CR * grid.tau_loc[I2V(iip,jjt,kkr)] \
                              + _4_CR * grid.tau_loc[I2V(iip,jjt+1,kkr)] \
                              -         grid.tau_loc[I2V(iip,jjt+2,kkr)] ) / _2_CR / dt;

    } else if (jjt == nt-2) {
        wt1 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.tau_loc[I2V(iip,  jjt,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,jjt-1,kkr)] \
                                                           +      grid.tau_loc[I2V(iip,jjt-2,kkr)]) ) \

                                         / (eps + my_square(       grid.tau_loc[I2V(iip,jjt+1,kkr)] \
                                                            -_2_CR*grid.tau_loc[I2V(iip,  jjt,kkr)] \
                                                            +      grid.tau_loc[I2V(iip,jjt-1,kkr)]) )) );

        pt1 = (_1_CR - wt1) * (           grid.tau_loc[I2V(iip,jjt+1,kkr)] \
                                -         grid.tau_loc[I2V(iip,jjt-1,kkr)]) / _2_CR / dt \
                     + wt1  * ( + _3_CR * grid.tau_loc[I2V(iip,jjt,kkr)] \
                                - _4_CR * grid.tau_loc[I2V(iip,jjt-1,kkr)] \
                                +         grid.tau_loc[I2V(iip,jjt-2,kkr)] ) / _2_CR / dt;

        pt2 = (grid.tau_loc[I2V(iip,jjt+1,kkr)] - grid.tau_loc[I2V(iip,jjt,kkr)]) / dt;

    } else {
        wt1 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.tau_loc[I2V(iip,  jjt,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,jjt-1,kkr)] \
                                                           +      grid.tau_loc[I2V(iip,jjt-2,kkr)]) ) \

                                         / (eps + my_square(      grid.tau_loc[I2V(iip,jjt+1,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,  jjt,kkr)] \
                                                           +      grid.tau_loc[I2V(iip,jjt-1,kkr)]) )) );

        pt1 = (_1_CR - wt1) * (           grid.tau_loc[I2V(iip,jjt+1,kkr)] \
                                -         grid.tau_loc[I2V(iip,jjt-1,kkr)]) / _2_CR / dt \
                     + wt1  * ( + _3_CR * grid.tau_loc[I2V(iip,jjt,kkr)] \
                                - _4_CR * grid.tau_loc[I2V(iip,jjt-1,kkr)] \
                                +         grid.tau_loc[I2V(iip,jjt-2,kkr)] ) / _2_CR / dt;

        wt2 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.tau_loc[I2V(iip,  jjt,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,jjt+1,kkr)] \
                                                           +      grid.tau_loc[I2V(iip,jjt+2,kkr)]) ) \

                                         / (eps + my_square(      grid.tau_loc[I2V(iip,jjt-1,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,  jjt,kkr)] \
                                                           +      grid.tau_loc[I2V(iip,jjt+1,kkr)]) )) );

        pt2 = (_1_CR - wt2) * (           grid.tau_loc[I2V(iip,jjt+1,kkr)] \
                                -         grid.tau_loc[I2V(iip,jjt-1,kkr)]) / _2_CR / dt \
                     + wt2  * ( - _3_CR * grid.tau_loc[I2V(iip,jjt,kkr)] \
                                + _4_CR * grid.tau_loc[I2V(iip,jjt+1,kkr)] \
                                -         grid.tau_loc[I2V(iip,jjt+2,kkr)] ) / _2_CR / dt;

    }

    // direction r
    if (kkr == 1) {
        pr1 = (grid.tau_loc[I2V(iip,jjt,kkr  )] \
             - grid.tau_loc[I2V(iip,jjt,kkr-1)]) / dr;

        wr2 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.tau_loc[I2V(iip,jjt,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,jjt,kkr+1)] \
                                                           +      grid.tau_loc[I2V(iip,jjt,kkr+2)]) ) \

                                         / (eps + my_square(      grid.tau_loc[I2V(iip,jjt,kkr-1)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,  jjt,kkr)] \
                                                           +      grid.tau_loc[I2V(iip,jjt,kkr+1)]) )) );

        pr2 = (_1_CR - wr2) * (           grid.tau_loc[I2V(iip,jjt,kkr+1)] \
                                -         grid.tau_loc[I2V(iip,jjt,kkr-1)]) / _2_CR / dr \
                     + wr2  * ( - _3_CR * grid.tau_loc[I2V(iip,jjt,kkr)] \
                                + _4_CR * grid.tau_loc[I2V(iip,jjt,kkr+1)] \
                                -         grid.tau_loc[I2V(iip,jjt,kkr+2)] ) / _2_CR / dr;

    } else if (kkr == nr - 2) {
        wr1 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.tau_loc[I2V(iip,jjt,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,jjt,kkr-1)] \
                                                           +      grid.tau_loc[I2V(iip,jjt,kkr-2)]) ) \

                                         / (eps + my_square(      grid.tau_loc[I2V(iip,jjt,kkr+1)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,jjt,kkr)] \
                                                           +      grid.tau_loc[I2V(iip,jjt,kkr-1)]) )) );

        pr1 = (_1_CR - wr1) * (            grid.tau_loc[I2V(iip,jjt,kkr+1)] \
                                -          grid.tau_loc[I2V(iip,jjt,kkr-1)]) / _2_CR / dr \
                      + wr1  * ( + _3_CR * grid.tau_loc[I2V(iip,jjt,kkr)] \
                                 - _4_CR * grid.tau_loc[I2V(iip,jjt,kkr-1)] \
                                 +         grid.tau_loc[I2V(iip,jjt,kkr-2)] ) / _2_CR / dr;

        pr2 = (grid.tau_loc[I2V(iip,jjt,kkr+1)] - grid.tau_loc[I2V(iip,jjt,kkr)]) / dr;

    } else {
        wr1 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.tau_loc[I2V(iip,jjt,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,jjt,kkr-1)] \
                                                           +      grid.tau_loc[I2V(iip,jjt,kkr-2)]) ) \

                                         / (eps + my_square(      grid.tau_loc[I2V(iip,jjt,kkr+1)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,jjt,kkr)] \
                                                           +      grid.tau_loc[I2V(iip,jjt,kkr-1)]) )) );

        pr1 = (_1_CR - wr1) * (           grid.tau_loc[I2V(iip,jjt,kkr+1)] \
                                -         grid.tau_loc[I2V(iip,jjt,kkr-1)]) / _2_CR / dr \
                     + wr1  * ( + _3_CR * grid.tau_loc[I2V(iip,jjt,kkr)] \
                                - _4_CR * grid.tau_loc[I2V(iip,jjt,kkr-1)] \
                                +         grid.tau_loc[I2V(iip,jjt,kkr-2)] ) / _2_CR / dr;

        wr2 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.tau_loc[I2V(iip,jjt,kkr)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,jjt,kkr+1)] \
                                                           +      grid.tau_loc[I2V(iip,jjt,kkr+2)]) ) \

                                         / (eps + my_square(      grid.tau_loc[I2V(iip,jjt,kkr-1)] \
                                                           -_2_CR*grid.tau_loc[I2V(iip,  jjt,kkr)] \
                                                           +      grid.tau_loc[I2V(iip,jjt,kkr+1)]) )) );

        pr2 = (_1_CR - wr2) * (           grid.tau_loc[I2V(iip,jjt,kkr+1)] \
                                -         grid.tau_loc[I2V(iip,jjt,kkr-1)]) / _2_CR / dr \
                     + wr2  * ( - _3_CR * grid.tau_loc[I2V(iip,jjt,kkr)] \
                                + _4_CR * grid.tau_loc[I2V(iip,jjt,kkr+1)] \
                                -         grid.tau_loc[I2V(iip,jjt,kkr+2)] ) / _2_CR / dr;

    }

    // LF Hamiltonian
    Htau = calc_LF_Hamiltonian(grid, pp1, pp2, pt1, pt2, pr1, pr2, iip, jjt, kkr);

    // update tau
    grid.tau_loc[I2V(iip,jjt,kkr)] += coe * ((grid.fun_loc[I2V(iip,jjt,kkr)] - Htau) + (sigr*(pr2-pr1) + sigt*(pt2-pt1) + sigp*(pp2-pp1))/_2_CR);

}


void Iterator::calculate_stencil_adj(Grid& grid, int& iip, int& jjt, int& kkr){
    //CUSTOMREAL tmpr1 = (grid.r_loc_1d[kkr-1]+grid.r_loc_1d[kkr])/_2_CR;
    //CUSTOMREAL tmpr2 = (grid.r_loc_1d[kkr+1]+grid.r_loc_1d[kkr])/_2_CR;
    CUSTOMREAL tmpt1 = (grid.t_loc_1d[jjt-1]+grid.t_loc_1d[jjt])/_2_CR;
    CUSTOMREAL tmpt2 = (grid.t_loc_1d[jjt]  +grid.t_loc_1d[jjt+1])/_2_CR;

    CUSTOMREAL a1  = - (_1_CR+grid.zeta_loc[I2V(iip,jjt,kkr-1)]+grid.zeta_loc[I2V(iip,jjt,kkr)]) * (grid.T_loc[I2V(iip,jjt,kkr)]-grid.T_loc[I2V(iip,jjt,kkr-1)]) / dr;
    CUSTOMREAL a1m = (a1 - std::abs(a1))/_2_CR;
    CUSTOMREAL a1p = (a1 + std::abs(a1))/_2_CR;

    CUSTOMREAL a2  = - (_1_CR+grid.zeta_loc[I2V(iip,jjt,kkr)]+grid.zeta_loc[I2V(iip,jjt,kkr+1)]) * (grid.T_loc[I2V(iip,jjt,kkr+1)]-grid.T_loc[I2V(iip,jjt,kkr)]) / dr;
    CUSTOMREAL a2m = (a2 - std::abs(a2))/_2_CR;
    CUSTOMREAL a2p = (a2 + std::abs(a2))/_2_CR;

    CUSTOMREAL b1  = - (_1_CR-grid.xi_loc[ I2V(iip,jjt-1,kkr)]-grid.xi_loc[ I2V(iip,jjt,kkr)]) /  my_square(grid.r_loc_1d[kkr]) * (grid.T_loc[I2V(iip,jjt,kkr)]-grid.T_loc[I2V(iip,jjt-1,kkr)]) / dt \
                     - (      grid.eta_loc[I2V(iip,jjt-1,kkr)]+grid.eta_loc[I2V(iip,jjt,kkr)]) / (my_square(grid.r_loc_1d[kkr])*std::cos(tmpt1)) / (_4_CR*dp) \
                     * ((grid.T_loc[I2V(iip+1,jjt-1,kkr)]-grid.T_loc[I2V(iip-1,jjt-1,kkr)]) \
                      + (grid.T_loc[I2V(iip+1,jjt,kkr)]  -grid.T_loc[I2V(iip-1,jjt,kkr)]));
    CUSTOMREAL b1m = (b1 - std::abs(b1))/_2_CR;
    CUSTOMREAL b1p = (b1 + std::abs(b1))/_2_CR;

    CUSTOMREAL b2  = - (_1_CR-grid.xi_loc[ I2V(iip,jjt,kkr)]-grid.xi_loc[ I2V(iip,jjt+1,kkr)]) /  my_square(grid.r_loc_1d[kkr]) * (grid.T_loc[I2V(iip,jjt+1,kkr)]-grid.T_loc[I2V(iip,jjt,kkr)]) / dt \
                     - (      grid.eta_loc[I2V(iip,jjt,kkr)]+grid.eta_loc[I2V(iip,jjt+1,kkr)]) / (my_square(grid.r_loc_1d[kkr])*std::cos(tmpt2)) / (_4_CR*dp) \
                     * ((grid.T_loc[I2V(iip+1,jjt,kkr)]-grid.T_loc[I2V(iip-1,jjt,kkr)]) \
                      + (grid.T_loc[I2V(iip+1,jjt+1,kkr)]-grid.T_loc[I2V(iip-1,jjt+1,kkr)]));
    CUSTOMREAL b2m = (b2 - std::abs(b2))/_2_CR;
    CUSTOMREAL b2p = (b2 + std::abs(b2))/_2_CR;

    CUSTOMREAL c1  =  - (_1_CR+grid.xi_loc[I2V(iip-1,jjt,kkr)]+grid.xi_loc[I2V(iip,jjt,kkr)]) / (my_square(grid.r_loc_1d[kkr])*my_square(std::cos(grid.t_loc_1d[jjt]))) \
                      * (grid.T_loc[I2V(iip,jjt,kkr)]-grid.T_loc[I2V(iip-1,jjt,kkr)]) / dp \
                    - (grid.eta_loc[I2V(iip-1,jjt,kkr)]+grid.eta_loc[I2V(iip,jjt,kkr)]) / (my_square(grid.r_loc_1d[kkr])*std::cos(grid.t_loc_1d[jjt])) / (_4_CR*dt) \
                         * ((grid.T_loc[I2V(iip-1,jjt+1,kkr)]-grid.T_loc[I2V(iip-1,jjt-1,kkr)]) \
                          + (grid.T_loc[I2V(iip,  jjt+1,kkr)]-grid.T_loc[I2V(iip,  jjt-1,kkr)]));
    CUSTOMREAL c1m = (c1 - std::abs(c1))/_2_CR;
    CUSTOMREAL c1p = (c1 + std::abs(c1))/_2_CR;

    CUSTOMREAL c2  =  - (_1_CR+grid.xi_loc[I2V(iip,jjt,kkr)]+grid.xi_loc[I2V(iip+1,jjt,kkr)]) / (my_square(grid.r_loc_1d[kkr])*my_square(std::cos(grid.t_loc_1d[jjt]))) \
                      * (grid.T_loc[I2V(iip+1,jjt,kkr)]-grid.T_loc[I2V(iip,jjt,kkr)]) / dp \
                      - (grid.eta_loc[I2V(iip,jjt,kkr)]+grid.eta_loc[I2V(iip+1,jjt,kkr)]) / (my_square(grid.r_loc_1d[kkr])*std::cos(grid.t_loc_1d[jjt])) / (_4_CR*dt) \
                           * ((grid.T_loc[I2V(iip,jjt+1,kkr)]-grid.T_loc[I2V(iip,jjt-1,kkr)]) \
                            + (grid.T_loc[I2V(iip+1,jjt+1,kkr)]-grid.T_loc[I2V(iip+1,jjt-1,kkr)]));
    CUSTOMREAL c2m = (c2 - std::abs(c2))/_2_CR;
    CUSTOMREAL c2p = (c2 + std::abs(c2))/_2_CR;

    // coe
    CUSTOMREAL coe = (a2p-a1m)/dr + (b2p-b1m)/dt + (c2p-c1m)/dp;

    if (isZeroAdj(coe)) {
        grid.tau_loc[I2V(iip,jjt,kkr)] = _0_CR;
    } else {
        // hamiltonian
        CUSTOMREAL Hadj = (a1p*grid.tau_loc[I2V(iip,jjt,kkr-1)] - a2m*grid.tau_loc[I2V(iip,jjt,kkr+1)]) / dr \
                        + (b1p*grid.tau_loc[I2V(iip,jjt-1,kkr)] - b2m*grid.tau_loc[I2V(iip,jjt+1,kkr)]) / dt \
                        + (c1p*grid.tau_loc[I2V(iip-1,jjt,kkr)] - c2m*grid.tau_loc[I2V(iip+1,jjt,kkr)]) / dp;

        grid.tau_loc[I2V(iip,jjt,kkr)] = (grid.tau_old_loc[I2V(iip,jjt,kkr)] + Hadj) / coe;
    }
}


void Iterator::calculate_stencil_1st_order_tele(Grid& grid, int& iip, int& jjt, int&kkr){
    sigr = SWEEPING_COEFF_TELE*std::sqrt(grid.fac_a_loc[I2V(iip, jjt, kkr)]);
    sigt = SWEEPING_COEFF_TELE*std::sqrt(grid.fac_b_loc[I2V(iip, jjt, kkr)]);
    sigp = SWEEPING_COEFF_TELE*std::sqrt(grid.fac_c_loc[I2V(iip, jjt, kkr)]);
    coe  = _1_CR/((sigr/dr)+(sigt/dt)+(sigp/dp));

    pp1 = (grid.T_loc[I2V(iip  , jjt  , kkr  )] - grid.T_loc[I2V(iip-1, jjt  , kkr  )])/dp;
    pp2 = (grid.T_loc[I2V(iip+1, jjt  , kkr  )] - grid.T_loc[I2V(iip  , jjt  , kkr  )])/dp;
    pt1 = (grid.T_loc[I2V(iip  , jjt  , kkr  )] - grid.T_loc[I2V(iip  , jjt-1, kkr  )])/dt;
    pt2 = (grid.T_loc[I2V(iip  , jjt+1, kkr  )] - grid.T_loc[I2V(iip  , jjt  , kkr  )])/dt;
    pr1 = (grid.T_loc[I2V(iip  , jjt  , kkr  )] - grid.T_loc[I2V(iip  , jjt  , kkr-1)])/dr;
    pr2 = (grid.T_loc[I2V(iip  , jjt  , kkr+1)] - grid.T_loc[I2V(iip  , jjt  , kkr  )])/dr;

    // LF Hamiltonian
    Htau = calc_LF_Hamiltonian_tele(grid, pp1, pp2, pt1, pt2, pr1, pr2, iip, jjt, kkr);

    grid.T_loc[I2V(iip, jjt, kkr)] += coe*(grid.fun_loc[I2V(iip, jjt, kkr)] - Htau) \
                                    + coe*(sigr*(pr2-pr1)/_2_CR + sigt*(pt2-pt1)/_2_CR + sigp*(pp2-pp1)/_2_CR);

}


void Iterator::calculate_stencil_3rd_order_tele(Grid& grid, int& iip, int& jjt, int&kkr){
    sigr = SWEEPING_COEFF_TELE*std::sqrt(grid.fac_a_loc[I2V(iip, jjt, kkr)]);
    sigt = SWEEPING_COEFF_TELE*std::sqrt(grid.fac_b_loc[I2V(iip, jjt, kkr)]);
    sigp = SWEEPING_COEFF_TELE*std::sqrt(grid.fac_c_loc[I2V(iip, jjt, kkr)]);
    coe  = _1_CR/((sigr/dr)+(sigt/dt)+(sigp/dp));

    // direction p
    if (iip == 1) {
        pp1 = (grid.T_loc[I2V(iip,  jjt,kkr)] \
             - grid.T_loc[I2V(iip-1,jjt,kkr)]) / dp;

        wp2 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.T_loc[I2V(iip,  jjt,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip+1,jjt,kkr)] \
                                                           +      grid.T_loc[I2V(iip+2,jjt,kkr)]) ) \

                                         / (eps + my_square(      grid.T_loc[I2V(iip-1,jjt,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,  jjt,kkr)] \
                                                           +      grid.T_loc[I2V(iip+1,jjt,kkr)]) )) );

        pp2 = (_1_CR - wp2) * (         grid.T_loc[I2V(iip+1,jjt,kkr)] \
                              -         grid.T_loc[I2V(iip-1,jjt,kkr)]) / _2_CR / dp \
                   + wp2  * ( - _3_CR * grid.T_loc[I2V(iip,jjt,kkr)] \
                              + _4_CR * grid.T_loc[I2V(iip+1,jjt,kkr)] \
                              -         grid.T_loc[I2V(iip+2,jjt,kkr)] ) / _2_CR / dp;

    } else if (iip == np-2) {
        wp1 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.T_loc[I2V(iip,  jjt,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip-1,jjt,kkr)] \
                                                           +      grid.T_loc[I2V(iip-2,jjt,kkr)]) ) \

                                         / (eps + my_square(      grid.T_loc[I2V(iip+1,jjt,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,  jjt,kkr)] \
                                                           +      grid.T_loc[I2V(iip-1,jjt,kkr)]) )) );

        pp1 = (_1_CR - wp1) * (           grid.T_loc[I2V(iip+1,jjt,kkr)] \
                                -         grid.T_loc[I2V(iip-1,jjt,kkr)]) / _2_CR / dp \
                     + wp1  * ( + _3_CR * grid.T_loc[I2V(iip,jjt,kkr)] \
                                - _4_CR * grid.T_loc[I2V(iip-1,jjt,kkr)] \
                                +         grid.T_loc[I2V(iip-2,jjt,kkr)] ) / _2_CR / dp;

        pp2 = (grid.T_loc[I2V(iip+1,jjt,kkr)] \
             - grid.T_loc[I2V(iip,  jjt,kkr)]) / dp;

    } else {
        wp1 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.T_loc[I2V(iip,  jjt,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip-1,jjt,kkr)] \
                                                           +      grid.T_loc[I2V(iip-2,jjt,kkr)]) ) \

                                         / (eps + my_square(      grid.T_loc[I2V(iip+1,jjt,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,  jjt,kkr)] \
                                                             +    grid.T_loc[I2V(iip-1,jjt,kkr)]) )) );

        pp1 = (_1_CR - wp1) * (         grid.T_loc[I2V(iip+1,jjt,kkr)] \
                              -         grid.T_loc[I2V(iip-1,jjt,kkr)]) / _2_CR / dp \
                   + wp1  * ( + _3_CR * grid.T_loc[I2V(iip,jjt,kkr)] \
                              - _4_CR * grid.T_loc[I2V(iip-1,jjt,kkr)] \
                              +         grid.T_loc[I2V(iip-2,jjt,kkr)] ) / _2_CR / dp;

        wp2 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.T_loc[I2V(iip,  jjt,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip+1,jjt,kkr)] \
                                                           +      grid.T_loc[I2V(iip+2,jjt,kkr)]) ) \

                                         / (eps + my_square(      grid.T_loc[I2V(iip-1,jjt,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,  jjt,kkr)] \
                                                           +      grid.T_loc[I2V(iip+1,jjt,kkr)]) )) );

        pp2 = (_1_CR - wp2) * (           grid.T_loc[I2V(iip+1,jjt,kkr)] \
                                -         grid.T_loc[I2V(iip-1,jjt,kkr)]) / _2_CR / dp \
                     + wp2  * ( - _3_CR * grid.T_loc[I2V(iip  ,jjt,kkr)] \
                                + _4_CR * grid.T_loc[I2V(iip+1,jjt,kkr)] \
                                -         grid.T_loc[I2V(iip+2,jjt,kkr)] ) / _2_CR / dp;

    }

    // direction t
    if (jjt == 1) {
        pt1 = (grid.T_loc[I2V(iip,jjt  ,kkr)] \
             - grid.T_loc[I2V(iip,jjt-1,kkr)]) / dt;

        wt2 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.T_loc[I2V(iip,  jjt,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,jjt+1,kkr)] \
                                                           +      grid.T_loc[I2V(iip,jjt+2,kkr)]) ) \

                                         / (eps + my_square(      grid.T_loc[I2V(iip,jjt-1,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,  jjt,kkr)] \
                                                           +      grid.T_loc[I2V(iip,jjt+1,kkr)]) )) );

        pt2 = (_1_CR - wt2) * (         grid.T_loc[I2V(iip,jjt+1,kkr)] \
                              -         grid.T_loc[I2V(iip,jjt-1,kkr)]) / _2_CR / dt \
                   + wt2  * ( - _3_CR * grid.T_loc[I2V(iip,jjt,kkr)] \
                              + _4_CR * grid.T_loc[I2V(iip,jjt+1,kkr)] \
                              -         grid.T_loc[I2V(iip,jjt+2,kkr)] ) / _2_CR / dt;

    } else if (jjt == nt-2) {
        wt1 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.T_loc[I2V(iip,  jjt,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,jjt-1,kkr)] \
                                                           +      grid.T_loc[I2V(iip,jjt-2,kkr)]) ) \

                                         / (eps + my_square(      grid.T_loc[I2V(iip,jjt+1,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,  jjt,kkr)] \
                                                           +      grid.T_loc[I2V(iip,jjt-1,kkr)]) )) );

        pt1 = (_1_CR - wt1) * (           grid.T_loc[I2V(iip,jjt+1,kkr)] \
                                -         grid.T_loc[I2V(iip,jjt-1,kkr)]) / _2_CR / dt \
                     + wt1  * ( + _3_CR * grid.T_loc[I2V(iip,jjt,kkr)] \
                                - _4_CR * grid.T_loc[I2V(iip,jjt-1,kkr)] \
                                +         grid.T_loc[I2V(iip,jjt-2,kkr)] ) / _2_CR / dt;

        pt2 = (grid.T_loc[I2V(iip,jjt+1,kkr)] \
             - grid.T_loc[I2V(iip,jjt,  kkr)]) / dt;

    } else {
        wt1 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.T_loc[I2V(iip,  jjt,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,jjt-1,kkr)] \
                                                           +      grid.T_loc[I2V(iip,jjt-2,kkr)]) ) \

                                         / (eps + my_square(      grid.T_loc[I2V(iip,jjt+1,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,  jjt,kkr)] \
                                                           +      grid.T_loc[I2V(iip,jjt-1,kkr)]) )) );

        pt1 = (_1_CR - wt1) * (           grid.T_loc[I2V(iip,jjt+1,kkr)] \
                                -         grid.T_loc[I2V(iip,jjt-1,kkr)]) / _2_CR / dt \
                     + wt1  * ( + _3_CR * grid.T_loc[I2V(iip,jjt,kkr)] \
                                - _4_CR * grid.T_loc[I2V(iip,jjt-1,kkr)] \
                                +         grid.T_loc[I2V(iip,jjt-2,kkr)] ) / _2_CR / dt;

        wt2 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.T_loc[I2V(iip,  jjt,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,jjt+1,kkr)] \
                                                           +      grid.T_loc[I2V(iip,jjt+2,kkr)]) ) \

                                         / (eps + my_square(      grid.T_loc[I2V(iip,jjt-1,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,  jjt,kkr)] \
                                                           +      grid.T_loc[I2V(iip,jjt+1,kkr)]) )) );

        pt2 = (_1_CR - wt2) * (           grid.T_loc[I2V(iip,jjt+1,kkr)] \
                                -         grid.T_loc[I2V(iip,jjt-1,kkr)]) / _2_CR / dt \
                     + wt2  * ( - _3_CR * grid.T_loc[I2V(iip,jjt,kkr)] \
                                + _4_CR * grid.T_loc[I2V(iip,jjt+1,kkr)] \
                                -         grid.T_loc[I2V(iip,jjt+2,kkr)] ) / _2_CR / dt;

    }

    // direction r
    if (kkr == 1) {
        pr1 = (grid.T_loc[I2V(iip,jjt,kkr  )] \
             - grid.T_loc[I2V(iip,jjt,kkr-1)]) / dr;

        wr2 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.T_loc[I2V(iip,jjt,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,jjt,kkr+1)] \
                                                           +      grid.T_loc[I2V(iip,jjt,kkr+2)]) ) \

                                         / (eps + my_square(      grid.T_loc[I2V(iip,jjt,kkr-1)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,  jjt,kkr)] \
                                                           +      grid.T_loc[I2V(iip,jjt,kkr+1)]) )) );

        pr2 = (_1_CR - wr2) * (           grid.T_loc[I2V(iip,jjt,kkr+1)] \
                                -         grid.T_loc[I2V(iip,jjt,kkr-1)]) / _2_CR / dr \
                     + wr2  * ( - _3_CR * grid.T_loc[I2V(iip,jjt,kkr)] \
                                + _4_CR * grid.T_loc[I2V(iip,jjt,kkr+1)] \
                                -         grid.T_loc[I2V(iip,jjt,kkr+2)] ) / _2_CR / dr;

    } else if (kkr == nr - 2) {
        wr1 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.T_loc[I2V(iip,jjt,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,jjt,kkr-1)] \
                                                           +      grid.T_loc[I2V(iip,jjt,kkr-2)]) ) \

                                         / (eps + my_square(      grid.T_loc[I2V(iip,jjt,kkr+1)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,jjt,kkr)] \
                                                           +      grid.T_loc[I2V(iip,jjt,kkr-1)]) )) );

        pr1 = (_1_CR - wr1) * (            grid.T_loc[I2V(iip,jjt,kkr+1)] \
                                -          grid.T_loc[I2V(iip,jjt,kkr-1)]) / _2_CR / dr \
                      + wr1  * ( + _3_CR * grid.T_loc[I2V(iip,jjt,kkr)] \
                                 - _4_CR * grid.T_loc[I2V(iip,jjt,kkr-1)] \
                                 +         grid.T_loc[I2V(iip,jjt,kkr-2)] ) / _2_CR / dr;

        pr2 = (grid.T_loc[I2V(iip,jjt,kkr+1)] \
             - grid.T_loc[I2V(iip,jjt,kkr)]) / dr;

    } else {
        wr1 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.T_loc[I2V(iip,jjt,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,jjt,kkr-1)] \
                                                           +      grid.T_loc[I2V(iip,jjt,kkr-2)]) ) \

                                         / (eps + my_square(      grid.T_loc[I2V(iip,jjt,kkr+1)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,jjt,kkr)] \
                                                           +      grid.T_loc[I2V(iip,jjt,kkr-1)]) )) );

        pr1 = (_1_CR - wr1) * (           grid.T_loc[I2V(iip,jjt,kkr+1)] \
                                -         grid.T_loc[I2V(iip,jjt,kkr-1)]) / _2_CR / dr \
                     + wr1  * ( + _3_CR * grid.T_loc[I2V(iip,jjt,kkr)] \
                                - _4_CR * grid.T_loc[I2V(iip,jjt,kkr-1)] \
                                +         grid.T_loc[I2V(iip,jjt,kkr-2)] ) / _2_CR / dr;

        wr2 = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      grid.T_loc[I2V(iip,jjt,kkr)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,jjt,kkr+1)] \
                                                           +      grid.T_loc[I2V(iip,jjt,kkr+2)]) ) \

                                         / (eps + my_square(      grid.T_loc[I2V(iip,jjt,kkr-1)] \
                                                           -_2_CR*grid.T_loc[I2V(iip,  jjt,kkr)] \
                                                           +      grid.T_loc[I2V(iip,jjt,kkr+1)]) )) );

        pr2 = (_1_CR - wr2) * (           grid.T_loc[I2V(iip,jjt,kkr+1)] \
                                -         grid.T_loc[I2V(iip,jjt,kkr-1)]) / _2_CR / dr \
                     + wr2  * ( - _3_CR * grid.T_loc[I2V(iip,jjt,kkr)] \
                                + _4_CR * grid.T_loc[I2V(iip,jjt,kkr+1)] \
                                -         grid.T_loc[I2V(iip,jjt,kkr+2)] ) / _2_CR / dr;

    }

    // LF Hamiltonian
    Htau = calc_LF_Hamiltonian_tele(grid, pp1, pp2, pt1, pt2, pr1, pr2, iip, jjt, kkr);

    // update tau
    grid.T_loc[I2V(iip,jjt,kkr)] += coe * (grid.fun_loc[I2V(iip,jjt,kkr)] - Htau) \
        + coe * (sigr*(pr2-pr1)/_2_CR + sigt*(pt2-pt1)/_2_CR + sigp*(pp2-pp1)/_2_CR);

    // check value
    // if(kkr==1 && jjt==1 && iip==1){
    //     std::cout   << sigr << " " << sigt << " " << sigp << " " << coe << " "
    //                 << grid.fun_loc[I2V(iip,jjt,kkr)] << " " << Htau << " " << grid.T_loc[I2V(iip,jjt,kkr)] << " "
    //                 << std::endl;
    // }
}


inline CUSTOMREAL Iterator::calc_LF_Hamiltonian(Grid& grid, \
                                         CUSTOMREAL& pp1, CUSTOMREAL& pp2, \
                                         CUSTOMREAL& pt1, CUSTOMREAL& pt2, \
                                         CUSTOMREAL& pr1, CUSTOMREAL& pr2, \
                                         int& iip, int& jjt, int& kkr) {
    // LF Hamiltonian
    return sqrt(
              grid.fac_a_loc[I2V(iip,jjt,kkr)] * my_square(grid.T0r_loc[I2V(iip,jjt,kkr)] * grid.tau_loc[I2V(iip,jjt,kkr)] + grid.T0v_loc[I2V(iip,jjt,kkr)] * (pr1+pr2)/_2_CR) \
    +         grid.fac_b_loc[I2V(iip,jjt,kkr)] * my_square(grid.T0t_loc[I2V(iip,jjt,kkr)] * grid.tau_loc[I2V(iip,jjt,kkr)] + grid.T0v_loc[I2V(iip,jjt,kkr)] * (pt1+pt2)/_2_CR) \
    +         grid.fac_c_loc[I2V(iip,jjt,kkr)] * my_square(grid.T0p_loc[I2V(iip,jjt,kkr)] * grid.tau_loc[I2V(iip,jjt,kkr)] + grid.T0v_loc[I2V(iip,jjt,kkr)] * (pp1+pp2)/_2_CR) \
    -   _2_CR*grid.fac_f_loc[I2V(iip,jjt,kkr)] * (grid.T0t_loc[I2V(iip,jjt,kkr)] * grid.tau_loc[I2V(iip,jjt,kkr)] + grid.T0v_loc[I2V(iip,jjt,kkr)] * (pt1+pt2)/_2_CR) \
                                               * (grid.T0p_loc[I2V(iip,jjt,kkr)] * grid.tau_loc[I2V(iip,jjt,kkr)] + grid.T0v_loc[I2V(iip,jjt,kkr)] * (pp1+pp2)/_2_CR) \
    );

}


inline CUSTOMREAL Iterator::calc_LF_Hamiltonian_tele(Grid& grid, \
                                         CUSTOMREAL& pp1, CUSTOMREAL& pp2, \
                                         CUSTOMREAL& pt1, CUSTOMREAL& pt2, \
                                         CUSTOMREAL& pr1, CUSTOMREAL& pr2, \
                                         int& iip, int& jjt, int& kkr) {
    // LF Hamiltonian for teleseismic source
    return std::sqrt(
              grid.fac_a_loc[I2V(iip,jjt,kkr)] * my_square((pr1+pr2)/_2_CR) \
    +         grid.fac_b_loc[I2V(iip,jjt,kkr)] * my_square((pt1+pt2)/_2_CR) \
    +         grid.fac_c_loc[I2V(iip,jjt,kkr)] * my_square((pp1+pp2)/_2_CR) \
    -   _2_CR*grid.fac_f_loc[I2V(iip,jjt,kkr)] * ((pt1+pt2)/_2_CR) \
                                               * ((pp1+pp2)/_2_CR) \
    );

}


void Iterator::calculate_boundary_nodes(Grid& grid){
    CUSTOMREAL v0, v1;

    //plane
    for (int jjt = 0; jjt < nt; jjt++){
        #pragma omp simd
        for (int iip = 0; iip < np; iip++){
            v0 = _2_CR * grid.tau_loc[I2V(iip,jjt,1)] - grid.tau_loc[I2V(iip,jjt,2)];
            v1 = grid.tau_loc[I2V(iip,jjt,2)];
            grid.tau_loc[I2V(iip,jjt,0)] = std::max({v0,v1});
            v0 = _2_CR * grid.tau_loc[I2V(iip,jjt,nr-2)] - grid.tau_loc[I2V(iip,jjt,nr-3)];
            v1 = grid.tau_loc[I2V(iip,jjt,nr-3)];
            grid.tau_loc[I2V(iip,jjt,nr-1)] = std::max({v0,v1});
        }
    }

    for (int kkr = 0; kkr < nr; kkr++){
        #pragma omp simd
        for (int iip = 0; iip < np; iip++){
            v0 = _2_CR * grid.tau_loc[I2V(iip,1,kkr)] - grid.tau_loc[I2V(iip,2,kkr)];
            v1 = grid.tau_loc[I2V(iip,2,kkr)];
            grid.tau_loc[I2V(iip,0,kkr)] = std::max({v0,v1});
            v0 = _2_CR * grid.tau_loc[I2V(iip,nt-2,kkr)] - grid.tau_loc[I2V(iip,nt-3,kkr)];
            v1 = grid.tau_loc[I2V(iip,nt-3,kkr)];
            grid.tau_loc[I2V(iip,nt-1,kkr)] = std::max({v0,v1});
        }
    }

    for (int kkr = 0; kkr < nr; kkr++){
        #pragma omp simd
        for (int jjt = 0; jjt < nt; jjt++){
            v0 = _2_CR * grid.tau_loc[I2V(1,jjt,kkr)] - grid.tau_loc[I2V(2,jjt,kkr)];
            v1 = grid.tau_loc[I2V(2,jjt,kkr)];
            grid.tau_loc[I2V(0,jjt,kkr)] = std::max({v0,v1});
            v0 = _2_CR * grid.tau_loc[I2V(np-2,jjt,kkr)] - grid.tau_loc[I2V(np-3,jjt,kkr)];
            v1 = grid.tau_loc[I2V(np-3,jjt,kkr)];
            grid.tau_loc[I2V(np-1,jjt,kkr)] = std::max({v0,v1});
        }
    }

}


void Iterator::calculate_boundary_nodes_tele(Grid& grid, int& iip, int& jjt, int& kkr){
    CUSTOMREAL v0, v1;

    // Bottom
    if (kkr == 0 && grid.k_first())
        if (grid.is_changed[I2V(iip,jjt,0)]){
            v0 = _2_CR * grid.T_loc[I2V(iip,jjt,1)] - grid.T_loc[I2V(iip,jjt,2)];
            v1 = grid.T_loc[I2V(iip,jjt,2)];
            grid.T_loc[I2V(iip,jjt,0)] = std::max({v0,v1});
        }

    // Top
    if (kkr == nr-1 && grid.k_last())
        if (grid.is_changed[I2V(iip,jjt,nr-1)]){
            v0 = _2_CR * grid.T_loc[I2V(iip,jjt,nr-2)] - grid.T_loc[I2V(iip,jjt,nr-3)];
            v1 = grid.T_loc[I2V(iip,jjt,nr-3)];
            grid.T_loc[I2V(iip,jjt,nr-1)] = std::max({v0,v1});
        }

    // South
    if (jjt == nt-1 && grid.j_first())
        if (grid.is_changed[I2V(iip,0,kkr)]){
            v0 = _2_CR * grid.T_loc[I2V(iip,1,kkr)] - grid.T_loc[I2V(iip,2,kkr)];
            v1 = grid.T_loc[I2V(iip,2,kkr)];
            grid.T_loc[I2V(iip,0,kkr)] = std::max({v0,v1});
        }
    // North
    if (jjt == 0 && grid.j_last())
        if (grid.is_changed[I2V(iip,nt-1,kkr)]){
            v0 = _2_CR * grid.T_loc[I2V(iip,nt-2,kkr)] - grid.T_loc[I2V(iip,nt-3,kkr)];
            v1 = grid.T_loc[I2V(iip,nt-3,kkr)];
            grid.T_loc[I2V(iip,nt-1,kkr)] = std::max({v0,v1});
        }

    // West
    if (iip == 0 && grid.i_first())
        if (grid.is_changed[I2V(0,jjt,kkr)]){
            v0 = _2_CR * grid.T_loc[I2V(1,jjt,kkr)] - grid.T_loc[I2V(2,jjt,kkr)];
            v1 = grid.T_loc[I2V(2,jjt,kkr)];
            grid.T_loc[I2V(0,jjt,kkr)] = std::max({v0,v1});
        }
    // East
    if (iip == np-1 && grid.i_last())
        if (grid.is_changed[I2V(np-1,jjt,kkr)]){
            v0 = _2_CR * grid.T_loc[I2V(np-2,jjt,kkr)] - grid.T_loc[I2V(np-3,jjt,kkr)];
            v1 = grid.T_loc[I2V(np-3,jjt,kkr)];
            grid.T_loc[I2V(np-1,jjt,kkr)] = std::max({v0,v1});
        }
}


void Iterator::calculate_boundary_nodes_tele_adj(Grid& grid, int& iip, int& jjt, int& kkr){

    // West
    if (iip == 0 && grid.i_first()) {
        if (!grid.is_changed[I2V(0,jjt,kkr)]) {
            if (grid.tau_loc[I2V(2,jjt,kkr)] >= 0)
                grid.tau_loc[I2V(0,jjt,kkr)] = std::max(_0_CR, _2_CR*grid.tau_loc[I2V(1,jjt,kkr)] - grid.tau_loc[I2V(2,jjt,kkr)]);
            else
                grid.tau_loc[I2V(0,jjt,kkr)] = std::min(_0_CR, _2_CR*grid.tau_loc[I2V(1,jjt,kkr)] - grid.tau_loc[I2V(2,jjt,kkr)]);
        } else {
            grid.tau_loc[I2V(0,jjt,kkr)] = _0_CR;
        }
    }

    // East
    if (iip == np-1 && grid.i_last()) {
        if (!grid.is_changed[I2V(np-1,jjt,kkr)]) {
            if (grid.tau_loc[I2V(np-3,jjt,kkr)] >= 0)
                grid.tau_loc[I2V(np-1,jjt,kkr)] = std::max(_0_CR, _2_CR*grid.tau_loc[I2V(np-2,jjt,kkr)] - grid.tau_loc[I2V(np-3,jjt,kkr)]);
            else
                grid.tau_loc[I2V(np-1,jjt,kkr)] = std::min(_0_CR, _2_CR*grid.tau_loc[I2V(np-2,jjt,kkr)] - grid.tau_loc[I2V(np-3,jjt,kkr)]);
        } else {
            grid.tau_loc[I2V(np-1,jjt,kkr)] = _0_CR;
        }
    }

    // South
    if (jjt == 0 && grid.j_first()) {
        if (!grid.is_changed[I2V(iip,0,kkr)]) {
            if (grid.tau_loc[I2V(iip,2,kkr)] >= 0)
                grid.tau_loc[I2V(iip,0,kkr)] = std::max(_0_CR, _2_CR*grid.tau_loc[I2V(iip,1,kkr)] - grid.tau_loc[I2V(iip,2,kkr)]);
            else
                grid.tau_loc[I2V(iip,0,kkr)] = std::min(_0_CR, _2_CR*grid.tau_loc[I2V(iip,1,kkr)] - grid.tau_loc[I2V(iip,2,kkr)]);
        } else {
            grid.tau_loc[I2V(iip,0,kkr)] = _0_CR;
        }
    }

    // North
    if (jjt == nt-1 && grid.j_last()) {
        if (!grid.is_changed[I2V(iip,nt-1,kkr)]) {
            if (grid.tau_loc[I2V(iip,nt-3,kkr)] >= 0)
                grid.tau_loc[I2V(iip,nt-1,kkr)] = std::max(_0_CR, _2_CR*grid.tau_loc[I2V(iip,nt-2,kkr)] - grid.tau_loc[I2V(iip,nt-3,kkr)]);
            else
                grid.tau_loc[I2V(iip,nt-1,kkr)] = std::min(_0_CR, _2_CR*grid.tau_loc[I2V(iip,nt-2,kkr)] - grid.tau_loc[I2V(iip,nt-3,kkr)]);
        } else {
            grid.tau_loc[I2V(iip,nt-1,kkr)] = _0_CR;
        }
    }

    // Bottom
    if (kkr == 0 && grid.k_first()) {
        if (!grid.is_changed[I2V(iip,jjt,0)]) {
            if (grid.tau_loc[I2V(iip,jjt,2)] >= 0)
                grid.tau_loc[I2V(iip,jjt,0)] = std::max(_0_CR, _2_CR*grid.tau_loc[I2V(iip,jjt,1)] - grid.tau_loc[I2V(iip,jjt,2)]);
            else
                grid.tau_loc[I2V(iip,jjt,0)] = std::min(_0_CR, _2_CR*grid.tau_loc[I2V(iip,jjt,1)] - grid.tau_loc[I2V(iip,jjt,2)]);
        } else {
            grid.tau_loc[I2V(iip,jjt,0)] = _0_CR;
        }
    }

    // Top
    if (kkr == nr-1 && grid.k_last())
        grid.tau_loc[I2V(iip,jjt,nr-1)] = _0_CR;

}


void Iterator::set_sweep_direction(int iswp) {
    // set sweep direction
    if (iswp == 0) {
        r_dirc = -1;
        t_dirc = -1;
        p_dirc = -1;
    } else if (iswp == 1) {
        r_dirc = -1;
        t_dirc = -1;
        p_dirc =  1;
    } else if (iswp == 2) {
        r_dirc = -1;
        t_dirc =  1;
        p_dirc = -1;
    } else if (iswp == 3) {
        r_dirc = -1;
        t_dirc =  1;
        p_dirc =  1;
    } else if (iswp == 4) {
        r_dirc =  1;
        t_dirc = -1;
        p_dirc = -1;
    } else if (iswp == 5) {
        r_dirc =  1;
        t_dirc = -1;
        p_dirc =  1;
    } else if (iswp == 6) {
        r_dirc =  1;
        t_dirc =  1;
        p_dirc = -1;
    } else if (iswp == 7) {
        r_dirc =  1;
        t_dirc =  1;
        p_dirc =  1;
    }
}

