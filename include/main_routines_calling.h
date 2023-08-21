#ifndef MAIN_ROUTINES_CALLING_H
#define MAIN_ROUTINES_CALLING_H

#include <iostream>
#include <memory>
#include "mpi_funcs.h"
#include "config.h"
#include "utils.h"
#include "input_params.h"
#include "grid.h"
#include "io.h"
#include "main_routines_inversion_mode.h"
#include "model_optimization_routines.h"
#include "main_routines_earthquake_relocation.h"
#include "iterator_selector.h"
#include "iterator.h"
#include "iterator_legacy.h"
#include "iterator_level.h"
#include "source.h"
#include "receiver.h"
#include "kernel.h"
#include "model_update.h"
#include "lbfgs.h"


// prepare header line of objective_funciton.txt
inline void prepare_header_line(InputParams &IP, std::ofstream &out_main) {
    // prepare output for iteration status
    if(myrank == 0 && id_sim ==0){
        out_main.open(output_dir + "/objective_function.txt");
        if (optim_method == GRADIENT_DESCENT || optim_method == HALVE_STEPPING_MODE){
            bool have_abs    = ( IP.data_type.find("abs")    != IP.data_type.end() );
            bool have_cs_dif = ( IP.data_type.find("cs_dif") != IP.data_type.end() );
            bool have_cr_dif = ( IP.data_type.find("cr_dif") != IP.data_type.end() );
            bool have_tele   = ( IP.data_type.find("tele")   != IP.data_type.end() );


            out_main << std::setw(8) << std::right << "# iter,";
            if (optim_method == HALVE_STEPPING_MODE)
                out_main << std::setw(8) << std::right << "subiter,";
            std::string tmp = "obj(";
            tmp.append(std::to_string(IP.N_data));
            tmp.append("),");
            out_main << std::setw(20) << tmp;
            if (have_abs){
                tmp = "obj_abs(";
                tmp.append(std::to_string(IP.N_abs_local_data));
                tmp.append("),");
                out_main << std::setw(20) << tmp;
            }
            if (have_cs_dif){
                tmp = "obj_cs_dif(";
                if (IP.get_is_srcrec_swap())
                    tmp.append(std::to_string(IP.N_cr_dif_local_data));
                else
                    tmp.append(std::to_string(IP.N_cs_dif_local_data));
                tmp.append("),");
                out_main << std::setw(20) << tmp;
            }
            if (have_cr_dif){
                tmp = "obj_cr_dif(";
                if (IP.get_is_srcrec_swap())
                    tmp.append(std::to_string(IP.N_cs_dif_local_data));
                else
                    tmp.append(std::to_string(IP.N_cr_dif_local_data));
                tmp.append("),");
                out_main << std::setw(20) << tmp;
            }
            if (have_tele){
                tmp = "obj_tele(";
                tmp.append(std::to_string(IP.N_teleseismic_data));
                tmp.append("),");
                out_main << std::setw(20) << tmp;
            }
            out_main << std::setw(20) << "misfit,";
            if (have_abs){
                out_main << std::setw(20) << "misfit_abs,";
            }
            if (have_cs_dif){
                out_main << std::setw(20) << "misfit_cs_dif,";
            }
            if (have_cr_dif){
                out_main << std::setw(20) << "misfit_cr_dif,";
            }
            if (have_tele){
                out_main << std::setw(20) << "misfit_tele,";
            }
            out_main << std::setw(20) << "step_length," << std::endl;
        } else if (optim_method == LBFGS_MODE)
            out_main << std::setw(6)  << "it,"        << std::setw(6)  << "subit,"  << std::setw(16) << "step_length," << std::setw(16) << "qpt," << std::setw(16) << "v_obj_new," \
                     << std::setw(16) << "v_obj_reg," << std::setw(16) << "q_new,"  << std::setw(16) << "q_k,"       << std::setw(16) << "td,"  << std::setw(16) << "tg," \
                     << std::setw(16) << "c1*q_k,"    << std::setw(16) << "c2*q_k," << std::setw(6)  << "step ok,"    << std::endl;

    }
}


//
// run forward-only or inversion mode
//
inline void run_forward_only_or_inversion(InputParams &IP, Grid &grid, IO_utils &io) {

    // for check if the current source is the first source
    bool first_src = true;

    if(myrank == 0)
        std::cout << "id_sim: " << id_sim << ", size of src_map: " << IP.src_map.size() << std::endl;

    // prepare objective_function file
    std::ofstream out_main; // close() is not mandatory
    prepare_header_line(IP, out_main);


    if (subdom_main && id_sim==0) {
        //io.change_xdmf_obj(0); // change xmf file for next src
        io.change_group_name_for_model();

        // write out model info
        io.write_vel(grid, 0);
        io.write_xi( grid, 0);
        io.write_eta(grid, 0);
        //io.write_zeta(grid, i_inv); // TODO

        if (IP.get_verbose_output_level()){
            io.write_a(grid,   0);
            io.write_b(grid,   0);
            io.write_c(grid,   0);
            io.write_f(grid,   0);
            io.write_fun(grid, 0);
        }

        // output model_parameters_inv_0000.dat
        if (IP.get_if_output_model_dat())
            io.write_concerning_parameters(grid, 0, IP);
    }



    // output station correction file (only for teleseismic differential data)
    IP.write_station_correction_file(0);

    synchronize_all_world();

    /////////////////////
    // loop for inversion
    /////////////////////

    bool line_search_mode = false; // if true, run_simulation_one_step skips adjoint simulation and only calculates objective function value

    // objective function for all src
    CUSTOMREAL v_obj = 0.0, old_v_obj = 0.0;
    std::vector<CUSTOMREAL> v_obj_misfit;

    for (int i_inv = 0; i_inv < IP.get_max_iter_inv(); i_inv++) {

        if(myrank == 0 && id_sim ==0){
            std::cout << "iteration " << i_inv << " starting ... " << std::endl;
        }

        old_v_obj = v_obj;

        // prepare inverstion iteration group in xdmf file
        io.prepare_grid_inv_xdmf(i_inv);

        ///////////////////////////////////////////////////////
        // run (forward and adjoint) simulation for each source
        ///////////////////////////////////////////////////////

        // run forward and adjoint simulation and calculate current objective function value and sensitivity kernel for all sources
        line_search_mode = false;
        // skip for the mode with sub-iteration
        if (i_inv > 0 && optim_method != GRADIENT_DESCENT) {
        } else {
            bool is_read_time = false; 
            v_obj_misfit = run_simulation_one_step(IP, grid, io, i_inv, first_src, line_search_mode, is_read_time);
            v_obj = v_obj_misfit[0];
        }

        // wait for all processes to finish
        synchronize_all_world();

        // check if v_obj is nan
        if (std::isnan(v_obj)) {
            if (myrank == 0)
                std::cout << "v_obj is nan, stop inversion" << std::endl;

            // stop inversion
            break;
        }

        // output src rec file with the result arrival times
        if (IP.get_if_output_in_process_data()){
            IP.write_src_rec_file(i_inv,0);
        } else if (i_inv == IP.get_max_iter_inv()-1 || i_inv==0) {
            IP.write_src_rec_file(i_inv,0);
        }

        ///////////////
        // model update
        ///////////////
        if (IP.get_run_mode() == DO_INVERSION) {
            if (optim_method == GRADIENT_DESCENT)
                model_optimize(IP, grid, io, i_inv, v_obj, old_v_obj, v_obj_misfit, first_src, out_main);
            else if (optim_method == LBFGS_MODE)
                model_optimize_lbfgs(IP, grid, io, i_inv, v_obj, first_src, out_main);
            else if (optim_method == HALVE_STEPPING_MODE)
                model_optimize_halve_stepping(IP, grid, io, i_inv, v_obj, first_src, out_main);
        }

        // output station correction file (only for teleseismic differential data)
        IP.write_station_correction_file(i_inv + 1);

        // output updated model
        if (subdom_main && id_sim==0) {
            //io.change_xdmf_obj(0); // change xmf file for next src
            io.change_group_name_for_model();

            // write out model info
            if (IP.get_if_output_in_process() || i_inv >= IP.get_max_iter_inv() - 2){
                io.write_vel(grid, i_inv+1);
                io.write_xi( grid, i_inv+1);
                io.write_eta(grid, i_inv+1);
            }
            //io.write_zeta(grid, i_inv); // TODO

            if (IP.get_verbose_output_level()){
                io.write_a(grid,   i_inv+1);
                io.write_b(grid,   i_inv+1);
                io.write_c(grid,   i_inv+1);
                io.write_f(grid,   i_inv+1);
                io.write_fun(grid, i_inv+1);
            }

            // output model_parameters_inv_0000.dat
            if (IP.get_if_output_model_dat() \
            && (IP.get_if_output_in_process() || i_inv >= IP.get_max_iter_inv() - 2))
                io.write_concerning_parameters(grid, i_inv + 1, IP);

        }

        // writeout temporary xdmf file
        io.update_xdmf_file();

        // wait for all processes to finish
        synchronize_all_world();

    } // end loop inverse

    // close xdmf file
    io.finalize_data_output_file();

}


// run earthquake relocation mode
inline void run_earthquake_relocation(InputParams& IP, Grid& grid, IO_utils& io) {

    Receiver recs;
    int nrec_total = IP.rec_map_all.size();
    // broadcast
    broadcast_i_single_inter_and_intra_sim(nrec_total, 0);

    // calculate traveltime for each receiver (swapped from source) and write in output file
    calculate_traveltime_for_all_src_rec(IP, grid, io);

    //std::cout << "cpk-sp, id_sim: " << id_sim << ", myrank: " << myrank << ", Nrec: " << IP.rec_map.size() << std::endl;
    //for(auto iter = IP.rec_map.begin(); iter != IP.rec_map.end(); iter++){
    //    std::cout << "cpk-sp, id_sim: " << id_sim << ", myrank: " << myrank << ", name: " << iter->first << std::endl;
    //}

    // prepare output for iteration status
    std::ofstream out_main;
    if(myrank == 0 && id_sim ==0){
        out_main.open(output_dir + "/objective_function_reloc.txt");

        bool have_abs    = ( IP.data_type.find("abs")    != IP.data_type.end() );
        bool have_cs_dif = ( IP.data_type.find("cs_dif") != IP.data_type.end() );
        bool have_cr_dif = ( IP.data_type.find("cr_dif") != IP.data_type.end() );
        bool have_tele   = ( IP.data_type.find("tele")   != IP.data_type.end() );

        out_main << std::setw(8) << std::right << "# iter,";
        out_main << std::setw(16) << std::right << "N_reloc,";
        out_main << std::setw(16) << std::right << "N_located,";
        
        std::string tmp = "obj(";
        tmp.append(std::to_string(IP.N_data));
        tmp.append("),");
        out_main << std::setw(20) << tmp;
        if (have_abs){
            tmp = "obj_abs(";
            tmp.append(std::to_string(IP.N_abs_local_data));
            tmp.append("),");
            out_main << std::setw(20) << tmp;
        }
        if (have_cs_dif){
            tmp = "obj_cs_dif(";
            tmp.append(std::to_string(IP.N_cs_dif_local_data));
            tmp.append("),");
            out_main << std::setw(20) << tmp;
        }
        if (have_cr_dif){
            tmp = "obj_cr_dif(";
            tmp.append(std::to_string(IP.N_cr_dif_local_data));
            tmp.append("),");
            out_main << std::setw(20) << tmp;
        }
        if (have_tele){
            tmp = "obj_tele(";
            tmp.append(std::to_string(IP.N_teleseismic_data));
            tmp.append("),");
            out_main << std::setw(20) << tmp;
        }
        out_main << std::endl;
    }

    // objective function and its gradient
    CUSTOMREAL v_obj      = 999999999.0;
    CUSTOMREAL v_obj_old  = 0.0;
    CUSTOMREAL v_obj_grad = 0.0;
    CUSTOMREAL v_obj_abs      = 999999999.0;
    CUSTOMREAL v_obj_cr       = 999999999.0;
    CUSTOMREAL v_obj_cs       = 999999999.0;

    int        i_iter     = 0;

    // iterate
    while (true) {

        v_obj_old  = v_obj;
        v_obj      = 0.0;
        v_obj_grad = 0.0;
        v_obj_abs  = 0.0;
        v_obj_cr   = 0.0;
        v_obj_cs   = 0.0;
        // determine which earthquake should be located
        IP.name_for_reloc.clear();
        for(auto iter = IP.rec_map.begin(); iter != IP.rec_map.end(); iter++){
            if (!iter->second.is_stop)
                IP.name_for_reloc.push_back(iter->first);
        }

        // calculate gradient of objective function at sources
        calculate_gradient_objective_function(IP, grid, io, i_iter);

        // update source location
        recs.update_source_location(IP, grid);

        synchronize_all_world();

        // calculate sum of objective function and gradient
        for (auto iter = IP.rec_map.begin(); iter != IP.rec_map.end(); iter++) {
            v_obj      += iter->second.vobj_src_reloc;
            v_obj_grad += iter->second.vobj_grad_norm_src_reloc;
            v_obj_abs  += iter->second.vobj_src_reloc_data[0];
            v_obj_cr   += iter->second.vobj_src_reloc_data[1];
            v_obj_cs   += iter->second.vobj_src_reloc_data[2];
        }

        // reduce
        //allreduce_cr_sim_single_inplace(v_obj);
        //allreduce_cr_sim_single_inplace(v_obj_grad);

        // check convergence
        int count_loc = 0;
        bool finished = false;

        if (subdom_main && id_subdomain==0) {
            if (IP.name_for_reloc.size() == 0){
                std::cout << "Finished relocation because all receivers have been located." << std::endl;
                finished = true;
            }

            if (i_iter >= N_ITER_MAX_SRC_RELOC){
                std::cout << "Finished relocation because iteration number exceeds the maximum " << N_ITER_MAX_SRC_RELOC << std::endl;
                finished = true;
            }
            allreduce_bool_single_inplace_sim(finished); //LAND
        }

        synchronize_all_world();

        // check if all processes have finished
        broadcast_bool_inter_and_intra_sim(finished, 0);

        synchronize_all_world();

        // output location information
        if(id_sim == 0 && myrank == 0){
            // number of receiver which have been completed
            int n_relocated = IP.rec_map.size() - IP.name_for_reloc.size();

            // write objective function
            std::cout << "iteration: " << i_iter << ", objective function: "              << v_obj
                                                 << ", objective function old: "          << v_obj_old
                                                 << ", norm grad of relocating: "         << v_obj_grad // BUG: strangely high when simultaneous run
                                                 << ", average norm grad of relocating: " << v_obj_grad/IP.name_for_reloc.size()
                                                 << ", v_obj/n_src: "                     << v_obj/nrec_total
                                                 << ", diff_v/v_obj_old "                 << std::abs(v_obj-v_obj_old)/v_obj_old << std::endl;
            std::cout << "Earthquakes require location: " << n_relocated << " / " << nrec_total << " completed." << std::endl;

            // the last 10 sources under location
            if (nrec_total - count_loc < 10){
                std::cout << "Last 10 sources under location. names: ";
                for (int i = 0; i < (int)IP.name_for_reloc.size(); i++){
                    std::cout << IP.name_for_reloc[i] << ", ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        // write objective functions
        if (myrank==0 && id_sim==0) {
            out_main << std::setw(8) << std::right << i_iter + 1;
            out_main << "," << std::setw(15) << (int)IP.rec_map.size();
            out_main << "," << std::setw(15) << std::right << count_loc;
            out_main << "," << std::setw(19) << std::right << v_obj;
            out_main << "," << std::setw(19) << std::right << v_obj_abs;
            out_main << "," << std::setw(19) << std::right << 0.0;
            out_main << "," << std::setw(19) << std::right << v_obj_cr;
            out_main << "," << std::endl;
        }

        if (finished)
            break;

        // new iteration
        i_iter++;


    }


    // modify the receiver's location
    IP.modify_swapped_source_location();
    // write out new src_rec_file
    IP.write_src_rec_file(0,0);
    // close xdmf file
    io.finalize_data_output_file();


}

// prepare header line of objective_funciton.txt
inline void prepare_header_line_mode_3(InputParams &IP, std::ofstream &out_main) {
    // prepare output for iteration status
    if(myrank == 0 && id_sim ==0){
        out_main.open(output_dir + "/objective_function.txt");
        if (optim_method == GRADIENT_DESCENT || optim_method == HALVE_STEPPING_MODE){
            // bool have_abs    = ( IP.data_type.find("abs")    != IP.data_type.end() );
            // bool have_cs_dif = ( IP.data_type.find("cs_dif") != IP.data_type.end() );
            // bool have_cr_dif = ( IP.data_type.find("cr_dif") != IP.data_type.end() );
            // bool have_tele   = ( IP.data_type.find("tele")   != IP.data_type.end() );


            out_main << std::setw(8) << std::right << "# iter,";
            out_main << std::setw(13) << std::right << " type,";
            
            // if (optim_method == HALVE_STEPPING_MODE) 
            //     out_main << std::setw(8) << std::right << "subiter,";        (TODO in the future)
            std::string tmp = "obj(";
            tmp.append(std::to_string(IP.N_data));
            tmp.append("),");
            out_main << std::setw(20) << tmp;

            tmp = "obj_abs(";
            tmp.append(std::to_string(IP.N_abs_local_data));
            tmp.append("),");
            out_main << std::setw(20) << tmp;
            
            tmp = "obj_cs_dif(";
            tmp.append(std::to_string(IP.N_cs_dif_local_data));
            tmp.append("),");
            out_main << std::setw(20) << tmp;
        
            tmp = "obj_cr_dif(";
            tmp.append(std::to_string(IP.N_cr_dif_local_data));
            tmp.append("),");
            out_main << std::setw(20) << tmp;

            tmp = "obj_tele(";
            tmp.append(std::to_string(IP.N_teleseismic_data));
            tmp.append("),");
            out_main << std::setw(20) << tmp;

            out_main << std::setw(25) << "res(mean/std),";

            out_main << std::setw(25) << "res_abs(mean/std),";
        
            out_main << std::setw(25) << "res_cs_dif(mean/std),";

            out_main << std::setw(25) << "res_cr_dif(mean/std),";

            out_main << std::setw(25) << "res_tele(mean/std),";

            out_main << std::setw(20) << "step_length," << std::endl;

        } else if (optim_method == LBFGS_MODE)
            out_main << std::setw(6)  << "it,"        << std::setw(6)  << "subit,"  << std::setw(16) << "step_length," << std::setw(16) << "qpt," << std::setw(16) << "v_obj_new," \
                     << std::setw(16) << "v_obj_reg," << std::setw(16) << "q_new,"  << std::setw(16) << "q_k,"       << std::setw(16) << "td,"  << std::setw(16) << "tg," \
                     << std::setw(16) << "c1*q_k,"    << std::setw(16) << "c2*q_k," << std::setw(6)  << "step ok,"    << std::endl;

    }
}


// run earthquake relocation mode
inline void run_inversion_and_relocation(InputParams& IP, Grid& grid, IO_utils& io) {

    /////////////////////
    // preparation of model update
    /////////////////////

    // for check if the current source is the first source
    bool first_src = true;

    if(myrank == 0)
        std::cout << "id_sim: " << id_sim << ", size of src_map: " << IP.src_map.size() << std::endl;

    // prepare objective_function file
    std::ofstream out_main; // close() is not mandatory
    prepare_header_line_mode_3(IP, out_main);

    if (subdom_main && id_sim==0) {
        // prepare inverstion iteration group in xdmf file
        io.prepare_grid_inv_xdmf(0);

        //io.change_xdmf_obj(0); // change xmf file for next src
        io.change_group_name_for_model();

        // write out model info
        io.write_vel(grid, 0);
        io.write_xi( grid, 0);
        io.write_eta(grid, 0);
        //io.write_zeta(grid, i_inv); // TODO

        if (IP.get_verbose_output_level()){
            io.write_a(grid,   0);
            io.write_b(grid,   0);
            io.write_c(grid,   0);
            io.write_f(grid,   0);
            io.write_fun(grid, 0);
        }

        // output model_parameters_inv_0000.dat
        if (IP.get_if_output_model_dat())
            io.write_concerning_parameters(grid, 0, IP);
    }

    // output station correction file (only for teleseismic differential data)
    IP.write_station_correction_file(0);

    synchronize_all_world();

    bool line_search_mode = false; // if true, run_simulation_one_step skips adjoint simulation and only calculates objective function value

    /////////////////////
    // preparation of relocation
    /////////////////////

    Receiver recs;
    int nrec_total = IP.rec_map_all.size();
    // broadcast
    broadcast_i_single_inter_and_intra_sim(nrec_total, 0);


    /////////////////////
    // main loop for model update and relocation
    /////////////////////

    CUSTOMREAL v_obj = 0.0, old_v_obj = 0.0;
    std::vector<CUSTOMREAL> v_obj_misfit;

    int model_update_step = 0;
    int relocation_step = 0;

    for (int i_loop = 0; i_loop < IP.get_max_loop(); i_loop++){
        /////////////////////
        // stage 1, update model parameters 
        /////////////////////

        for (int one_loop_i_inv = 0; one_loop_i_inv < IP.get_model_update_N_iter(); one_loop_i_inv++){
            // the actural step index of model update
            int i_inv = i_loop * IP.get_model_update_N_iter() + one_loop_i_inv;

            if(myrank == 0 && id_sim ==0){
                std::cout   << std::endl; 
                std::cout   << "loop " << i_loop+1 << ", model update iteration " << one_loop_i_inv+1 
                            << " ( the " << i_inv+1 << "-th model update) starting ... " << std::endl;  
                std::cout   << std::endl; 
            }

            old_v_obj = v_obj;

            // prepare inverstion iteration group in xdmf file
            io.prepare_grid_inv_xdmf(i_inv);

            ///////////////////////////////////////////////////////
            // run (forward and adjoint) simulation for each source
            ///////////////////////////////////////////////////////

            // run forward and adjoint simulation and calculate current objective function value and sensitivity kernel for all sources
            line_search_mode = false;
            // skip for the mode with sub-iteration
            if (i_inv > 0 && optim_method != GRADIENT_DESCENT) {
            } else {
                bool is_read_time;
                if (i_loop > 0 && one_loop_i_inv == 0)
                    is_read_time = true; 
                else
                    is_read_time = false; 
                v_obj_misfit = run_simulation_one_step(IP, grid, io, i_inv, first_src, line_search_mode, is_read_time);
                v_obj = v_obj_misfit[0];
            }

            // wait for all processes to finish
            synchronize_all_world();

            // check if v_obj is nan
            if (std::isnan(v_obj)) {
                if (myrank == 0)
                    std::cout << "v_obj is nan, stop inversion" << std::endl;
                // stop inversion
                break;
            }

            ///////////////
            // model update
            ///////////////

            if (optim_method == GRADIENT_DESCENT)
                model_optimize(IP, grid, io, i_inv, v_obj, old_v_obj, v_obj_misfit, first_src, out_main);
            else if (optim_method == LBFGS_MODE)
                model_optimize_lbfgs(IP, grid, io, i_inv, v_obj, first_src, out_main);
            else if (optim_method == HALVE_STEPPING_MODE)
                model_optimize_halve_stepping(IP, grid, io, i_inv, v_obj, first_src, out_main);

            // output updated model
            if (subdom_main && id_sim==0) {
                //io.change_xdmf_obj(0); // change xmf file for next src
                io.change_group_name_for_model();

                // write out model info
                if (IP.get_if_output_in_process() || i_inv >= IP.get_max_loop()*IP.get_model_update_N_iter() - 2){
                    io.write_vel(grid, i_inv+1);
                    io.write_xi( grid, i_inv+1);
                    io.write_eta(grid, i_inv+1);
                }
                //io.write_zeta(grid, i_inv); // TODO

                if (IP.get_verbose_output_level()){
                    io.write_a(grid,   i_inv+1);
                    io.write_b(grid,   i_inv+1);
                    io.write_c(grid,   i_inv+1);
                    io.write_f(grid,   i_inv+1);
                    io.write_fun(grid, i_inv+1);
                }

                // output model_parameters_inv_0000.dat
                if (IP.get_if_output_model_dat() \
                && (IP.get_if_output_in_process() || i_inv >= IP.get_max_loop()*IP.get_model_update_N_iter() - 2))
                    io.write_concerning_parameters(grid, i_inv + 1, IP);         
            }  

            // writeout temporary xdmf file
            io.update_xdmf_file();

            // write src rec files
            if (IP.get_if_output_in_process_data()){
                IP.write_src_rec_file(model_update_step,relocation_step);
            } else if (i_inv == IP.get_max_loop() * IP.get_model_update_N_iter()-1 || i_inv==0) {
                IP.write_src_rec_file(model_update_step,relocation_step);
            }           
            model_update_step += 1;

            // wait for all processes to finish
            synchronize_all_world();
        }   // end model update in one loop




        /////////////////////
        // stage 2, update earthquake locations 
        /////////////////////

        if(myrank == 0 && id_sim ==0){
            std::cout   << std::endl; 
            std::cout   << "loop " << i_loop+1 << ", computing traveltime field for relocation" << std::endl;  
            std::cout   << std::endl; 
        }

        // calculate traveltime for each receiver (swapped from source) and write in output file
        calculate_traveltime_for_all_src_rec(IP, grid, io);

        

        // initilize all earthquakes
        for(auto iter = IP.rec_map.begin(); iter != IP.rec_map.end(); iter++){
            iter->second.is_stop = false;
        }

        // objective function and its gradient
        CUSTOMREAL v_obj_grad = 0.0;
        CUSTOMREAL v_obj_abs      = 999999999.0;
        CUSTOMREAL v_obj_cr       = 999999999.0;
        CUSTOMREAL v_obj_cs       = 999999999.0;
        CUSTOMREAL res            = 999999999.0;
        CUSTOMREAL res_sq         = 999999999.0;  
        CUSTOMREAL res_abs        = 999999999.0;
        CUSTOMREAL res_abs_sq     = 999999999.0; 
        CUSTOMREAL res_cs         = 999999999.0;
        CUSTOMREAL res_cs_sq      = 999999999.0; 
        CUSTOMREAL res_cr         = 999999999.0;
        CUSTOMREAL res_cr_sq      = 999999999.0; 

        // iterate
        for (int one_loop_i_iter = 0; one_loop_i_iter < IP.get_relocation_N_iter(); one_loop_i_iter++){
            int i_iter = i_loop * IP.get_model_update_N_iter() + one_loop_i_iter;

            if(myrank == 0 && id_sim ==0){
                std::cout   << std::endl; 
                std::cout   << "loop " << i_loop+1 << ", relocation iteration " << one_loop_i_iter+1 
                            << " ( the " << i_iter+1 << "-th relocation) starting ... " << std::endl;  
                std::cout   << std::endl; 
            }

            old_v_obj  = v_obj;
            v_obj      = 0.0;
            v_obj_grad = 0.0;
            v_obj_abs  = 0.0;
            v_obj_cs   = 0.0;
            v_obj_cr   = 0.0;
            res        = 0.0;
            res_sq     = 0.0;            
            res_abs    = 0.0;
            res_abs_sq = 0.0; 
            res_cs     = 0.0;
            res_cs_sq  = 0.0; 
            res_cr     = 0.0;
            res_cr_sq  = 0.0; 
            // determine which earthquake should be located
            IP.name_for_reloc.clear();
            for(auto iter = IP.rec_map.begin(); iter != IP.rec_map.end(); iter++){
                if (!iter->second.is_stop)
                    IP.name_for_reloc.push_back(iter->first);
            }

            // calculate gradient of objective function at sources
            calculate_gradient_objective_function(IP, grid, io, i_iter);

            // update source location
            recs.update_source_location(IP, grid);

            synchronize_all_world();

            // calculate sum of objective function and gradient
            for (auto iter = IP.rec_map.begin(); iter != IP.rec_map.end(); iter++) {
                v_obj       += iter->second.vobj_src_reloc;
                v_obj_grad  += iter->second.vobj_grad_norm_src_reloc;
                v_obj_abs   += iter->second.vobj_src_reloc_data[0];
                v_obj_cs    += iter->second.vobj_src_reloc_data[1];
                v_obj_cr    += iter->second.vobj_src_reloc_data[2];
                res         += iter->second.vobj_src_reloc_data[3] + iter->second.vobj_src_reloc_data[5] + iter->second.vobj_src_reloc_data[7];
                res_sq      += iter->second.vobj_src_reloc_data[4] + iter->second.vobj_src_reloc_data[6] + iter->second.vobj_src_reloc_data[8]; 
                res_abs     += iter->second.vobj_src_reloc_data[3];
                res_abs_sq  += iter->second.vobj_src_reloc_data[4];
                res_cs      += iter->second.vobj_src_reloc_data[5];
                res_cs_sq   += iter->second.vobj_src_reloc_data[6];
                res_cr      += iter->second.vobj_src_reloc_data[7];
                res_cr_sq   += iter->second.vobj_src_reloc_data[8];
            }

            synchronize_all_world();

            // output location information
            if(id_sim == 0 && myrank == 0){
            //     // number of receiver which have been completed
            //     int n_relocated = IP.rec_map.size() - IP.name_for_reloc.size();
                std::cout << "objective function: "              << v_obj << std::endl;
                // write objective function
                // std::cout << "iteration: " << i_iter << ", objective function: "              << v_obj
            //                                         << ", objective function old: "          << old_v_obj
            //                                         << ", norm grad of relocating: "         << v_obj_grad // BUG: strangely high when simultaneous run
            //                                         << ", average norm grad of relocating: " << v_obj_grad/IP.name_for_reloc.size()
            //                                         << ", v_obj/n_src: "                     << v_obj/nrec_total
            //                                         << ", diff_v/old_v_obj "                 << std::abs(v_obj-old_v_obj)/old_v_obj << std::endl;
            //     std::cout << "Earthquakes require location: " << n_relocated << " / " << nrec_total << " completed." << std::endl;
            }

            // write objective functions
            if (myrank==0 && id_sim==0) {
                out_main << std::setw(5) << std::right << i_iter + 1 << ",";
                out_main << std::setw(13) << "relocation";
                out_main << "," << std::setw(19) << std::right << v_obj;
                out_main << "," << std::setw(19) << std::right << v_obj_abs;
                out_main << "," << std::setw(19) << std::right << v_obj_cs;
                out_main << "," << std::setw(19) << std::right << v_obj_cr;
                out_main << "," << std::setw(19) << std::right << 0.0;
                CUSTOMREAL mean;
                CUSTOMREAL std;
                std::string tmp;
                // res
                
                if (IP.N_data > 0){
                    mean = res/IP.N_data;
                    std  = sqrt(res_sq/IP.N_data - my_square(mean));
                    tmp = std::to_string(mean);
                    tmp.append("/");
                    tmp.append(std::to_string(std));
                    out_main << "," << std::setw(24) << tmp;
                } else {
                    out_main << "," << std::setw(24) << "0.0/0.0";
                }
                // res_abs
                if (IP.N_abs_local_data > 0){
                    mean = res_abs/IP.N_abs_local_data;
                    std  = sqrt(res_abs_sq/IP.N_abs_local_data - my_square(mean));
                    tmp = std::to_string(mean);
                    tmp.append("/");
                    tmp.append(std::to_string(std));
                    out_main << "," << std::setw(24) << tmp;
                } else {
                    out_main << "," << std::setw(24) << "0.0/0.0";
                }
                // res_cs (swapped cr)
                if (IP.N_cr_dif_local_data > 0){
                    mean = res_cr/IP.N_cr_dif_local_data;
                    std  = sqrt(res_cr_sq/IP.N_cr_dif_local_data - my_square(mean));
                    tmp = std::to_string(mean);
                    tmp.append("/");
                    tmp.append(std::to_string(std));
                    out_main << "," << std::setw(24) << tmp;
                    
                } else {
                    out_main << "," << std::setw(24) << "0.0/0.0";
                }
                // res_cr (swapped cs)
                if (IP.N_cs_dif_local_data > 0){
                    mean = res_cs/IP.N_cs_dif_local_data;
                    std  = sqrt(res_cs_sq/IP.N_cs_dif_local_data - my_square(mean));
                    tmp = std::to_string(mean);
                    tmp.append("/");
                    tmp.append(std::to_string(std));
                    out_main << "," << std::setw(24) << tmp;
                } else {
                    out_main << "," << std::setw(24) << "0.0/0.0";
                } 
                // res_tele
                out_main << "," << std::setw(24) << "0.0/0.0";
                out_main << "," << std::endl;
            }

            // write out new src_rec_file
            if (IP.get_if_output_in_process_data()){
                IP.write_src_rec_file(model_update_step,relocation_step);
            } else if (i_iter == IP.get_max_loop() * IP.get_relocation_N_iter()-1 || i_iter==0) {
                IP.write_src_rec_file(model_update_step,relocation_step);
            }   
            // modify the receiver's location for output
            IP.modify_swapped_source_location();

            relocation_step += 1;
        }
 
        grid.rejunenate_abcf();     // (a,b/r^2,c/(r^2*cos^2),f/(r^2*cos)) -> (a,b,c,f)
    }

    // close xdmf file
    io.finalize_data_output_file();

}



#endif // MAIN_ROUTINES_CALLING_H