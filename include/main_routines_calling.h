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
            if (IP.get_is_srcrec_swap())
                tmp.append(std::to_string(IP.N_cr_dif_local_data));
            else
                tmp.append(std::to_string(IP.N_cs_dif_local_data));
            tmp.append("),");
            out_main << std::setw(20) << tmp;

            tmp = "obj_cr_dif(";
            if (IP.get_is_srcrec_swap())
                tmp.append(std::to_string(IP.N_cs_dif_local_data));
            else
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
            out_main << std::setw(6)  << "it,"        << std::setw(6)  << "subit,"  << std::setw(16) << "step_length," << std::setw(16) << "q_0," << std::setw(16) << "q_t," \
                     << std::setw(16) << "v_obj_reg," << std::setw(16) << "qp_0,"  << std::setw(16) << "qp_t,"       << std::setw(16) << "td,"  << std::setw(16) << "tg," \
                     << std::setw(16) << "c1*qp_0,"    << std::setw(16) << "c2*qp_0," << std::setw(6)  << "step ok,"    << std::endl;

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

    if (id_sim==0) {
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

    /////////////////////
    // loop for inversion
    /////////////////////

    bool line_search_mode = false; // if true, run_simulation_one_step skips adjoint simulation and only calculates objective function value

    // objective function for all src
    CUSTOMREAL v_obj = 0.0, old_v_obj = 0.0;
    std::vector<CUSTOMREAL> v_obj_misfit(20, 0.0);

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
                model_optimize(IP, grid, io, i_inv, v_obj, old_v_obj, first_src, out_main);
            else if (optim_method == HALVE_STEPPING_MODE)
                v_obj_misfit = model_optimize_halve_stepping(IP, grid, io, i_inv, v_obj, first_src, out_main);
            else if (optim_method == LBFGS_MODE) {
                bool found_next_step = model_optimize_lbfgs(IP, grid, io, i_inv, v_obj, first_src, out_main);
                if (!found_next_step)
                    goto end_of_inversion;
            }
        }

        // output station correction file (only for teleseismic differential data)
        IP.write_station_correction_file(i_inv + 1);

        // output updated model
        if (id_sim==0) {
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

        write_objective_function(IP, i_inv, v_obj_misfit, out_main, "model update");

        // writeout temporary xdmf file
        io.update_xdmf_file();

        // wait for all processes to finish
        synchronize_all_world();

    } // end loop inverse

end_of_inversion:

    // close xdmf file
    io.finalize_data_output_file();

}


// run earthquake relocation mode
inline void run_earthquake_relocation(InputParams& IP, Grid& grid, IO_utils& io) {

    Receiver recs;

    // calculate traveltime for each receiver (swapped from source) and write in output file
    calculate_traveltime_for_all_src_rec(IP, grid, io);

    // prepare output for iteration status
    std::ofstream out_main; // close() is not mandatory
    prepare_header_line(IP, out_main);

    // objective function and its gradient
    CUSTOMREAL v_obj      = 999999999.0;

    int        i_iter     = 0;

    std::vector<CUSTOMREAL> v_obj_misfit;

    // iterate
    while (true) {

        v_obj      = 0.0;

        // calculate gradient of objective function at sources
        v_obj_misfit = calculate_gradient_objective_function(IP, grid, io, i_iter);
        v_obj = v_obj_misfit[0];

        // update source location
        recs.update_source_location(IP, grid);

        synchronize_all_world();

        // check convergence
        bool finished = false;

        if (subdom_main && id_subdomain==0) {
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
            // write objective function
            std::cout << "iteration: " << i_iter << ", objective function: "              << v_obj << std::endl;
        }

        // write objective functions
        write_objective_function(IP, i_iter, v_obj_misfit, out_main, "relocation");
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
    prepare_header_line(IP, out_main);

    if (id_sim==0) {
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

    /////////////////////
    // main loop for model update and relocation
    /////////////////////

    CUSTOMREAL v_obj = 0.0, old_v_obj = 10000000000.0;
    std::vector<CUSTOMREAL> v_obj_misfit(20, 0.0);

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
                model_optimize(IP, grid, io, i_inv, v_obj, old_v_obj, first_src, out_main);
            else if (optim_method == HALVE_STEPPING_MODE)
                v_obj_misfit = model_optimize_halve_stepping(IP, grid, io, i_inv, v_obj, first_src, out_main);
            else if (optim_method == LBFGS_MODE) {
                bool found_next_step = model_optimize_lbfgs(IP, grid, io, i_inv, v_obj, first_src, out_main);

                if (!found_next_step)
                    break;
            }

            // define old_v_obj
            old_v_obj = v_obj;

            // output objective function
            write_objective_function(IP, i_inv, v_obj_misfit, out_main, "model update");

            // since model is update. The written traveltime field should be discraded.
            // initialize is_T_written_into_file
            for (int i_src = 0; i_src < IP.n_src_this_sim_group; i_src++){
                const std::string name_sim_src = IP.get_src_name(i_src);

                if (proc_store_srcrec) // only proc_store_srcrec has the src_map object
                    IP.src_map[name_sim_src].is_T_written_into_file = false;
            }

            // output updated model
            if (id_sim==0) {
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

            } // end output updated model

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
        if (proc_store_srcrec){
            for(auto iter = IP.rec_map.begin(); iter != IP.rec_map.end(); iter++){
                iter->second.is_stop = false;
            }
        }

        // iterate
        for (int one_loop_i_iter = 0; one_loop_i_iter < IP.get_relocation_N_iter(); one_loop_i_iter++){
            int i_iter = i_loop * IP.get_relocation_N_iter() + one_loop_i_iter;

            if(myrank == 0 && id_sim ==0){
                std::cout   << std::endl;
                std::cout   << "loop " << i_loop+1 << ", relocation iteration " << one_loop_i_iter+1
                            << " ( the " << i_iter+1 << "-th relocation) starting ... " << std::endl;
                std::cout   << std::endl;
            }


            // calculate gradient of objective function at sources
            v_obj_misfit = calculate_gradient_objective_function(IP, grid, io, i_iter);
            v_obj = v_obj_misfit[0];

            // update source location
            recs.update_source_location(IP, grid);

            synchronize_all_world();


            // output location information
            if(id_sim == 0 && myrank == 0){
                std::cout << "objective function: "              << v_obj << std::endl;
            }

            // write objective functions
            write_objective_function(IP, i_iter, v_obj_misfit, out_main, "relocation");

            // write out new src_rec_file
            if (IP.get_if_output_in_process_data()){
                IP.write_src_rec_file(model_update_step,relocation_step);
            } else if (i_iter == IP.get_max_loop() * IP.get_relocation_N_iter()-1 || i_iter==0) {
                IP.write_src_rec_file(model_update_step,relocation_step);
            }

            // modify the receiver's location for output
            IP.modify_swapped_source_location();

            relocation_step += 1;
        } // end relocation loop

        grid.rejunenate_abcf();     // (a,b/r^2,c/(r^2*cos^2),f/(r^2*cos)) -> (a,b,c,f)

    } // end loop for model update and relocation

    // close xdmf file
    io.finalize_data_output_file();

}



#endif // MAIN_ROUTINES_CALLING_H