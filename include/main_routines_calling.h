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


// run forward-only or inversion mode
inline void run_forward_only_or_inversion(InputParams &IP, Grid &grid, IO_utils &io) {

    // for check if the current source is the first source
    bool first_src = true;

    if(myrank == 0)
        std::cout << "size of src_list: " << IP.src_ids_this_sim.size() << std::endl;

    // prepare output for iteration status
    std::ofstream out_main;
    if(myrank == 0 && id_sim ==0){
        out_main.open(output_dir + "objective_function.txt");
        if (optim_method == GRADIENT_DESCENT)
            out_main << std::setw(6) << "iter," << std::setw(16) << "v_obj," << std::setw(16) << "step_size," << std::endl;
        else if (optim_method == LBFGS_MODE)
            out_main << std::setw(6)  << "it,"        << std::setw(6)  << "subit,"  << std::setw(16) << "step_size," << std::setw(16) << "qpt," << std::setw(16) << "v_obj_new," \
                     << std::setw(16) << "v_obj_reg," << std::setw(16) << "q_new,"  << std::setw(16) << "q_k,"       << std::setw(16) << "td,"  << std::setw(16) << "tg," \
                     << std::setw(16) << "c1*q_k,"    << std::setw(16) << "c2*q_k," << std::setw(6)  << "step ok"    << std::endl;
        else if (optim_method == HALVE_STEPPING_MODE)
            out_main << std::setw(6)  << "it,"       << std::setw(6)  << "subit,"     << std::setw(16) << "step_size," \
                     << std::setw(16) << "diff_obj," << std::setw(16) << "v_obj_new," << std::setw(16) << "v_obj_old," << std::endl;

    }

    if (subdom_main && id_sim==0 && IP.get_is_output_model_dat()==1) {
        io.write_concerning_parameters(grid, 0);
    }

    synchronize_all_world();

    /////////////////////
    // loop for inversion
    /////////////////////

    bool line_search_mode = false; // if true, run_simulation_one_step skips adjoint simulation and only calculates objective function value

    // objective function for all src
    CUSTOMREAL v_obj = 0.0, old_v_obj = 0.0;

    for (int i_inv = 0; i_inv < IP.get_max_iter_inv(); i_inv++) {

        if(myrank == 0 && id_sim ==0){
            std::cout << "iteration " << i_inv << " starting ... " << std::endl;
        }

        old_v_obj = v_obj;

        ///////////////////////////////////////////////////////
        // run (forward and adjoint) simulation for each source
        ///////////////////////////////////////////////////////

        // run forward and adjoint simulation and calculate current objective function value and sensitivity kernel for all sources
        line_search_mode = false;
        // skip for the  mode with sub-iteration
        if (i_inv > 0 && optim_method != GRADIENT_DESCENT) {
        } else {
            v_obj = run_simulation_one_step(IP, grid, io, i_inv, first_src, line_search_mode);
        }

        // wait for all processes to finish
        synchronize_all_world();

        // output src rec file with the result arrival times
        IP.write_src_rec_file(i_inv);

        ///////////////
        // model update
        ///////////////

        if (IP.get_run_mode() == DO_INVERSION) {
            if (optim_method == GRADIENT_DESCENT)
                model_optimize(IP, grid, io, i_inv, v_obj, old_v_obj, first_src, out_main);
            else if (optim_method == LBFGS_MODE)
                model_optimize_lbfgs(IP, grid, io, i_inv, v_obj, first_src, out_main);
            else if (optim_method == HALVE_STEPPING_MODE)
                model_optimize_halve_stepping(IP, grid, io, i_inv, v_obj, first_src, out_main);
        }

        // output updated model
        if (subdom_main && id_sim==0) {
            if (IP.get_is_output_source_field()){
                io.change_xdmf_obj(0); // change xmf file for next src

                // write out model info
                io.write_fun(grid, i_inv);
                io.write_xi(grid, i_inv);
                io.write_eta(grid, i_inv);
                io.write_b(grid, i_inv);
                io.write_c(grid, i_inv);
                io.write_f(grid, i_inv);
            }

            if (IP.get_is_output_model_dat())        // output model_parameters_inv_0000.dat
                io.write_concerning_parameters(grid, i_inv + 1);
        }

        // writeout temporary xdmf file
        if (IP.get_is_output_source_field())
            io.update_xdmf_file(IP.src_ids_this_sim.size());

        // wait for all processes to finish
        synchronize_all_world();

    } // end loop inverse

    // close xdmf file
    if (IP.get_is_output_source_field())
        io.finalize_data_output_file(IP.src_ids_this_sim.size());

}


// run earthquake relocation mode
inline void run_earthquake_relocation(InputParams& IP, Grid& grid, IO_utils& io) {

    // this routine is not supporting simultaneous run
    if (n_sims > 1) {
        std::cout << "Earthquake relocation mode is not supporting simultaneous run" << std::endl;
        exit(1);
    }

    Receiver recs;

    // calculate traveltime for each receiver (swapped from source) and write in output file
    calculate_traveltime_for_all_src_rec(IP, grid, io);

    // create a unique receiver list among all sources
    // while creating this list, each receiver object stores the id of correspoinding receiver in this unique list
    std::vector<SrcRec> unique_rec_list = create_unique_rec_list(IP);

    // objective function and its gradient
    CUSTOMREAL v_obj = 0.0, v_obj_old = 0.0;
    CUSTOMREAL v_obj_grad = 0.0;
    int i_iter = 0;

    // iterate
    while (true) {

        v_obj_old = v_obj;
        v_obj = 0.0;
        v_obj_grad = 0.0;

        // calculate gradient of objective function at sources
        calculate_gradient_objective_function(IP, grid, io, unique_rec_list);

        // update source location
        for (long unsigned int i_src = 0; i_src < IP.src_ids_this_sim.size(); i_src++){
            id_sim_src = IP.src_ids_this_sim[i_src];
            recs.update_source_location(IP, grid, unique_rec_list);
        }

        // calculate sum of objective function and gradient
        for (auto& rec : unique_rec_list) {
            v_obj      += rec.vobj_src_reloc;
            v_obj_grad += rec.vobj_grad_norm_src_reloc;
        }

        // write objective functions
        if(myrank == 0){
            // write objective function
            std::cout << "iteration: " << i_iter << " objective function: " << v_obj \
                                                 << " v_obj_grad: " << v_obj_grad \
                                                 << " v_obj/n_src: " << v_obj/unique_rec_list.size() \
                                                 << " diff_v/v_obj_old " << std::abs(v_obj-v_obj_old)/v_obj_old << std::endl;
        }


        // check convergence
        if (i_iter > N_ITER_MAX_SRC_RELOC || v_obj/unique_rec_list.size() < TOL_SRC_RELOC)
            break;

        i_iter++;
    }

    // modify the receiver's location


    // write out new src_rec_file
    IP.write_src_rec_file(0);

    // close xdmf file
    if (IP.get_is_output_source_field())
        io.finalize_data_output_file(IP.src_ids_this_sim.size());
}


#endif // MAIN_ROUTINES_CALLING_H