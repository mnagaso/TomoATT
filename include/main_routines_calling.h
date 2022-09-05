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
        out_main.open("objective_function.txt");
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

    /////////////////////
    // loop for inversion
    /////////////////////

    bool line_search_mode = false; // if true, run_simulation_one_step skips adjoint simulation and only calculates objective function value

    // objective function for all src
    CUSTOMREAL v_obj = 0.0, old_v_obj = 0.0;

    for (int i_inv = 0; i_inv < IP.get_max_iter_inv(); i_inv++) {

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
        IP.write_src_rec_file();

        ///////////////
        // model update
        ///////////////

        if (IP.get_run_mode()==DO_INVERSION) {
            if (optim_method == GRADIENT_DESCENT)
                model_optimize(IP, grid, io, i_inv, v_obj, old_v_obj, first_src, out_main);
            else if (optim_method == LBFGS_MODE)
                model_optimize_lbfgs(IP, grid, io, i_inv, v_obj, first_src, out_main);
            else if (optim_method == HALVE_STEPPING_MODE)
                model_optimize_halve_stepping(IP, grid, io, i_inv, v_obj, first_src, out_main);
        }

        // output updated model
        if (subdom_main && id_sim==0) {
            io.change_xdmf_obj(0); // change xmf file for next src

            // write out model info
            io.write_fun(grid, i_inv);
            io.write_xi(grid, i_inv);
            io.write_eta(grid, i_inv);
            io.write_b(grid, i_inv);
            io.write_c(grid, i_inv);
            io.write_f(grid, i_inv);
        }

        // writeout temporary xdmf file
        io.update_xdmf_file(IP.src_ids_this_sim.size());

    } // end loop inverse

    // close xdmf file
    io.finalize_data_output_file(IP.src_ids_this_sim.size());

}


// run earthquake relocation mode
inline void run_earthquake_relocation(InputParams& IP, Grid& grid, IO_utils& io) {

//    if(myrank == 0)
//        std::cout << "size of src_list: " << IP.src_ids_this_sim.size() << std::endl;
//
//    // prepare output for iteration status
//    std::ofstream out_main;
//    if(myrank == 0 && id_sim ==0)
//        out_main.open("objective_function.txt");
//
//    // calculate traveltime field with input model. (stored in grid.T_loc)
//    bool first_src = true;
//    int i_inv_dummy = 0;
//    bool line_search_mode_dummy = false;
//    CUSTOMREAL v_obj = run_simulation_one_step(IP, grid, io, i_inv_dummy, first_src, line_search_mode_dummy);

    // calculate traveltime at each receiver (swapped from source)

    // calculate the approximated optimal origin time.

    // calculate gradient of travel time field at each source location.

    // calculate gradient of objective function at each source location.

    // run step-size-controlled gradient descent for searching optimal source location.

}


#endif // MAIN_ROUTINES_CALLING_H