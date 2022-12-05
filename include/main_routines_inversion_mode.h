#ifndef MAIN_ROUTINES_INVERSION_MODE_H
#define MAIN_ROUTINES_INVERSION_MODE_H

#include <iostream>
#include <memory>
#include "mpi_funcs.h"
#include "config.h"
#include "utils.h"
#include "input_params.h"
#include "grid.h"
#include "io.h"
#include "iterator_selector.h"
#include "iterator.h"
#include "iterator_legacy.h"
#include "iterator_level.h"
#include "source.h"
#include "receiver.h"
#include "kernel.h"
#include "model_update.h"
#include "lbfgs.h"


// run forward and adjoint simulation and calculate current objective function value and sensitivity kernel if requested
inline CUSTOMREAL run_simulation_one_step(InputParams& IP, Grid& grid, IO_utils& io, int i_inv, bool& first_src, bool line_search_mode){

    CUSTOMREAL v_obj = _0_CR;

    // initialize kernel arrays
    if (IP.get_run_mode() == DO_INVERSION)
        grid.initialize_kernels();

    // reinitialize factors
    grid.reinitialize_abcf();

    ///////////////////////
    // loop for each source
    ///////////////////////

    for (long unsigned int i_src = 0; i_src < IP.src_ids_this_sim.size(); i_src++) {

        // load the global id of this src
        id_sim_src = IP.src_ids_this_sim[i_src]; // local src id to global src id

        if (myrank == 0)
            std::cout << "source id: " << id_sim_src << ", forward modeling starting..." << std::endl;

        // set group name to be used for output in h5
        io.change_group_name_for_source();

        // output initial field
        if(first_src && IP.get_is_output_source_field()) {
            if (subdom_main) {
                // write true solution
                if (if_test){
                    io.write_true_solution(grid);
                }
                // write initial velocity model
                //io.write_velocity_model_h5(grid);
            }
            first_src = false;
        }

        // (re) initialize source object and set to grid
        Source src(IP, grid);

        // initialize iterator object
        bool first_init = (i_inv == 0 && i_src==0);

        /////////////////////////
        // run forward simulation
        /////////////////////////

        // initialize iterator object
        std::unique_ptr<Iterator> It;

        if (!hybrid_stencil_order){
            select_iterator(IP, grid, src, io, first_init, It, false);
            It->run_iteration_forward(IP, grid, io, first_init);
        } else {
            // hybrid stencil mode
            std::cout << "\nrunnning in hybrid stencil mode\n" << std::endl;

            // run 1st order forward simulation
            std::unique_ptr<Iterator> It_pre;
            IP.set_stencil_order(1);
            IP.set_conv_tol(IP.get_conv_tol()*100.0);
            select_iterator(IP, grid, src, io, first_init, It_pre, false);
            It_pre->run_iteration_forward(IP, grid, io, first_init);

            // run 3rd order forward simulation
            IP.set_stencil_order(3);
            IP.set_conv_tol(IP.get_conv_tol()/100.0);
            select_iterator(IP, grid, src, io, first_init, It, true);
            It->run_iteration_forward(IP, grid, io, first_init);
        }

        // output the result of forward simulation
        // ignored for inversion mode.
        if (subdom_main && !line_search_mode && IP.get_is_output_source_field()) {

            // output T (result timetable)
            io.write_T(grid, i_inv);

        }

        // calculate the arrival times at each receivers
        Receiver recs;
        recs.calculate_arrival_time(IP, grid);

        /////////////////////////
        // run adjoint simulation
        /////////////////////////

        if (IP.get_run_mode()==DO_INVERSION){
            // calculate adjoint source
            v_obj += recs.calculate_adjoint_source(IP);

            // run iteration for adjoint field calculation
            It->run_iteration_adjoint(IP, grid, io);

            // calculate sensitivity kernel
            calculate_sensitivity_kernel(grid, IP);

            if (subdom_main && !line_search_mode && IP.get_is_output_source_field()) {
                // adjoint field will be output only at the end of subiteration
                // output the result of adjoint simulation
                io.write_adjoint_field(grid,i_inv);
            }
       }

        // delete iterator object

    } // end loop sources


    // wait for all processes to finish
    synchronize_all_world();

    // allreduce sum_adj_src
    allreduce_cr_sim_single(v_obj, v_obj);

    // return current objective function value
    return v_obj;

}


#endif // MAIN_ROUTINES_INVERSION_MODE_H