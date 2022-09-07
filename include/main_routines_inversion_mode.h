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

        if (i_inv == 0 && !line_search_mode)
            io.init_data_output_file(); // initialize data output file

        io.change_xdmf_obj(i_src); // change xmf file for next src

        // output initial field
        if(first_src) {
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

        // get is_teleseismic flag
        bool is_teleseismic = IP.get_src_point(id_sim_src).is_teleseismic;

        // (re) initialize source object and set to grid
        Source src(IP, grid, is_teleseismic);

        // initialize iterator object
        bool first_init = (i_inv == 0 && i_src==0);

        // initialize iterator object
        std::unique_ptr<Iterator> It;
        select_iterator(IP, grid, src, io, first_init, is_teleseismic, It);

        /////////////////////////
        // run forward simulation
        /////////////////////////

        It->run_iteration_forward(IP, grid, io, first_init);

        // output the result of forward simulation
        // ignored for inversion mode.
        if (subdom_main && !line_search_mode) { // && IP.get_run_mode()!=1) {
            // output T0v
            io.write_T0v(grid,i_inv); // initial Timetable
            // output u (true solution)
            if (if_test)
                io.write_u(grid);  // true Timetable
            // output tau
            io.write_tau(grid, i_inv); // calculated deviation
            // output T (result timetable)
            io.write_T(grid, i_inv);
            // output residual (residual = true_solution - result)
            if (if_test)
                io.write_residual(grid); // this will over write the u_loc, so we need to call write_u_h5 first
        }

        // calculate the arrival times at each receivers
        Receiver recs;
        recs.calculate_arrival_time(IP, grid);

        /////////////////////////
        // run adjoint simulation
        /////////////////////////

        if (IP.get_run_mode()==1){
            if (!is_teleseismic) {// for regional source
                // calculate adjoint source
                v_obj += recs.calculate_adjoint_source(IP);
           } else { // for teleseismic source
                // calculate adjoint source
                v_obj += recs.calculate_adjoint_source_teleseismic(IP);
            }

            // run iteration for adjoint field calculation
            It->run_iteration_adjoint(IP, grid, io);

            // calculate sensitivity kernel
            calculate_sensitivity_kernel(grid, IP);

            if (!line_search_mode){
                // adjoint field will be output only at the end of subiteration
                // output the result of adjoint simulation
                if (subdom_main) {
                    // adjoint field
                    io.write_adjoint_field(grid,i_inv);
                }
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