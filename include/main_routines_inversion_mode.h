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



inline void pre_run_forward_only(InputParams& IP, Grid& grid, IO_utils& io, int i_inv){

    if (IP.src_map_comm_src.size() == 0)
        return;

    int i_src = 0;

    for (auto iter = IP.src_map_comm_src.begin(); iter != IP.src_map_comm_src.end(); iter++){

        int         id_sim_src   = iter->second.id;
        std::string name_sim_src = iter->second.name;

        // check if this source is common receiver data

        if (myrank == 0)
            std::cout << "source id: " << id_sim_src << ", name: " << name_sim_src << ", computing common receiver differential traveltime starting..." << std::endl;

        bool is_teleseismic = IP.get_src_point(name_sim_src).is_out_of_region;

        Source src(IP, grid, is_teleseismic, name_sim_src);

        // initialize iterator object
        bool first_init = (i_inv == 0 && i_src==0);

        // initialize iterator object
        std::unique_ptr<Iterator> It;

        select_iterator(IP, grid, src, io, name_sim_src, first_init, is_teleseismic, It, false);
        It->run_iteration_forward(IP, grid, io, first_init);

        Receiver recs;
        recs.interpolate_and_store_arrival_times_at_rec_position(IP, grid, name_sim_src); // CHS: At this point, all the synthesised arrival times for all the co-located stations are recorded in syn_time_map_sr. When you need to use it later, you can just look it up.

        i_src++;
    }

    // wait for all processes to finish
    synchronize_all_world();

    // all processors share the calculated synthetic traveltime // MNMN this causes a dead lock when the number of sources are different among simultaneous runs
    //IP.reduce_syn_time_list(); // is it necessary?
//    synchronize_all_world(); // dead lock

    //
    std::cout << "synthetic traveltimes of common receiver data have been prepared." << std::endl;


    // check syn_time_map_sr
    // for(auto iter1 = IP.syn_time_map_sr.begin(); iter1 != IP.syn_time_map_sr.end(); iter1++){
    //     for (auto iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++){
    //         std::cout << "id_sim: " << id_sim << ", src name: " << iter1->first << ", rec name: " << iter2->first
    //                   << ", syn time is: " << iter2->second
    //                   << std::endl;
    //     }
    // }

}


inline void run_forward_and_adjoint(){

}


// run forward and adjoint simulation and calculate current objective function value and sensitivity kernel if requested
inline std::vector<CUSTOMREAL> run_simulation_one_step(InputParams& IP, Grid& grid, IO_utils& io, int i_inv, bool& first_src, bool line_search_mode){

    // obj, obj_abs, obj_cs_dif, obj_cr_dif, obj_tele, misfit, misfit_abs, misfit_cs_dif, misfit_cr_dif, misfit_tele
    std::vector<CUSTOMREAL> v_obj_misfit = std::vector<CUSTOMREAL>(20, 0.0);

    // initialize kernel arrays
    if (IP.get_run_mode() == DO_INVERSION)
        grid.initialize_kernels();

    // reinitialize factors
    grid.reinitialize_abcf();

    ///////////////////////////////////////////////////////////////////////
    //  compute the synthetic common receiver differential traveltime first
    ///////////////////////////////////////////////////////////////////////

    // let syn_time_list[name_src][name_rec] = 0.0
    //IP.initialize_syn_time();

    //prepare synthetic traveltime for all earthquakes
    pre_run_forward_only(IP, grid, io, i_inv);

    //
    // loop over all sources
    //
    int i_src = 0;

    // iterate IP.src_map
    for (auto iter = IP.src_map.begin(); iter != IP.src_map.end(); iter++){

        // load the global id of this src
        const int id_sim_src           = iter->second.id;
        const std::string name_sim_src = iter->first;

        // set simu group id and source name for output files/dataset names
        io.set_id_src(id_sim_src);
        io.set_name_src(name_sim_src);

        if (myrank == 0)
            std::cout << "source id: " << id_sim_src << ", name: " << name_sim_src << ", forward modeling starting..." << std::endl;

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

        /////////////////////////
        // run forward simulation
        /////////////////////////

        // get is_teleseismic flag
        bool is_teleseismic = IP.get_if_src_teleseismic(name_sim_src);

        // (re) initialize source object and set to grid
        Source src(IP, grid, is_teleseismic, name_sim_src);

        synchronize_all_world();

        // initialize iterator object
        bool first_init = (i_inv == 0 && i_src==0);

        // initialize iterator object
        std::unique_ptr<Iterator> It;

        if (!hybrid_stencil_order){
            select_iterator(IP, grid, src, io, name_sim_src, first_init, is_teleseismic, It, false);
            It->run_iteration_forward(IP, grid, io, first_init);
        } else {
            // hybrid stencil mode
            std::cout << "\nrunnning in hybrid stencil mode\n" << std::endl;

            // run 1st order forward simulation
            std::unique_ptr<Iterator> It_pre;
            IP.set_stencil_order(1);
            IP.set_conv_tol(IP.get_conv_tol()*100.0);
            select_iterator(IP, grid, src, io, name_sim_src, first_init, is_teleseismic, It_pre, false);
            It_pre->run_iteration_forward(IP, grid, io, first_init);

            // run 3rd order forward simulation
            IP.set_stencil_order(3);
            IP.set_conv_tol(IP.get_conv_tol()/100.0);
            select_iterator(IP, grid, src, io, name_sim_src, first_init, is_teleseismic, It, true);
            It->run_iteration_forward(IP, grid, io, first_init);
        }

        // output the result of forward simulation
        // ignored for inversion mode.
        if (subdom_main && !line_search_mode && IP.get_is_output_source_field()) {

            // output T (result timetable)
            io.write_T(grid, i_inv);

            // output T0v
            //io.write_T0v(grid,i_inv); // initial Timetable
            // output u (true solution)
            //if (if_test)
            //    io.write_u(grid);  // true Timetable
            // output tau
            //io.write_tau(grid, i_inv); // calculated deviation
            // output residual (residual = true_solution - result)
            //if (if_test)
            //    io.write_residual(grid); // this will over write the u_loc, so we need to call write_u_h5 first
        }

        // calculate the arrival times at each receivers
        Receiver recs;
        recs.interpolate_and_store_arrival_times_at_rec_position(IP, grid, name_sim_src);

        /////////////////////////
        // run adjoint simulation
        /////////////////////////

        if (IP.get_run_mode()==DO_INVERSION){

            // calculate adjoint source
            std::vector<CUSTOMREAL> v_obj_misfit_event = recs.calculate_adjoint_source(IP, name_sim_src);
            for(int i = 0; i < (int)v_obj_misfit_event.size(); i++){
                v_obj_misfit[i] += v_obj_misfit_event[i];
            }
            // run iteration for adjoint field calculation
            It->run_iteration_adjoint(IP, grid, io);

            // calculate sensitivity kernel
            calculate_sensitivity_kernel(grid, IP, name_sim_src);

            if (subdom_main && !line_search_mode && IP.get_is_output_source_field()) {
                // adjoint field will be output only at the end of subiteration
                // output the result of adjoint simulation
                io.write_adjoint_field(grid,i_inv);
            }

            // check adjoint source
            // for (auto iter = IP.rec_map_nv.begin(); iter != IP.rec_map_nv.end(); iter++){
            //     std::cout << "rec id: " << iter->second.id << ", rec name: " << iter->second.name << ", adjoint source: " << iter->second.adjoint_source << std::endl;
            // }

        } // end if run_mode == DO_INVERSION

        // wait for all processes to finish
        synchronize_all_world();

        // allreduce sum_adj_src
        for(int i = 0; i < (int)v_obj_misfit.size(); i++){
            allreduce_cr_sim_single_inplace(v_obj_misfit[i]);
        }

    } // end for i_src

    // synchronize all processes
    synchronize_all_world();

    // return current objective function value
    return v_obj_misfit;
}


#endif // MAIN_ROUTINES_INVERSION_MODE_H