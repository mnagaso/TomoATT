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

    // noted that src_map_comm_rec is the subset of src_map
    for (int i_src = 0; i_src < IP.n_src_comm_rec_this_sim_group; i_src++){

        std::string name_sim_src = IP.get_src_name_comm(i_src);
        int         id_sim_src   = IP.get_src_id(name_sim_src); // global source id

        // set simu group id and source name for output files/dataset names
        io.set_id_src(id_sim_src);
        io.set_name_src(name_sim_src);

        // set group name to be used for output in h5
        io.change_group_name_for_source();

        bool is_teleseismic = IP.get_if_src_teleseismic(name_sim_src);

        Source src(IP, grid, is_teleseismic, name_sim_src);

        // initialize iterator object
        bool first_init = (i_inv == 0 && i_src==0);

        // initialize iterator object
        std::unique_ptr<Iterator> It;

        select_iterator(IP, grid, src, io, name_sim_src, first_init, is_teleseismic, It, false);

        if (IP.src_map[name_sim_src].is_T_written_into_file){
            // load travel time field on grid.T_loc
            if (myrank == 0){
                std::cout << "reading source (" << i_src+1 << "/" << (int)IP.src_id2name.size() 
                        << "), for common receiver differntial traveltime. name: "
                        << name_sim_src << ", lat: " << IP.src_map[name_sim_src].lat
                        << ", lon: " << IP.src_map[name_sim_src].lon << ", dep: " << IP.src_map[name_sim_src].dep
                        << std::endl;
            }
            
            io.read_T(grid);
        } else {
            // We need to solve eikonal equation
            if (myrank == 0){
                std::cout << "calculating source (" << i_src+1 << "/" << (int)IP.src_id2name.size() 
                        << "), for common receiver differntial traveltime. name: "
                        << name_sim_src << ", lat: " << IP.src_map[name_sim_src].lat
                        << ", lon: " << IP.src_map[name_sim_src].lon << ", dep: " << IP.src_map[name_sim_src].dep
                        << std::endl;
            }
            // solve travel time field on grid.T_loc
            It->run_iteration_forward(IP, grid, io, first_init);

            // writeout travel time field
            if (subdom_main) {
                // output T (result timetable)
                io.write_T(grid, 0);
            }
            
            IP.src_map[name_sim_src].is_T_written_into_file = true;
        }

        Receiver recs;
        recs.interpolate_and_store_arrival_times_at_rec_position(IP, grid, name_sim_src); 
        // CHS: At this point, all the synthesised arrival times for all the co-located stations are recorded in syn_time_map_sr. When you need to use it later, you can just look it up.

    }

    // wait for all processes to finish
    synchronize_all_world();

    std::cout << "synthetic traveltimes of common receiver data have been prepared." << std::endl;

    // gather all the traveltime to the main process and distribute to all processes
    // for calculating the synthetic common receiver differential traveltime
    IP.gather_traveltimes_and_calc_syn_diff();


}


// run forward and adjoint simulation and calculate current objective function value and sensitivity kernel if requested
inline std::vector<CUSTOMREAL> run_simulation_one_step(InputParams& IP, Grid& grid, IO_utils& io, int i_inv, bool& first_src, bool line_search_mode, bool is_read_time){


    // initialize kernel arrays
    if (IP.get_run_mode() == DO_INVERSION || IP.get_run_mode() == INV_RELOC)
        grid.initialize_kernels();

    // reinitialize factors
    grid.reinitialize_abcf();

    ///////////////////////////////////////////////////////////////////////
    //  compute the synthetic common receiver differential traveltime first
    ///////////////////////////////////////////////////////////////////////

    // prepare synthetic traveltime for all earthquakes, if
    //  1. common receiver data exists;
    //  2. we use common receiver data to update model; (cr + not swap) or (cs + swap)
    //  3. we do inversion
    if ( src_pair_exists &&
              ((IP.get_use_cr() && !IP.get_is_srcrec_swap())
            || (IP.get_use_cs() && IP.get_is_srcrec_swap()) )
        && (IP.get_run_mode() == DO_INVERSION || IP.get_run_mode() == INV_RELOC)){
        pre_run_forward_only(IP, grid, io, i_inv);
    }
    //
    // loop over all sources
    //

     // iterate over sources
    for (int i_src = 0; i_src < IP.n_src_this_sim_group; i_src++){

        const std::string name_sim_src = IP.get_src_name(i_src);
        const int         id_sim_src   = IP.get_src_id(name_sim_src); // global source id

        // set simu group id and source name for output files/dataset names
        io.set_id_src(id_sim_src);
        io.set_name_src(name_sim_src);

        // set group name to be used for output in h5
        // TODO set id_sim_src and name_sim_src with this function
        io.change_group_name_for_source();

        // output initial field
        if(first_src && IP.get_if_output_source_field()) {
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

        // initialize iterator object
        bool first_init = (i_inv == 0 && i_src==0);

        // initialize iterator object
        std::unique_ptr<Iterator> It;

        if (!hybrid_stencil_order){
            select_iterator(IP, grid, src, io, name_sim_src, first_init, is_teleseismic, It, false);
            // // we can choose to read traveltime field from the computed traveltime field if (is_read_time == true and is_teleseismic == false)
            // // if ((is_read_time) && (!is_teleseismic)) {
            // if traveltime field has been wriiten into the file, we choose to read the traveltime data.
            if (IP.src_map[name_sim_src].is_T_written_into_file){
                if (myrank == 0){
                    std::cout << "reading source (" << i_src+1 << "/" << (int)IP.src_id2name.size() 
                            << "), for traveltime field. name: "
                            << name_sim_src << ", lat: " << IP.src_map[name_sim_src].lat
                            << ", lon: " << IP.src_map[name_sim_src].lon << ", dep: " << IP.src_map[name_sim_src].dep
                            << std::endl;
                }
                // load travel time field on grid.T_loc
                io.read_T(grid);
            } else {
                // We need to compute traveltime field, including such cases:
                //  case 1. is_read_time == false;
                //  case 2. is_read_time == true, but we compute teleseismic data;

                if (myrank == 0){
                    std::cout << "calculating source (" << i_src+1 << "/" << (int)IP.src_id2name.size() 
                            << "), for traveltime field. name: "
                            << name_sim_src << ", lat: " << IP.src_map[name_sim_src].lat
                            << ", lon: " << IP.src_map[name_sim_src].lon << ", dep: " << IP.src_map[name_sim_src].dep
                            << std::endl;
                }
                // solve travel time field on grid.T_loc
                It->run_iteration_forward(IP, grid, io, first_init);
            }
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
        if (subdom_main && !line_search_mode && IP.get_if_output_source_field()) {

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

        // if (myrank == 0){
        //     std::cout << "calculating adjoint field, source (" << i_src+1 << "/" << (int)IP.src_id2name.size() << "), name: "
        //             << name_sim_src << ", lat: " << IP.src_map[name_sim_src].lat
        //             << ", lon: " << IP.src_map[name_sim_src].lon << ", dep: " << IP.src_map[name_sim_src].dep
        //             << std::endl;
        // }
        if (IP.get_run_mode()==DO_INVERSION || IP.get_run_mode()==INV_RELOC){

            // calculate adjoint source
            recs.calculate_adjoint_source(IP, name_sim_src);

            // run iteration for adjoint field calculation
            It->run_iteration_adjoint(IP, grid, io);

            // calculate sensitivity kernel
            calculate_sensitivity_kernel(grid, IP, name_sim_src);

            if (subdom_main && !line_search_mode && IP.get_if_output_source_field()) {
                // adjoint field will be output only at the end of subiteration
                // output the result of adjoint simulation
                io.write_adjoint_field(grid,i_inv);
            }

            // check adjoint source
            // if (proc_store_srcrec){
            //     for (auto iter = IP.rec_map.begin(); iter != IP.rec_map.end(); iter++){
            //         std::cout << "rec id: " << iter->second.id << ", rec name: " << iter->second.name << ", adjoint source: " << iter->second.adjoint_source << std::endl;
            //     }
            // }

        } // end if run_mode == DO_INVERSION

        // wait for all processes to finish
        // this should not be called here, for the case that the simultaneous run group has different number of sources
        //synchronize_all_world();

    } // end for i_src

    // synchronize all processes
    synchronize_all_world();

    // gather all the traveltime to the main process and distribute to all processes
    // for calculating the synthetic common receiver differential traveltime
    if ( IP.get_run_mode()==ONLY_FORWARD ||                 // case 1. if we are doing forward modeling, traveltime is not prepared for computing cr_dif data. Now we need to compute it
        (!IP.get_use_cr() && !IP.get_is_srcrec_swap())  ||   // case 2-1, we do inversion, but we do not use cr data (cr + no swap)
        (!IP.get_use_cs() &&  IP.get_is_srcrec_swap())){     // case 2-2, we do inversion, but we do not use cr data (cs +    swap)
        IP.gather_traveltimes_and_calc_syn_diff();
    }

    // compute all residual and obj
    Receiver recs;
    std::vector<CUSTOMREAL> obj_residual = recs.calculate_obj_and_residual(IP);


    // return current objective function value
    return obj_residual;
}


#endif // MAIN_ROUTINES_INVERSION_MODE_H