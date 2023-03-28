#ifndef MAIN_ROUTINES_EARTHQUAKE_RELOCATION_H
#define MAIN_ROUTINES_EARTHQUAKE_RELOCATION_H

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


// calculate traveltime field for all sources-receivers and write to file
void calculate_traveltime_for_all_src_rec(InputParams& IP, Grid& grid, IO_utils& io){

    // check if this run is in source receiver swap mode
    // if not, stop the program
    if (IP.get_is_srcrec_swap() == false){
        std::cout << "Error: Source relocation mode must run with swap_src_rec = 1" << std::endl;
        exit(1);
    }

    // initialize kernel arrays
    //if (IP.get_run_mode() == DO_INVERSION)
    //    grid.initialize_kernels();

    // reinitialize factors
    grid.reinitialize_abcf();

    // prepare inverstion iteration group in xdmf file
    io.prepare_grid_inv_xdmf(0);

    ///////////////////////
    // loop for each source
    ///////////////////////

    // iterate over sources
    for (int i_src = 0; i_src < (int)IP.src_id2name.size(); i_src++){

        const std::string name_sim_src = IP.src_id2name[i_src];
        const int         id_sim_src   = IP.src_map[name_sim_src].id; // global source id

        // set simu group id and source name for output files/dataset names
        io.set_id_src(id_sim_src);
        io.set_name_src(name_sim_src);

        // set group name to be used for output in h5
        io.change_group_name_for_source();

        // check if the source is teleseismic or not
        // because teleseismic source is not supported in this mode
        bool is_teleseismic = IP.get_if_src_teleseismic(name_sim_src);

        if (is_teleseismic){
            std::cout << "Error: Teleseismic source is not supported in source relocation mode." << std::endl;
            exit(1);
        }

        Source src(IP, grid, is_teleseismic, name_sim_src);

        // initialize iterator object
        bool first_init = (i_src==0);

        // initialize iterator object
        std::unique_ptr<Iterator> It;
        select_iterator(IP, grid, src, io, name_sim_src, first_init, is_teleseismic, It, false);

        /////////////////////////
        // run forward simulation
        /////////////////////////

        std::cout << "calculateing source (" << i_src+1 << "/" << (int)IP.src_id2name.size() << "), name: " 
                  << name_sim_src << ", lat: " << IP.src_map[name_sim_src].lat
                  << ", lon: " << IP.src_map[name_sim_src].lon << ", dep: " << IP.src_map[name_sim_src].dep
                  << std::endl;

        It->run_iteration_forward(IP, grid, io, first_init);

        // writeout travel time field
        if (subdom_main) {
            // output T (result timetable)
            io.write_T(grid, 0);
        }
    }
}


void calculate_gradient_objective_function(InputParams& IP, Grid& grid, IO_utils& io, int i_iter){

    Receiver recs; // here the source is swapped to receiver
    // initialize source parameters (obj, kernel, )
    recs.init_vars_src_reloc(IP);

    // iterate over sources
    for (int i_src = 0; i_src < (int)IP.src_id2name.size(); i_src++){

        const std::string name_sim_src = IP.src_id2name[i_src];
        const int         id_sim_src   = IP.src_map[name_sim_src].id; // global source id

        // set simu group id and source name for output files/dataset names
        io.set_id_src(id_sim_src);
        io.set_name_src(name_sim_src);

        // change target group to be read
        io.change_group_name_for_source();
        // io.change_group_name_for_source();

        // load travel time field on grid.T_loc
        io.read_T(grid);

        // calculate travel time at the actual source location
        recs.interpolate_and_store_arrival_times_at_rec_position(IP, grid, name_sim_src);

        // calculate approximated orptimal origin time
        recs.calculate_optimal_origin_time(IP, name_sim_src);

    }

    // divide optimal origin time by summed weight
    if (is_ortime_local_search == 0) {
        recs.divide_optimal_origin_time_by_summed_weight(IP);
    } else {
        // sum grad_tau of all simulation groups
        for(auto iter = IP.rec_map.begin(); iter != IP.rec_map.end(); iter++){
            allreduce_cr_sim_single_inplace(iter->second.grad_tau);
        }
    }
    // compute the objective function
    recs.calculate_obj_reloc(IP, i_iter);

    

    // iterate over sources
    for (int i_src = 0; i_src < (int)IP.src_id2name.size(); i_src++){

        const std::string name_sim_src = IP.src_id2name[i_src];
        const int         id_sim_src   = IP.src_map[name_sim_src].id; // global source id

        // set simu group id and source name for output files/dataset names
        io.set_id_src(id_sim_src);
        io.set_name_src(name_sim_src);

        // reset the file name to be read
        io.change_group_name_for_source();

        // load travel time field on grid.T_loc
        io.read_T(grid);

        // calculate gradient at the actual source location
        recs.calculate_T_gradient(IP, grid, name_sim_src);

        // calculate gradient of objective function
        recs.calculate_grad_obj_src_reloc(IP, name_sim_src);
    }

    // sum grad_chi_k of all simulation groups
    for(auto iter = IP.rec_map.begin(); iter != IP.rec_map.end(); iter++){
        allreduce_cr_sim_single_inplace(iter->second.grad_chi_k);
        allreduce_cr_sim_single_inplace(iter->second.grad_chi_j);
        allreduce_cr_sim_single_inplace(iter->second.grad_chi_i);
    }

    //synchronize_all_world(); // not necessary here because allreduce is already synchronizing communication
}



#endif // MAIN_ROUTINES_EARTHQUAKE_RELOCATION_H