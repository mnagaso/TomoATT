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

    for (long unsigned int i_src = 0; i_src < IP.src_names_this_sim.size(); i_src++) {

        // load the global id of this src
        id_sim_src = IP.src_ids_this_sim[i_src]; // local src id to global src id
        name_sim_src = IP.src_names_this_sim[i_src]; // local src id to global src id
        

        // set group name to be used for output in h5
        io.change_group_name_for_source();

        // check if the source is teleseismic or not
        // because teleseismic source is not supported in this mode
        bool is_teleseismic = IP.get_src_point_nv(name_sim_src).is_out_of_region;

        if (is_teleseismic){
            std::cout << "Error: Teleseismic source is not supported in source relocation mode." << std::endl;
            exit(1);
        }

        Source src(IP, grid, is_teleseismic);

        // initialize iterator object
        bool first_init = (i_src==0);

        // initialize iterator object
        std::unique_ptr<Iterator> It;
        select_iterator(IP, grid, src, io, first_init, is_teleseismic, It, false);

        /////////////////////////
        // run forward simulation
        /////////////////////////

        std::cout << "forward modeling start, name: " << name_sim_src << ", id: " << id_sim_src << ", lat: " << IP.src_list_nv[name_sim_src].lat 
                  << ", lon: " << IP.src_list_nv[name_sim_src].lon << ", dep: " << IP.src_list_nv[name_sim_src].dep
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
    // recs.init_vars_src_reloc(IP, unique_rec_list);

    // initialize source parameters (obj, kernel, )
    // if(id_sim == 0 && myrank == 0)
    //     std::cout << "ckp1, init_vars_src_reloc ... " << std::endl;
    recs.init_vars_src_reloc(IP);

    // iterate over all sources for calculating optimal origin time
    // for (long unsigned int i_src = 0; i_src < IP.src_ids_this_sim.size(); i_src++) {
    // if(id_sim == 0 && myrank == 0)
    //     std::cout << "ckp2, compute optimal ortime step 1, load time field... " << std::endl;

    for (long unsigned int i_src = 0; i_src < IP.src_names_this_sim.size(); i_src++) {    

        // load the global id of this src
        id_sim_src = IP.src_ids_this_sim[i_src]; // local src id to global src id
        name_sim_src = IP.src_names_this_sim[i_src]; // local src name to global src id
        
        // change target group to be read
        // io.change_group_name_for_source_nv();
        io.change_group_name_for_source();       
        // load travel time field on grid.T_loc
        io.read_T(grid);

        // calculate travel time at the actual source location
        recs.store_arrival_time(IP, grid, name_sim_src);

        // calculate approximated orptimal origin time
        recs.calculate_optimal_origin_time(IP);
    }
    
    // if(id_sim == 0 && myrank == 0)
    //     std::cout << "ckp3, compute optimal ortime step 2 ... " << std::endl;
    // divide optimal origin time by summed weight
    if (IP.is_ortime_local_search == 0) {
        recs.divide_optimal_origin_time_by_summed_weight(IP);
    } else {
        for(int i = 0; i < (int)IP.name_for_reloc.size(); i++){
            std::string name_rec = IP.name_for_reloc[i];
            allreduce_cr_sim_single_inplace(IP.rec_list_nv[name_rec].grad_tau);
        }
        synchronize_all_world();
    }
        
    // compute the objective function
    recs.calculate_obj_reloc(IP, i_iter);

    // print the location info
    // for(auto iter = IP.rec_list_nv.begin(); iter != IP.rec_list_nv.end(); iter++){
    //     std::cout << "optimal_ortime is: " << iter->second.tau_opt << std::endl;
    // }

    // iterate over all sources for calculating gradient of objective function
    // if(id_sim == 0 && myrank == 0)
    //     std::cout << "ckp4, compute gradient ... " << std::endl;

    for (long unsigned int i_src = 0; i_src < IP.src_names_this_sim.size(); i_src++) {
        // load the global id of this src
        id_sim_src = IP.src_ids_this_sim[i_src]; // local src id to global src id
        name_sim_src = IP.src_names_this_sim[i_src]; // local src id to global src id

        // reset the file name to be read
        io.change_group_name_for_source();

        // load travel time field on grid.T_loc
        io.read_T(grid);

        // calculate gradient at the actual source location
        recs.calculate_T_gradient(IP, grid);

        // calculate gradient of objective function
        recs.calculate_grad_obj_src_reloc(IP);

    }

    
    // for(auto iter = IP.rec_list_nv.begin(); iter != IP.rec_list_nv.end(); iter++){
    //     allreduce_cr_sim_single_inplace(iter->second.grad_chi_k);
    //     allreduce_cr_sim_single_inplace(iter->second.grad_chi_j);
    //     allreduce_cr_sim_single_inplace(iter->second.grad_chi_i);
    // }
    // synchronize_all_world();

    // if(id_sim == 0 && myrank == 0)
    //     std::cout << "ckp5, reduce gradient from all sources ... " << std::endl;

    // sum grad_obj_src_reloc of all simulation groups
    for(int i = 0; i < (int)IP.name_for_reloc.size(); i++){
        std::string name_rec = IP.name_for_reloc[i];
        allreduce_cr_sim_single_inplace(IP.rec_list_nv[name_rec].grad_chi_k);
        allreduce_cr_sim_single_inplace(IP.rec_list_nv[name_rec].grad_chi_j);
        allreduce_cr_sim_single_inplace(IP.rec_list_nv[name_rec].grad_chi_i);
    }
    synchronize_all_world();

    // if(id_sim == 0 && myrank == 0){
    //     for(auto iter = IP.rec_list_nv.begin(); iter != IP.rec_list_nv.end(); iter++){
    //         std::cout << "grad_chi_k: " << iter->second.grad_chi_k
    //                   << ", grad_chi_j: " << iter->second.grad_chi_j
    //                   << ", grad_chi_i: " << iter->second.grad_chi_i
    //                   << std::endl;
    //     }
    // }
}


#endif // MAIN_ROUTINES_EARTHQUAKE_RELOCATION_H