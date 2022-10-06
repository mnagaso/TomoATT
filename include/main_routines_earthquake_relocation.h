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
        std::cout << "Error: Source relocation mode may run with swap_src_rec = 1" << std::endl;
        exit(1);
    }

    // initialize kernel arrays
    //if (IP.get_run_mode() == DO_INVERSION)
    //    grid.initialize_kernels();

    // reinitialize factors
    grid.reinitialize_abcf();

    ///////////////////////
    // loop for each source
    ///////////////////////

    for (long unsigned int i_src = 0; i_src < IP.src_ids_this_sim.size(); i_src++) {

        // load the global id of this src
        id_sim_src = IP.src_ids_this_sim[i_src]; // local src id to global src id

        io.init_data_output_file(); // initialize data output file
        io.change_xdmf_obj(i_src); // change xmf file for next src

        // check if the source is teleseismic or not
        // because teleseismic source is not supported in this mode
        bool is_teleseismic = IP.get_src_point(id_sim_src).is_teleseismic;

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

        It->run_iteration_forward(IP, grid, io, first_init);

        // writeout travel time field
        if (subdom_main) {
            // output T (result timetable)
            io.write_T(grid, 0);
        }

    }

}


void calculate_gradient_objective_function(InputParams& IP, Grid& grid, IO_utils& io, std::vector<SrcRec>& unique_rec_list){

    // initialize source parameters
    Receiver recs; // here the source is swapped to receiver
    recs.init_vars_src_reloc(IP, unique_rec_list);

    // iterate over all sources for calculating optimal origin time
    for (long unsigned int i_src = 0; i_src < IP.src_ids_this_sim.size(); i_src++) {

        // reset the file name to be read
        io.change_xdmf_obj(i_src);

        // load the global id of this src
        id_sim_src = IP.src_ids_this_sim[i_src]; // local src id to global src id

        // load travel time field on grid.T_loc
        io.read_T(grid);

        // calculate travel time at the actual source location
        recs.calculate_arrival_time(IP, grid);

        // calculate approximated orptimal origin time
        recs.calculate_optimal_origin_time(IP, unique_rec_list);
    }

    // divide optimal origin time by summed weight
    recs.divide_optimal_origin_time_by_summed_weight(IP, unique_rec_list);

    // iterate over all sources for calculating gradient of objective function
    for (long unsigned int i_src = 0; i_src < IP.src_ids_this_sim.size(); i_src++) {

        // reset the file name to be read
        io.change_xdmf_obj(i_src);

        // load the global id of this src
        id_sim_src = IP.src_ids_this_sim[i_src]; // local src id to global src id

        // load travel time field on grid.T_loc
        io.read_T(grid);

        // calculate gradient at the actual source location
        recs.calculate_T_gradient(IP, grid);

        // calculate gradient of objective function
        recs.calculate_grad_obj_src_reloc(IP, unique_rec_list);

    }
}


std::vector<SrcRec> create_unique_rec_list(InputParams& IP) {

    // unique vector of SrcRec objects
    std::vector<SrcRec> unique_rec_list;

    // loop over all sources
    for (long unsigned int i_src = 0; i_src < IP.src_ids_this_sim.size(); i_src++) {
        // load the global id of this src
        id_sim_src = IP.src_ids_this_sim[i_src]; // local src id to global src id

        // get receiver (swapped source) list for this source
        std::vector<SrcRec>& receivers = IP.get_rec_points(id_sim_src);

        // loop over all receivers
        for (auto& rec : receivers) {
            // first element
            if (unique_rec_list.size() == 0){
                rec.id_unique_list = 0;
                unique_rec_list.push_back(rec);
            } else {
                bool if_exists = false;
                for (long unsigned int i_rec = 0; i_rec < unique_rec_list.size(); i_rec++) {
                    if (unique_rec_list[i_rec].id_src == rec.id_src) {
                        if_exists = true;
                        // store the element id on the unique list
                        rec.id_unique_list = i_rec;
                        break;
                    }
                }

                // add element to unique list if the same station name is not exist
                if (!if_exists) {
                    rec.id_unique_list = unique_rec_list.size();
                    unique_rec_list.push_back(rec);
                }
            }
        }

    }

    return unique_rec_list;

}

#endif // MAIN_ROUTINES_EARTHQUAKE_RELOCATION_H