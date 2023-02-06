#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include "mpi_funcs.h"
#include "config.h"
#include "utils.h"
#include "input_params.h"
#include "grid.h"
#include "io.h"
#include "source.h"
#include "iterator.h"
#include "iterator_selector.h"
#include "eikonal_solver_2d.h"

#ifdef USE_CUDA
#include "cuda_initialize.cuh"
#endif

// TOMOATT main function
int main(int argc, char *argv[])
{
    // parse options
    parse_options(argc, argv);

    // here we set i_inv = 0; but later we can modify the code to set from command line
    int i_inv = 0;

    // initialize mpi
    initialize_mpi();

    stdout_by_rank_zero("------------------------------------------------------");
    stdout_by_rank_zero("start TOMOATT Kirchhoff migration mode.");
    stdout_by_rank_zero("------------------------------------------------------");

    // read input file
    InputParams IP(input_file);

    // create output directory
    create_output_dir(output_dir);

#ifdef USE_CUDA
    // initialize cuda
    if(use_gpu) initialize_cuda();
#endif

#ifdef USE_SIMD
    // check SIMD type
    if (myrank==0)
        print_simd_type();
#endif

    // check the number of mpi processes and ndiv setting is consistent
    check_total_nprocs_and_ndiv();

    // split mpi communicator for simultaneous run, subdomain division and sweep parallelization
    split_mpi_comm();

    // assign source for each simultaneous run group
    IP.prepare_src_list();

    // initialize file IO object
    IO_utils io(IP); // create IO object for main and non-main process as well

    // initialize compute grids
    Grid grid(IP, io); // member objects are created in only the main process of subdomain groups

    if (subdom_main) {
        // output grid data (grid data is only output in the main simulation)
        io.write_grid(grid);
    }

    // preapre teleseismic boundary conditions (do nothinng if no teleseismic source is defined)
    prepare_teleseismic_boundary_conditions(IP, grid, io);

    synchronize_all_world();

    // initialize factors
    grid.reinitialize_abcf();

    // prepare iteration group in xdmf file
    io.prepare_grid_inv_xdmf(i_inv);

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

        // get is_teleseismic flag
        bool is_teleseismic = IP.get_if_src_teleseismic(id_sim_src);

        // (re) initialize source object and set to grid
        Source src(IP, grid, is_teleseismic);

        // initialize iterator object
        bool first_init = (i_src==0);

        /////////////////////////
        // run forward simulation
        /////////////////////////

        // initialize iterator object
        std::unique_ptr<Iterator> It;

        select_iterator(IP, grid, src, io, first_init, is_teleseismic, It, false);
        It->run_iteration_forward(IP, grid, io, first_init);

        // output the result of forward simulation
        // ignored for inversion mode.
        if (subdom_main){
            // output T (result timetable)
            io.write_T_merged(grid, IP, i_inv);
        }

    } // end loop sources

    // wait for all processes to finish
    synchronize_all_world();

    // close xdmf file
    io.finalize_data_output_file();

    // write src_rec file
    IP.write_src_rec_file(i_inv);

    // wait for all processes to finish
    synchronize_all_world();

    //
    // # TODO : here for backward simulation
    //

    i_inv = 1;

    // swap the src-receiver pair for kirchhoff migration
    IP.prepare_src_list_for_backward_run();

    // prepare output directory for backward run
    output_dir = "OUTPUT_FILES_back";
    create_output_dir(output_dir);

    // initialize file IO object
    IO_utils io_back(IP); // create IO object for main and non-main process as well

    // create grid data file
    if (subdom_main) {
        // output grid data (grid data is only output in the main simulation)
        io_back.write_grid(grid);
    }


    // prepare iteration group in xdmf file
    io_back.prepare_grid_inv_xdmf(i_inv);

    /////////////////////////////////
    // loop for each (swapped/backward) source
    /////////////////////////////////
    for (long unsigned int i_src = 0; i_src < IP.src_ids_this_sim.size(); i_src++) {

        // load the global id of this src
        id_sim_src = IP.src_ids_this_sim[i_src]; // local src id to global src id

        if (myrank == 0)
            std::cout << "source id: " << id_sim_src << ", backward modeling starting..." << std::endl;

        // set group name to be used for output in h5
        io_back.change_group_name_for_source();

        // get is_teleseismic flag
        bool is_teleseismic = IP.get_if_src_teleseismic(id_sim_src);

        // abort if teleseismic source is used for backward simulation
        if (is_teleseismic){
            if (myrank == 0)
                std::cout << "Error: teleseismic source is used for backward simulation. Abort." << std::endl;
            exit(1);
        }

        // (re) initialize source object and set to grid
        Source src(IP, grid, is_teleseismic);

        // initialize iterator object
        bool first_init = (i_src==0);

        /////////////////////////
        // run backward simulation
        /////////////////////////

        // initialize iterator object
        std::unique_ptr<Iterator> It_back;

        select_iterator(IP, grid, src, io_back, first_init, is_teleseismic, It_back, false);
        It_back->run_iteration_forward(IP, grid, io_back, first_init);

        // output the result of backward simulation
        // ignored for inversion mode.
        if (subdom_main){
            // output T (result timetable)
            io_back.write_T_merged(grid, IP, i_inv);
        }

    } // end loop sources


    // close xdmf file
    io_back.finalize_data_output_file();

    // finalize cuda
#ifdef USE_CUDA
    if (use_gpu) finalize_cuda();
#endif

    // finalize mpi
    finalize_mpi();

    stdout_by_rank_zero("------------------------------------------------------");
    stdout_by_rank_zero("end TOMOATT Kirchhoff migration mode.");
    stdout_by_rank_zero("------------------------------------------------------");

    return 0;
}