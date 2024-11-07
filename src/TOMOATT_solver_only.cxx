//# TODO: this file need to modify for new srcrec handling routines


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
#include "main_routines_inversion_mode.h"

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
    stdout_by_rank_zero("start TOMOATT only traveltime field calculation mode.");
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
    IP.prepare_src_map();

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
    Source src;

    for (long unsigned int i_src = 0; i_src < IP.n_src_comm_rec_this_sim_group; i_src++) {

        // check if this is the first source
        bool first_init = (i_src==0);

        // check source info
        std::string name_sim_src   = IP.get_src_name_comm(i_src);             // source name
        int         id_sim_src     = IP.get_src_id(name_sim_src);             // global source id
        bool        is_teleseismic = IP.get_if_src_teleseismic(name_sim_src); // get is_teleseismic flag

        // set simu group id and source name for output files/dataset names
        io.reset_source_info(id_sim_src, name_sim_src);

        if (myrank == 0)
            std::cout << "source id: " << id_sim_src << ", forward modeling starting..." << std::endl;

        // (re) initialize source object and set to grid
        src.set_source_position(IP, grid, is_teleseismic, name_sim_src);

        /////////////////////////
        // run forward simulation
        /////////////////////////

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
        io.write_T_merged(grid, IP, i_inv);

    } // end loop sources

    // wait for all processes to finish
    synchronize_all_world();

    // close xdmf file
    io.finalize_data_output_file();

    // finalize cuda
#ifdef USE_CUDA
    if (use_gpu) finalize_cuda();
#endif

    // finalize mpi
    finalize_mpi();

    stdout_by_rank_zero("------------------------------------------------------");
    stdout_by_rank_zero("end TOMOATT only traveltime calculation mode.");
    stdout_by_rank_zero("------------------------------------------------------");

    return 0;
}