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

#ifdef USE_CUDA
#include "cuda_initialize.cuh"
#endif

// TOMOATT main function
int main(int argc, char *argv[])
{
    // change flag
    store_tau = true;


    // parse options
    parse_options(argc, argv);

    // here we set i_inv = 0; but later we can modify the code to set from command line
    int i_inv = 0;

    // initialize mpi
    initialize_mpi();

    stdout_by_rank_zero("------------------------------------------------------");
    stdout_by_rank_zero("start TOMOATT data preparation for NN model training.");
    stdout_by_rank_zero("------------------------------------------------------");

    // read input file
    InputParams IP(input_file);

    // create output directory
    create_output_dir(output_dir);

#ifdef USE_CUDA
    // end program
    if(use_gpu) {
        std::cout << "GPU mode is not supported in this version." << std::endl;
        std::cout << "In parameter file, please set use_gpu: 0." << std::endl;
        exit(1);
    }
#endif

#ifdef USE_SIMD
    // check SIMD type
    if (myrank==0)
        print_simd_type();
#endif


    //
    // This is initial try. Thus domain decomposition is not implemented yet.
    //
    if (ndiv_i != 1 || ndiv_j != 1 || ndiv_k != 1) {
        std::cout << "Domain decomposition is not implemented yet." << std::endl;
        exit(1);
    }

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

    // iterate over sources
    for (int i_src = 0; i_src < (int)IP.src_id2name.size(); i_src++){

        const std::string name_sim_src = IP.src_id2name[i_src];
        const int         id_sim_src   = IP.src_map[name_sim_src].id; // global source id

        if (myrank == 0)
            std::cout << "source name: " << name_sim_src << ", forward modeling starting..." << std::endl;

        // set group name to be used for output in h5
        io.set_id_src(id_sim_src);
        io.set_name_src(name_sim_src);
        io.change_group_name_for_source();

        // get is_teleseismic flag
        bool is_teleseismic = IP.get_if_src_teleseismic(name_sim_src);

        // (re) initialize source object and set to grid
        Source src(IP, grid, is_teleseismic, name_sim_src);

        // initialize iterator object
        bool first_init = (i_src==0);

        /////////////////////////
        // run forward simulation
        /////////////////////////

        // initialize iterator object
        std::unique_ptr<Iterator> It;

        select_iterator(IP, grid, src, io, name_sim_src, first_init, is_teleseismic, It, false);
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

    // finalize cuda
//#ifdef USE_CUDA
//    if (use_gpu) finalize_cuda();
//#endif

    // finalize mpi
    finalize_mpi();

    stdout_by_rank_zero("------------------------------------------------------");
    stdout_by_rank_zero("end TOMOATT only traveltime calculation mode.");
    stdout_by_rank_zero("------------------------------------------------------");

    return 0;
}