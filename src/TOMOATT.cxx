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
#include "receiver.h"
#include "kernel.h"
#include "model_update.h"
#include "main_routines_calling.h"
#include "eikonal_solver_2d.h"

#ifdef USE_CUDA
#include "cuda_initialize.cuh"
#endif

//#ifdef USE_BLAS
//#include "cblas.h"
//#endif

// TOMOATT main function
int main(int argc, char *argv[])
{
    // parse options
    parse_options(argc, argv);

    // initialize mpi
    initialize_mpi();

    stdout_by_main("------------------------------------------------------");
    stdout_by_main("start TOMOATT forward or inversion calculation.");
    stdout_by_main("------------------------------------------------------");

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

    // output inversion grid file (by main process)
    grid.write_inversion_grid_file();

    // initialize inversion grids (by other process)
    grid.setup_inversion_grids(IP);

    if (subdom_main) {
        // output grid data (grid data is only output in the main simulation)
        io.write_grid(grid);
    }

    // preapre teleseismic boundary conditions (do nothinng if no teleseismic source is defined)
    prepare_teleseismic_boundary_conditions(IP, grid, io);

    synchronize_all_world();

    //
    // run main calculation routines depending on specified mode
    //
    if (IP.get_run_mode() == ONLY_FORWARD || IP.get_run_mode() == DO_INVERSION){
        run_forward_only_or_inversion(IP, grid, io);
    } else if (IP.get_run_mode() == TELESEIS_PREPROCESS) {
        // terminate the program here if the run mode is TELESEIS_PREPROCESS
        // because prepare_teleseismic_boundary_conditions is already called.
    } else if (IP.get_run_mode() == SRC_RELOCATION) {
        run_earthquake_relocation(IP, grid, io);
    } else {
        std::cerr << "Error: invalid run mode is specified." << std::endl;
        exit(1);
    }

    // output final state of the model
    if (IP.get_is_output_final_model()) {
        io.write_final_model(grid, IP);
    }


    // finalize cuda
#ifdef USE_CUDA
    if (use_gpu) finalize_cuda();
#endif

    // finalize mpi
    finalize_mpi();

    stdout_by_main("------------------------------------------------------");
    stdout_by_main("end TOMOATT solver.");
    stdout_by_main("------------------------------------------------------");

    return 0;
}