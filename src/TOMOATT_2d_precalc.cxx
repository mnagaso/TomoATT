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
    stdout_by_rank_zero("start TOMOATT 2d pre calculation mode.");
    stdout_by_rank_zero("------------------------------------------------------");

    // read input file
    InputParams IP(input_file);

    // create output directory
    create_output_dir(output_dir);

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

    // finalize mpi
    finalize_mpi();

    stdout_by_rank_zero("------------------------------------------------------");
    stdout_by_rank_zero("end TOMOATT 2d pre calculation mode.");
    stdout_by_rank_zero("------------------------------------------------------");

    return 0;
}