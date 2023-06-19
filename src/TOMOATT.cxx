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

#include "timer.h"

//#ifdef USE_BLAS
//#include "cblas.h"
//#endif

// TOMOATT main function
int main(int argc, char *argv[])
{
    // parse options
    parse_options(argc, argv);

Timer timer("initialize_mpi");
    // initialize mpi
    initialize_mpi();
timer.stop_timer();

    stdout_by_rank_zero("------------------------------------------------------");
    stdout_by_rank_zero("start TOMOATT forward or inversion calculation.");
    stdout_by_rank_zero("------------------------------------------------------");

Timer timer2("read input file");
    // read input file
    InputParams IP(input_file);
timer2.stop_timer();

Timer timer3("create output directory");
    // create output directory
    create_output_dir(output_dir);
timer3.stop_timer();

#ifdef USE_CUDA
    // initialize cuda
    if(use_gpu) initialize_cuda();
#endif

#ifdef USE_SIMD
    // check SIMD type
    if (myrank==0)
        print_simd_type();
#endif

Timer timer4("check total number of mpi processes and ndiv");
    // check the number of mpi processes and ndiv setting is consistent
    check_total_nprocs_and_ndiv();
timer4.stop_timer();

Timer timer5("split mpi communicator");
    // split mpi communicator for simultaneous run, subdomain division and sweep parallelization
    split_mpi_comm();
timer5.stop_timer();

Timer timer6("prepare src list");
    // assign source for each simultaneous run group
    IP.prepare_src_map();
timer6.stop_timer();

Timer timer7("io");
    // initialize file IO object
    IO_utils io(IP); // create IO object for main and non-main process as well
timer7.stop_timer();

Timer timer8("grid");
    // initialize compute grids
    Grid grid(IP, io); // member objects are created in only the main process of subdomain groups
timer8.stop_timer();

Timer timer9("write inversion grid file");
    // output inversion grid file (by main process)
    grid.write_inversion_grid_file();
timer9.stop_timer();

Timer timer10("setup inversion grids");
    // initialize inversion grids (by other process)
    grid.setup_inversion_grids(IP);
timer10.stop_timer();

Timer timer11("write grid file");
    if (subdom_main) {
        // output grid data (grid data is only output in the main simulation)
        io.write_grid(grid);
    }
timer11.stop_timer();


    // preapre teleseismic boundary conditions (do nothinng if no teleseismic source is defined)
    prepare_teleseismic_boundary_conditions(IP, grid, io);      // not ready for new version of src rec data

    synchronize_all_world();

Timer timer12("run forward only or inversion");
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
timer12.stop_timer();

Timer timer13("output final model");
    // output final state of the model
    if (IP.get_is_output_final_model()) {
        io.write_final_model(grid, IP);
    }
timer13.stop_timer();


    // finalize cuda
#ifdef USE_CUDA
    if (use_gpu) finalize_cuda();
#endif

Timer timer14("finalize mpi");
    // finalize mpi
    finalize_mpi();
timer14.stop_timer();

    stdout_by_rank_zero("------------------------------------------------------");
    stdout_by_rank_zero("end TOMOATT solver.");
    stdout_by_rank_zero("------------------------------------------------------");

    return 0;
}