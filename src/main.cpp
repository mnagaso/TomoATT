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
#include "iterator.h"
#include "source.h"
#include "receiver.h"
#include "kernel.h"
#include "model_update.h"
#include "main_routines_calling.h"
#include "eikonal_solver_2d.h"

#ifdef USE_CUDA
#include "cuda_initialize.cuh"
#endif

int main(int argc, char *argv[])
{
    // parse options
    parse_options(argc, argv);

    // initialize mpi
    initialize_mpi();

    // create output directory
    create_output_dir();

#ifdef USE_CUDA
    // initialize cuda
    if(use_gpu) initialize_cuda();
#endif

    stdout_by_main("------------------------------------------------------");
    stdout_by_main("start TOMOATT solver.");
    stdout_by_main("------------------------------------------------------");

    // read input file
    InputParams IP(input_file);

    // check the number of mpi processes and ndiv setting is consistent
    check_total_nprocs_and_ndiv();

    // split mpi communicator for simultaneous run, subdomain division and sweep parallelization
    split_mpi_comm();

    // assign source for each simultaneous run group
    IP.prepare_src_list();

    // initialize file IO object
    IO_utils io; // create IO object for main and non-main process as well

    // initialize compute grids
    Grid grid(IP, io); // member objects are created in only the main process of subdomain groups

    // initialize inversion grids
    grid.setup_inversion_grids(IP);

    if (subdom_main) {
        // output grid data (grid data is only output in the main simulation)
        io.write_grid(grid);
    }

    // preapre teleseismic boundary conditions
    prepare_teleseismic_boundary_conditions(IP, grid, io);

    synchronize_all_world();

    // for check if the current source is the first source
    bool first_src = true;

    if(myrank == 0)
        std::cout << "size of src_list: " << IP.src_ids_this_sim.size() << std::endl;

    // objective function value output
    std::ofstream out_main;
    if(myrank == 0 && id_sim ==0){
        out_main.open("objective_function.txt");
        out_main << "i_inv, obj_value, step_size" << std::endl;
    }

    /////////////////////
    // loop for inversion
    /////////////////////

    bool line_search_mode = false; // if true, run_simulation_one_step skips adjoint simulation and only calculates objective function value

    // objective function for all src
    CUSTOMREAL v_obj = 0.0, old_v_obj = 0.0;

    for (int i_inv = 0; i_inv < IP.get_max_iter_inv(); i_inv++) {

        old_v_obj = v_obj;

        ///////////////////////////////////////////////////////
        // run (forward and adjoint) simulation for each source
        ///////////////////////////////////////////////////////

        // run forward and adjoint simulation and calculate current objective function value and sensitivity kernel for all sources
        line_search_mode = false;
        // skip for the  mode with sub-iteration
        if (i_inv > 0 && optim_method != GRADIENT_DESCENT) {
        } else {
            v_obj = run_simulation_one_step(IP, grid, io, i_inv, first_src, line_search_mode);
        }

        // wait for all processes to finish
        synchronize_all_world();

        // output src rec file with the result arrival times
        IP.write_src_rec_file();

        // change stepsize
        if (optim_method == GRADIENT_DESCENT) {
            if (i_inv > 0 && v_obj < old_v_obj)
                step_size_init = std::min(0.01, step_size_init*1.03);
            else if (i_inv > 0 && v_obj >= old_v_obj)
                step_size_init = std::max(0.00001, step_size_init*0.97);
        }

        // output objective function
        if (myrank==0 && id_sim==0) out_main << i_inv << ", " << v_obj << ", " << step_size_init << std::endl;

        ///////////////
        // model update
        ///////////////

        if (IP.get_do_inversion()==1){
            if (optim_method == GRADIENT_DESCENT)
                model_optimize(IP, grid, io, i_inv, v_obj, first_src, out_main);
            else if (optim_method == LBFGS_MODE)
                model_optimize_lbfgs(IP, grid, io, i_inv, v_obj, first_src, out_main);
            else if (optim_method == HALVE_STEPPING_MODE)
                model_optimize_halve_stepping(IP, grid, io, i_inv, v_obj, first_src, out_main);
        }

        // output updated model
        if (subdom_main && id_sim==0) {
            io.change_xdmf_obj(0); // change xmf file for next src

            // write out model info
            io.write_fun(grid, i_inv);
            io.write_xi(grid, i_inv);
            io.write_eta(grid, i_inv);
            io.write_b(grid, i_inv);
            io.write_c(grid, i_inv);
            io.write_f(grid, i_inv);
        }

        // writeout temporary xdmf file
        io.update_xdmf_file(IP.src_ids_this_sim.size());

    } // end loop inverse

    // close xdmf file
    io.finalize_data_output_file(IP.src_ids_this_sim.size());

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