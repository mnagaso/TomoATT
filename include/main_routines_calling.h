#ifndef MAIN_ROUTINES_CALLING_H
#define MAIN_ROUTINES_CALLING_H

#include <iostream>
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
#include "lbfgs.h"

// run forward and adjoint simulation and calculate current objective function value and sensitivity kernel if requested
CUSTOMREAL run_simulation_one_step(InputParams& IP, Grid& grid, IO_utils& io, int i_inv, bool& first_src, bool line_search_mode){

    CUSTOMREAL v_obj = _0_CR;

    // initialize kernel arrays
    if (IP.get_do_inversion()==1)
        grid.initialize_kernels();

    // reinitialize factors
    grid.reinitialize_abcf();

    ///////////////////////
    // loop for each source
    ///////////////////////

    for (long unsigned int i_src = 0; i_src < IP.src_ids_this_sim.size(); i_src++) {

        // load the global id of this src
        id_sim_src = IP.src_ids_this_sim[i_src]; // local src id to global src id

        if (i_inv == 0 && !line_search_mode)
            io.init_data_output_file(); // initialize data output file

        io.change_xdmf_obj(i_src); // change xmf file for next src

        // output initial field
        if(first_src) {
            if (subdom_main) {
                // write true solution
                if (if_test){
                    io.write_true_solution(grid);
                }
                // write initial velocity model
                //io.write_velocity_model_h5(grid);
            }
            first_src = false;
        }

        // get is_teleseismic flag
        bool is_teleseismic = IP.get_src_point(id_sim_src).is_teleseismic;

        // (re) initialize source object and set to grid
        Source src(IP, grid, is_teleseismic);

        // initialize iterator object
        bool first_init = (i_inv == 0 && i_src==0);
        Iterator It(IP, grid, src, io, first_init, is_teleseismic);

        /////////////////////////
        // run forward simulation
        /////////////////////////

        if(!is_teleseismic)
            It.run_iteration_forward(IP, grid, io, first_init);
        else
            It.run_iteration_forward_teleseismic(IP, grid, io, first_init);

        // output the result of forward simulation
        // ignored for inversion mode.
        if (subdom_main && !line_search_mode) { // && IP.get_do_inversion()!=1) {
            // output T0v
            io.write_T0v(grid,i_inv); // initial Timetable
            // output u (true solution)
            if (if_test)
                io.write_u(grid);  // true Timetable
            // output tau
            io.write_tau(grid, i_inv); // calculated deviation
            // output T (result timetable)
            io.write_T(grid, i_inv);
            // output residual (residual = true_solution - result)
            if (if_test)
                io.write_residual(grid); // this will over write the u_loc, so we need to call write_u_h5 first
        }

        // calculate the arrival times at each receivers
        Receiver recs;
        recs.calculate_arrival_time(IP, grid);

        /////////////////////////
        // run adjoint simulation
        /////////////////////////

        if (IP.get_do_inversion()==1){
            if (!is_teleseismic) {// for regional source
                // calculate adjoint source
                v_obj += recs.calculate_adjoint_source(IP);
           } else { // for teleseismic source
                // calculate adjoint source
                v_obj += recs.calculate_adjoint_source_teleseismic(IP);
            }

            // run iteration for adjoint field calculation
            It.run_iteration_adjoint(IP, grid, io);

            // calculate sensitivity kernel
            calculate_sensitivity_kernel(grid, IP);

            if (!line_search_mode){
                // adjoint field will be output only at the end of subiteration
                // output the result of adjoint simulation
                if (subdom_main) {
                    // adjoint field
                    io.write_adjoint_field(grid,i_inv);
                }
            }
       }


    } // end loop sources


    // wait for all processes to finish
    synchronize_all_world();

    // allreduce sum_adj_src
    allreduce_cr_sim_single(v_obj, v_obj);

    // return current objective function value
    return v_obj;

}


// do model update
void model_optimize(InputParams& IP, Grid& grid, IO_utils& io, int i_inv, CUSTOMREAL& v_obj_inout, bool& first_src, std::ofstream& out_main) {

    CUSTOMREAL step_size = step_size_init; // step size init is global variable

    // sum kernels among all simultaneous runs
    sumup_kernels(grid);

    // smooth kernels
    smooth_kernels(grid);

    // update the model with the initial step size
    set_new_model(grid, step_size);


    if (subdom_main) {
        // store kernel only in the first src datafile
        io.change_xdmf_obj(0); // change xmf file for next src

        // output updated velocity models
        io.write_Ks(grid, i_inv);
        io.write_Keta(grid, i_inv);
        io.write_Kxi(grid, i_inv);

        // output descent direction
        io.write_Ks_update(grid, i_inv);
        io.write_Keta_update(grid, i_inv);
        io.write_Kxi_update(grid, i_inv);
    }

    // writeout temporary xdmf file
    io.update_xdmf_file(IP.src_ids_this_sim.size());

    synchronize_all_world();



}


void model_optimize_halve_stepping(InputParams& IP, Grid& grid, IO_utils& io, int i_inv, CUSTOMREAL& v_obj_inout, bool& first_src, std::ofstream& out_main) {

    CUSTOMREAL step_size = step_size_init; // step size init is global variable
    CUSTOMREAL diff_obj = - 9999999999;
    CUSTOMREAL v_obj_old = v_obj_inout;
    CUSTOMREAL v_obj_new;
    int sub_iter_count = 0;

    // sum kernels among all simultaneous runs
    sumup_kernels(grid);

    // smooth kernels
    smooth_kernels(grid);

    // backup the initial model
    if(subdom_main) grid.back_up_fun_xi_eta_bcf();

    // update the model with the initial step size
    set_new_model(grid, step_size);


    if (subdom_main) {
        // store kernel only in the first src datafile
        io.change_xdmf_obj(0); // change xmf file for next src

        // output updated velocity models
        io.write_Ks(grid, i_inv);
        io.write_Keta(grid, i_inv);
        io.write_Kxi(grid, i_inv);

        // output descent direction
        io.write_Ks_update(grid, i_inv);
        io.write_Keta_update(grid, i_inv);
        io.write_Kxi_update(grid, i_inv);
    }

    // writeout temporary xdmf file
    io.update_xdmf_file(IP.src_ids_this_sim.size());

    synchronize_all_world();

    // loop till
    while (sub_iter_count < max_sub_iterations) {
        // check the new objective function value
        v_obj_new = run_simulation_one_step(IP, grid, io, i_inv, first_src, true); // run simulations with line search mode

        // if the new objective function value is larger than the old one, make the step width to be half of the previous one
        diff_obj = v_obj_new - v_obj_old;

        if (diff_obj > _0_CR) {
            // print status
            if(myrank == 0 && id_sim ==0)
                out_main << "iteration: " << i_inv << " subiteration: " << sub_iter_count << " step_size: " << step_size << \
                            " diff_obj: " << diff_obj << " v_obj_new: " << v_obj_new << " v_obj_old: " << v_obj_old << std::endl;

            if (subdom_main) grid.restore_fun_xi_eta_bcf();
            step_size /= _2_CR;
            set_new_model(grid, step_size);

            sub_iter_count++;
        } else {
            // if the new objective function value is smaller than the old one, make the step width to be twice of the previous one
            goto end_of_sub_iteration;
        }

    }

end_of_sub_iteration:
    // out log
    if(myrank == 0 && id_sim ==0)
        out_main << "iteration: " << i_inv << " subiteration: " << sub_iter_count << " step_size: " << step_size << \
                    " diff_obj: " << diff_obj << " v_obj_new: " << v_obj_new << " v_obj_old: " << v_obj_old << " accepted." << std::endl;

    v_obj_inout = v_obj_new;

    // write adjoint field
    int next_i_inv = i_inv + 1;
    if (subdom_main)
        io.write_adjoint_field(grid,next_i_inv);


}


// do model update
void model_optimize_lbfgs(InputParams& IP, Grid& grid, IO_utils& io, int i_inv, CUSTOMREAL& v_obj_inout, bool& first_src, std::ofstream& out_main) {

    int        subiter_count = 0;              // subiteration count
    CUSTOMREAL q_k           = _0_CR;          // store p_k * grad(f_k)
    CUSTOMREAL q_k_new       = _0_CR;          // store p_k * grad(f_k+1)
    CUSTOMREAL v_obj_cur     = v_obj_inout;    // store objective function value
    CUSTOMREAL td            = _0_CR;          // wolfes right step
    CUSTOMREAL tg            = _0_CR;          // wolfes left step
    CUSTOMREAL step_size     = step_size_init; // step size init is global variable
    CUSTOMREAL v_obj_reg     = _0_CR;          // regularization term
    CUSTOMREAL v_obj_new     = v_obj_cur;      // objective function value at new model


    // smooth kernels and calculate descent direction
    calc_descent_direction(grid, i_inv);

    // smooth descent direction
    smooth_descent_direction(grid);

    // regularization delta_chi -> delta_chi' = grad + coef * delta L(m)
    // add gradient of regularization term
    add_regularization_grad(grid);
    // add regularization term to objective function
    v_obj_reg = add_regularization_obj(grid);
    v_obj_new += v_obj_reg;

    // compute initial q_k for line search = current_gradient * descent_direction
    q_k = compute_q_k(grid);

    // backup the initial model
    if(subdom_main) grid.back_up_fun_xi_eta_bcf();

    if (subdom_main) {
        // store kernel only in the first src datafile
        io.change_xdmf_obj(0); // change xmf file for next src

        // output updated velocity models
        io.write_Ks(grid, i_inv);
        io.write_Keta(grid, i_inv);
        io.write_Kxi(grid, i_inv);

        // output descent direction
        io.write_Ks_update(grid, i_inv);
        io.write_Keta_update(grid, i_inv);
        io.write_Kxi_update(grid, i_inv);
    }

    // writeout temporary xdmf file
    io.update_xdmf_file(IP.src_ids_this_sim.size());

    synchronize_all_world();

    // do line search for finding a good step size
    while (optim_method == LBFGS_MODE) {

        // update the model
        if(subdom_main) grid.restore_fun_xi_eta_bcf();
        set_new_model(grid, step_size);

        // check current objective function value #BUG: strange v_obj at the first sub iteration
        v_obj_new = run_simulation_one_step(IP, grid, io, i_inv, first_src, true); // run simulations with line search mode

        // update gradient
        sumup_kernels(grid);

        // smooth kernels
        smooth_kernels(grid);

        // add gradient of regularization term
        add_regularization_grad(grid);
        // add regularization term to objective function
        CUSTOMREAL v_obj_reg = add_regularization_obj(grid);
        v_obj_new += v_obj_reg;

        // calculate q_k_new
        q_k_new = compute_q_k(grid);

        // check if the current step size satisfies the wolfes conditions
        CUSTOMREAL store_step_size = step_size;
        bool wolfe_cond_ok = check_wolfe_cond(grid, v_obj_cur, v_obj_new, q_k, q_k_new, td, tg, step_size);

        // log out
        if(myrank == 0 && id_sim ==0)
            out_main << "iteration: " << i_inv << " subiteration: " << subiter_count << " step_size: " << store_step_size << \
                                                                                        " qpt: " << (v_obj_new-v_obj_cur)/store_step_size << \
                                                                                        " v_obj_new: " << v_obj_new << " v_obj_reg: " << v_obj_reg << \
                                                                                        " q_new: " << q_k_new << " q_k: " << q_k << \
                                                                                        " td, tg: " << td << ", " << tg << \
                                                                                        " wolfe_c1*q_k, wolfe_c2*q_k: " << wolfe_c1*q_k << ", " << wolfe_c2*q_k << \
                                                                                        " step ok: " << wolfe_cond_ok << std::endl;

        if (wolfe_cond_ok) {
            // if yes, update the model and break
            v_obj_inout = v_obj_new;

            // store current model and gradient
            store_model_and_gradient(grid, i_inv+1);
            //step_size_init = step_size; // use the

            goto end_of_subiteration;
        } else if (subiter_count > max_sub_iterations){
            // reached max subiter
            // exit
            std::cout << "reached max subiterations" << std::endl;
            finalize_mpi();
            exit(1);
        } else {
            // wolfe conditions not satisfied

            v_obj_cur = v_obj_new;
            subiter_count++;
        }
    }

end_of_subiteration:
    if (myrank == 0)
        std::cout << "Wolfe conditions satisfied at iteration " << subiter_count << std::endl;

}


#endif // MAIN_ROUTINES_CALLING_H