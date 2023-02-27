#ifndef MODEL_OPTIMIZATION_ROUTINES_H
#define MODEL_OPTIMIZATION_ROUTINES_H

#include <iostream>
#include <memory>
#include "mpi_funcs.h"
#include "config.h"
#include "utils.h"
#include "input_params.h"
#include "grid.h"
#include "io.h"
#include "main_routines_calling.h"
#include "iterator_selector.h"
#include "iterator.h"
#include "iterator_legacy.h"
#include "iterator_level.h"
#include "source.h"
#include "receiver.h"
#include "kernel.h"
#include "model_update.h"
#include "lbfgs.h"



// do model update by gradient descent
inline void model_optimize(InputParams& IP, Grid& grid, IO_utils& io, int i_inv, \
                    CUSTOMREAL& v_obj_inout, CUSTOMREAL& old_v_obj, CUSTOMREAL& v_misfit_inout, bool& first_src, std::ofstream& out_main) {

    // change stepsize
    if (i_inv > 0 && v_obj_inout < old_v_obj) {
        // step_size_init = std::min(0.01, step_size_init*1.03);
        step_size_init = std::min(1.0, step_size_init);
        step_size_init_sc = std::min(1.0, step_size_init_sc);
    } else if (i_inv > 0 && v_obj_inout >= old_v_obj) {
        // step_size_init = std::max(0.00001, step_size_init*0.97);
        step_size_init = std::max(0.00001, step_size_init*step_size_decay);
        step_size_init_sc = std::max(0.00001, step_size_init_sc*step_size_decay);
    }
    // output objective function
    if (myrank==0 && id_sim==0) out_main << std::setw(5) << i_inv \
                                  << "," << std::setw(15) << v_obj_inout \
                                  << "," << std::setw(15) << v_misfit_inout \
                                  << "," << std::setw(15) << step_size_init << std::endl;

    // sum kernels among all simultaneous runs
    sumup_kernels(grid);

    // smooth kernels
    smooth_kernels(grid, IP);

    // update the model with the initial step size
    set_new_model(grid, step_size_init);

    // make station correction
    IP.station_correction_update(step_size_init_sc);

    // # TODO: only the first simultanoue run group should output the model. but now ever group outputs the model.
    if (subdom_main && IP.get_is_verbose_output()) {
        // store kernel only in the first src datafile
        io.change_group_name_for_model();

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
    if (IP.get_is_verbose_output())
        io.update_xdmf_file();

    synchronize_all_world();



}


inline void model_optimize_halve_stepping(InputParams& IP, Grid& grid, IO_utils& io, int i_inv, CUSTOMREAL& v_obj_inout, bool& first_src, std::ofstream& out_main) {

    CUSTOMREAL step_size = step_size_init; // step size init is global variable
    CUSTOMREAL diff_obj  = - 9999999999;
    CUSTOMREAL v_obj_old = v_obj_inout;
    CUSTOMREAL v_obj_new = 0.0;
    std::vector<CUSTOMREAL> v_obj_misfit_new = std::vector<CUSTOMREAL>(2);
    int sub_iter_count   = 0;

    // sum kernels among all simultaneous runs
    sumup_kernels(grid);

    // smooth kernels
    smooth_kernels(grid, IP);

    // backup the initial model
    if(subdom_main) grid.back_up_fun_xi_eta_bcf();

    // update the model with the initial step size
    set_new_model(grid, step_size);


    if (subdom_main && IP.get_is_verbose_output()) {
        // store kernel only in the first src datafile
        io.change_group_name_for_model();

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
    if (IP.get_is_verbose_output())
        io.update_xdmf_file();

    synchronize_all_world();

    // loop till
    while (sub_iter_count < max_sub_iterations) {
        // check the new objective function value
        // v_obj_new = run_simulation_one_step(IP, grid, io, i_inv, first_src, true); // run simulations with line search mode
        v_obj_misfit_new = run_simulation_one_step(IP, grid, io, i_inv, first_src, true); // run simulations with line search mode
        v_obj_new = v_obj_misfit_new[0];
        // if the new objective function value is larger than the old one, make the step width to be half of the previous one
        diff_obj = v_obj_new - v_obj_old;

        if (diff_obj > _0_CR) {
            // print status
            if(myrank == 0 && id_sim ==0)
                out_main \
                           << std::setw(5) << i_inv \
                    << "," << std::setw(5) << sub_iter_count \
                    << "," << std::setw(15) << step_size \
                    << "," << std::setw(15) << diff_obj \
                    << "," << std::setw(15) << v_obj_new \
                    << "," << std::setw(15) << v_obj_old << std::endl;

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
        out_main \
                << std::setw(5) << i_inv \
         << "," << std::setw(5) << sub_iter_count \
         << "," << std::setw(15) << step_size \
         << "," << std::setw(15) << diff_obj \
         << "," << std::setw(15) << v_obj_new \
         << "," << std::setw(15) << v_obj_old << " accepted." << std::endl;

    v_obj_inout = v_obj_new;

    // write adjoint field
    int next_i_inv = i_inv + 1;
    if (subdom_main && IP.get_is_output_source_field())
        io.write_adjoint_field(grid,next_i_inv);


}


// do model update
inline void model_optimize_lbfgs(InputParams& IP, Grid& grid, IO_utils& io, int i_inv, CUSTOMREAL& v_obj_inout, bool& first_src, std::ofstream& out_main) {

    int        subiter_count = 0;              // subiteration count
    CUSTOMREAL q_k           = _0_CR;          // store p_k * grad(f_k)
    CUSTOMREAL q_k_new       = _0_CR;          // store p_k * grad(f_k+1)
    CUSTOMREAL v_obj_cur     = v_obj_inout;    // store objective function value
    CUSTOMREAL td            = _0_CR;          // wolfes right step
    CUSTOMREAL tg            = _0_CR;          // wolfes left step
    CUSTOMREAL step_size     = step_size_init; // step size init is global variable
    CUSTOMREAL v_obj_reg     = _0_CR;          // regularization term
    CUSTOMREAL v_obj_new     = v_obj_cur;      // objective function value at new model
    std::vector<CUSTOMREAL> v_obj_misfit_new = std::vector<CUSTOMREAL>(2);

    // smooth kernels and calculate descent direction
    calc_descent_direction(grid, i_inv, IP);

    // smooth descent direction
    //smooth_descent_direction(grid);

    // compute initial q_k for line search = initial_gradient * descent_direction
    q_k = compute_q_k(grid);

    // regularization delta_chi -> delta_chi' = grad + coef * delta L(m)
    // add gradient of regularization term
    add_regularization_grad(grid);
    // add regularization term to objective function
    v_obj_reg = add_regularization_obj(grid);
    v_obj_new += v_obj_reg;

    // store initial gradient
    if (i_inv == 0) {
        store_model_and_gradient(grid, i_inv);
    }

    // backup the initial model
    if(subdom_main) grid.back_up_fun_xi_eta_bcf();

    if (subdom_main && IP.get_is_verbose_output()) {
        // store kernel only in the first src datafile
        io.change_group_name_for_model();

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
    if (IP.get_is_verbose_output())
        io.update_xdmf_file();

    synchronize_all_world();

    //initial_guess_step(grid, step_size, 1.0);
    bool init_bfgs = false;

    // do line search for finding a good step size
    while (optim_method == LBFGS_MODE) {
        // decide initial step size
        if (i_inv == 0 && subiter_count == 0) {
            initial_guess_step(grid, step_size, step_size_init);
            init_bfgs = true;
            step_size_lbfgs=step_size_init; // store input step length
        }
        else if (i_inv == 1 && subiter_count == 0) {
            initial_guess_step(grid, step_size, step_size_lbfgs*LBFGS_RELATIVE_STEP_SIZE);
            //step_size *= 0.1;
        }
        //if (i_inv==0) init_bfgs=true;

        // update the model
        if(subdom_main) grid.restore_fun_xi_eta_bcf();
        set_new_model(grid, step_size, init_bfgs);

        // check current objective function value
        v_obj_misfit_new = run_simulation_one_step(IP, grid, io, i_inv, first_src, true);
        v_obj_new = v_obj_misfit_new[0];

        // update gradient
        sumup_kernels(grid);

        // smooth kernels
        smooth_kernels(grid, IP);

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
            out_main \
                     << std::setw(5)  << i_inv \
              << "," << std::setw(5)  << subiter_count \
              << "," << std::setw(15) << store_step_size \
              << "," << std::setw(15) << (v_obj_new-v_obj_cur)/store_step_size \
              << "," << std::setw(15) << v_obj_new \
              << "," << std::setw(15) << v_obj_reg \
              << "," << std::setw(15) << q_k_new \
              << "," << std::setw(15) << q_k \
              << "," << std::setw(15) << td \
              << "," << std::setw(15) << tg \
              << "," << std::setw(15) << wolfe_c1*q_k \
              << "," << std::setw(15) << wolfe_c2*q_k \
              << "," << std::setw(15) << wolfe_cond_ok << std::endl;

        if (wolfe_cond_ok) {
            // if yes, update the model and break
            v_obj_inout = v_obj_new;

            // store current model and gradient
            store_model_and_gradient(grid, i_inv+1);

            step_size_init = step_size; // use current step size for the next iteration

            goto end_of_subiteration;
        } else if (subiter_count > max_sub_iterations){
            // reached max subiter
            // exit
            std::cout << "reached max subiterations" << std::endl;
            finalize_mpi();
            exit(1);
        } else {
            // wolfe conditions not satisfied

            //v_obj_cur = v_obj_new;
            subiter_count++;
        }
    }

end_of_subiteration:
    if (myrank == 0)
        std::cout << "Wolfe conditions satisfied at iteration " << subiter_count << std::endl;

}



#endif // MODEL_OPTIMIZATION_ROUTINES_H