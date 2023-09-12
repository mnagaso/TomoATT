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
                    CUSTOMREAL& v_obj_inout, CUSTOMREAL& old_v_obj, bool& first_src, std::ofstream& out_main) {

    // change stepsize
    if (i_inv > 0 && v_obj_inout < old_v_obj) {
        // step_length_init = std::min(0.01, step_length_init*1.03);
        step_length_init    = std::min((CUSTOMREAL)1.0, step_length_init);
        step_length_init_sc = std::min((CUSTOMREAL)1.0, step_length_init_sc);
    } else if (i_inv > 0 && v_obj_inout >= old_v_obj) {
        // step_length_init = std::max(0.00001, step_length_init*0.97);
        step_length_init    = std::max((CUSTOMREAL)0.00001, step_length_init*step_length_decay);
        step_length_init_sc = std::max((CUSTOMREAL)0.00001, step_length_init_sc*step_length_decay);
    }

    // sum kernels among all simultaneous runs
    sumup_kernels(grid);

    // smooth kernels
    smooth_kernels(grid, IP);

    // update the model with the initial step size
    set_new_model(grid, step_length_init);

    // make station correction
    IP.station_correction_update(step_length_init_sc);

    // # TODO: only the first simultanoue run group should output the model. but now ever group outputs the model.
    if (subdom_main && IP.get_verbose_output_level()) {
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
    if (IP.get_verbose_output_level())
        io.update_xdmf_file();

    synchronize_all_world();
}


inline void write_objective_function(InputParams& IP, int i_inv, std::vector<CUSTOMREAL>& v_misfit_inout, std::ofstream& out_main, std::string type) {
    // output objective function
    if (myrank==0 && id_sim==0) {
        out_main << std::setw(5) << i_inv << ",";
        out_main << std::setw(13) << type;
        out_main << "," << std::setw(19) << v_misfit_inout[0];
        // if( IP.data_type.find("abs") != IP.data_type.end())
        out_main << "," << std::setw(19) << v_misfit_inout[1];
        // if( IP.data_type.find("cs_dif") != IP.data_type.end()){
        if (IP.get_is_srcrec_swap())
            out_main << "," << std::setw(19) << v_misfit_inout[3];
        else
            out_main << "," << std::setw(19) << v_misfit_inout[2];
        // }
        // if( IP.data_type.find("cr_dif") != IP.data_type.end()){
        if (IP.get_is_srcrec_swap())
            out_main << "," << std::setw(19) << v_misfit_inout[2];
        else
            out_main << "," << std::setw(19) << v_misfit_inout[3];
        // }
        // if( IP.data_type.find("tele") != IP.data_type.end()){
        out_main << "," << std::setw(19) << v_misfit_inout[4];
        // }
        // res
        CUSTOMREAL mean;
        CUSTOMREAL std;
        std::string tmp;
        if (IP.N_data > 0) {
            mean = v_misfit_inout[5]/IP.N_data;
            std  = sqrt(v_misfit_inout[6]/IP.N_data - my_square(mean));
            tmp = std::to_string(mean);
            tmp.append("/");
            tmp.append(std::to_string(std));
            out_main << "," << std::setw(24) << tmp;
            // std::cout   << ", v_misfit_inout[5]: " << v_misfit_inout[5]
            //                     << ", v_misfit_inout[6]: " << v_misfit_inout[6]
            //                     << ", N_data: " << IP.N_data
            //                     << std::endl;
        } else {
            out_main << "," << std::setw(24) << "0.0/0.0";
        }
        // res_abs
        if (IP.N_abs_local_data > 0) {
            mean = v_misfit_inout[7]/IP.N_abs_local_data;
            std  = sqrt(v_misfit_inout[8]/IP.N_abs_local_data - my_square(mean));
            tmp = std::to_string(mean);
            tmp.append("/");
            tmp.append(std::to_string(std));
            out_main << "," << std::setw(24) << tmp;
        } else {
            out_main << "," << std::setw(24) << "0.0/0.0";
        }
        if (IP.get_is_srcrec_swap()){
            if (IP.N_cr_dif_local_data > 0) {
                mean = v_misfit_inout[11]/IP.N_cr_dif_local_data;
                std  = sqrt(v_misfit_inout[12]/IP.N_cr_dif_local_data - my_square(mean));
                tmp = std::to_string(mean);
                tmp.append("/");
                tmp.append(std::to_string(std));
                out_main << "," << std::setw(24) << tmp;
            } else {
                out_main << "," << std::setw(24) << "0.0/0.0";
            }
            if (IP.N_cs_dif_local_data > 0) {
                mean = v_misfit_inout[9]/IP.N_cs_dif_local_data;
                std  = sqrt(v_misfit_inout[10]/IP.N_cs_dif_local_data - my_square(mean));
                tmp = std::to_string(mean);
                tmp.append("/");
                tmp.append(std::to_string(std));
                out_main << "," << std::setw(24) << tmp;
            } else {
                out_main << "," << std::setw(24) << "0.0/0.0";
            }
        } else {
            if (IP.N_cs_dif_local_data > 0) {
                mean = v_misfit_inout[9]/IP.N_cs_dif_local_data;
                std  = sqrt(v_misfit_inout[10]/IP.N_cs_dif_local_data - my_square(mean));
                tmp = std::to_string(mean);
                tmp.append("/");
                tmp.append(std::to_string(std));
                out_main << "," << std::setw(24) << tmp;
            } else {
                out_main << "," << std::setw(24) << "0.0/0.0";
            }
            if (IP.N_cr_dif_local_data > 0) {
                mean = v_misfit_inout[11]/IP.N_cr_dif_local_data;
                std  = sqrt(v_misfit_inout[12]/IP.N_cr_dif_local_data - my_square(mean));
                tmp = std::to_string(mean);
                tmp.append("/");
                tmp.append(std::to_string(std));
                out_main << "," << std::setw(24) << tmp;
            } else {
                out_main << "," << std::setw(24) << "0.0/0.0";
            }
        }

        if (IP.N_teleseismic_data > 0) {
            mean = v_misfit_inout[13]/IP.N_teleseismic_data;
            std  = sqrt(v_misfit_inout[14]/IP.N_teleseismic_data - my_square(mean));
            tmp = std::to_string(mean);
            tmp.append("/");
            tmp.append(std::to_string(std));
            out_main << "," << std::setw(24) << tmp;
        } else {
            out_main << "," << std::setw(24) << "0.0/0.0";
        }
        out_main << "," << std::setw(19) << step_length_init << "," << std::endl;
    }
}


inline void model_optimize_halve_stepping(InputParams& IP, Grid& grid, IO_utils& io, int i_inv, CUSTOMREAL& v_obj_inout, bool& first_src, std::ofstream& out_main) {

    CUSTOMREAL step_length = step_length_init; // step size init is global variable
    CUSTOMREAL diff_obj  = - 9999999999;
    CUSTOMREAL v_obj_old = v_obj_inout;
    CUSTOMREAL v_obj_new = 0.0;
    std::vector<CUSTOMREAL> v_obj_misfit_new = std::vector<CUSTOMREAL>(10, 0.0);
    int sub_iter_count   = 0;

    // sum kernels among all simultaneous runs
    sumup_kernels(grid);

    // smooth kernels
    smooth_kernels(grid, IP);

    // backup the initial model
    if(subdom_main) grid.back_up_fun_xi_eta_bcf();

    // update the model with the initial step size
    set_new_model(grid, step_length);


    if (subdom_main && IP.get_verbose_output_level()) {
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
    if (IP.get_verbose_output_level())
        io.update_xdmf_file();

    synchronize_all_world();

    // loop till
    while (sub_iter_count < max_sub_iterations) {
        // check the new objective function value
        // v_obj_new = run_simulation_one_step(IP, grid, io, i_inv, first_src, true); // run simulations with line search mode
        bool is_read_time = false;
        v_obj_misfit_new = run_simulation_one_step(IP, grid, io, i_inv, first_src, true, is_read_time); // run simulations with line search mode
        v_obj_new = v_obj_misfit_new[0];
        // if the new objective function value is larger than the old one, make the step width to be half of the previous one
        diff_obj = v_obj_new - v_obj_old;

        if ( diff_obj > _0_CR // if the objective function value is larger than the old one
          || std::abs(diff_obj/v_obj_old) > MAX_DIFF_RATIO_VOBJ) { // if the objective function reduced too much  ) {
            // print status
            if(myrank == 0 && id_sim ==0)
                out_main \
                           << std::setw(5) << i_inv \
                    << "," << std::setw(5) << sub_iter_count \
                    << "," << std::setw(15) << step_length \
                    << "," << std::setw(15) << diff_obj \
                    << "," << std::setw(15) << v_obj_new \
                    << "," << std::setw(15) << v_obj_old << std::endl;

            if (subdom_main) grid.restore_fun_xi_eta_bcf();
            step_length *= HALVE_STEP_RATIO;
            set_new_model(grid, step_length);

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
         << "," << std::setw(15) << step_length \
         << "," << std::setw(15) << diff_obj \
         << "," << std::setw(15) << v_obj_new \
         << "," << std::setw(15) << v_obj_old << " accepted." << std::endl;

    v_obj_inout = v_obj_new;
    step_length_init = step_length/(HALVE_STEP_RATIO)*HALVE_STEP_RESTORAION_RATIO; // use this step size for the next inversion

    // write adjoint field
    int next_i_inv = i_inv + 1;
    if (subdom_main && IP.get_if_output_source_field())
        io.write_adjoint_field(grid,next_i_inv);


}


// do model update
inline bool model_optimize_lbfgs(InputParams& IP, Grid& grid, IO_utils& io, int i_inv, CUSTOMREAL& v_obj_inout, bool& first_src, std::ofstream& out_main) {

    int        subiter_count = 0;              // subiteration count
    //CUSTOMREAL qp_0          = _0_CR;          // store initial p_k * grad(f_k) (descent_direction * gradient)
    CUSTOMREAL qp_t          = _0_CR;          // store p_k * grad(f_k)
    CUSTOMREAL td            = _0_CR;          // wolfes right step
    CUSTOMREAL tg            = _0_CR;          // wolfes left step
    CUSTOMREAL step_length   = step_length_init; // step size init is global variable
    CUSTOMREAL v_obj_reg     = _0_CR;          // regularization term == q_t
    //CUSTOMREAL q_0           = _0_CR;          // initial cost function
    CUSTOMREAL q_t           = _0_CR;          // current cost function
    std::vector<CUSTOMREAL> v_obj_misfit_new = std::vector<CUSTOMREAL>(10, 0.0);

    //
    // initialize
    //
    if (i_inv == 0) {
        // initialize regularization
        // regularization_penalty[:] = 1.0
        init_regulaization_penalty_with_one(grid);

        // calculate volume_domain = SquaredL2Norm(regularization_penalty)
        volume_domain = compute_volume_domain(grid);

        // volume_domain /= N_params
        volume_domain /= N_params;

        // weight_Tikonov = penalty_weight / volume_domain
        weight_Tikonov = regularization_weight / volume_domain;
        weight_Tikonov_ani = regularization_weight;
        std::cout << "DEBUG: weight_Tikonov: " << weight_Tikonov << std::endl;

        // calculate damp_weights
        compute_damp_weights(grid);

        // calculate gradient (run onestep forward+adjoint) not necessary to do here. already done
        //v_obj_misfit_new = run_simulation_one_step(IP, grid, io, i_inv, first_src, true);

        // Q0 = initial obj
        q_0 = v_obj_inout;

        // sum kernels among all simultaneous runs
        sumup_kernels(grid);

        // smooth gradient
        smooth_kernels(grid, IP);

        // calc regularization penalties
        calculate_regularization_penalty(grid);

        //// smooth calculated grad regularization penalty
        //smooth_gradient_regularization(grid);

        // calc (update) obj_regl = squaredL2Norm(regularization_penalty)
        v_obj_reg = calculate_regularization_obj(grid);

        // Q0 += 0.5*weight_Tiknov*obj_regl
        q_0 += v_obj_reg;

        // init_grad += weight_Tikonov * gadient_regularization_penalty
        add_regularization_grad(grid);

        // store model and gradient
        store_model_and_gradient(grid, i_inv);

        // log out
        if(myrank == 0 && id_sim ==0) {
            out_main \
                     << std::setw(5)  << i_inv \
              << "," << std::setw(5)  << subiter_count \
              << "," << std::setw(15) << step_length \
              << "," << std::setw(15) << q_0 \
              << "," << std::setw(15) << q_t \
              << "," << std::setw(15) << v_obj_reg \
              << "," << std::setw(15) << qp_0 \
              << "," << std::setw(15) << qp_t \
              << "," << std::setw(15) << td \
              << "," << std::setw(15) << tg \
              << "," << std::setw(15) << wolfe_c1*qp_0 \
              << "," << std::setw(15) << wolfe_c2*qp_0 \
              << "," << std::setw(15) << "init" << std::endl;
        }
    }

    // backup the initial model
    if(subdom_main) grid.back_up_fun_xi_eta_bcf();

    // update the model with the initial step size
    if (subdom_main && IP.get_verbose_output_level() && id_sim==0) {
        // store kernel only in the first src datafile
        io.change_group_name_for_model();

        // output updated velocity models
        io.write_Ks(grid, i_inv);
        io.write_Keta(grid, i_inv);
        io.write_Kxi(grid, i_inv);

        // output descent direction
        io.write_Ks_descent_dir(grid, i_inv);
        io.write_Keta_descent_dir(grid, i_inv);
        io.write_Kxi_descent_dir(grid, i_inv);
    }

    // writeout temporary xdmf file
    if (IP.get_verbose_output_level()&& id_sim==0)
        io.update_xdmf_file();

    //
    // iteration
    //
    if (i_inv > 0) {
        // sum kernels among all simultaneous runs
        sumup_kernels(grid);

        // smooth kernels
        smooth_kernels(grid, IP);
    }

    // compute descent direction
    calc_descent_direction(grid, i_inv, IP);

    // smooth descent direction
    smooth_descent_direction(grid);

    // calc qp_0 = inital grad * descent direction
    qp_0 = compute_q_k(grid);

    // loop wolfes subiterations
    bool wolfe_cond_ok = false;
    while (wolfe_cond_ok != true) {

        //// initial guess for isub = 0
        if (i_inv == 0 && subiter_count == 0)
            initial_guess_step(grid, step_length, 1.0);
        else if (i_inv == 1 && subiter_count == 0)
            initial_guess_step(grid, step_length, LBFGS_RELATIVE_step_length);

        // log out
        //if(myrank == 0 && id_sim ==0)
        //    out_main \
        //             << std::setw(5)  << i_inv \
        //      << "," << std::setw(5)  << subiter_count \
        //      << "," << std::setw(15) << step_length \
        //      << "," << std::setw(15) << q_0 \
        //      << "," << std::setw(15) << q_t \
        //      << "," << std::setw(15) << v_obj_reg \
        //      << "," << std::setw(15) << qp_0 \
        //      << "," << std::setw(15) << qp_t \
        //      << "," << std::setw(15) << td \
        //      << "," << std::setw(15) << tg \
        //      << "," << std::setw(15) << wolfe_c1*qp_0 \
        //      << "," << std::setw(15) << wolfe_c2*qp_0 \
        //      << "," << std::setw(15) << wolfe_cond_ok << std::endl;

        //std::cout << "DEBUG STEPLENGTH: " << step_length << std::endl;

        //// Update model
        set_new_model(grid, step_length);

        //// calculate gradient (run onestep forward+adjoint)
        v_obj_misfit_new = run_simulation_one_step(IP, grid, io, i_inv, first_src, true, false);

        //// Qt update (=current obj)
        q_t = v_obj_misfit_new[0];

        // sum kernels among all simultaneous runs
        sumup_kernels(grid);

        //// smooth gradient
        smooth_kernels(grid, IP);

        //// calculate regularization to grad
        calculate_regularization_penalty(grid);

        //// smooth calculated grad regularization penalty
        //smooth_gradient_regularization(grid);

        //// add regularization term to Qt (Qt += 0.5 * obj_regl )
        v_obj_reg = calculate_regularization_obj(grid);
        q_t += v_obj_reg;

        //// current grad += weight_Tikonov * gadient_regularization_penalty
        add_regularization_grad(grid);

        //// Qpt = Inner product(current_grad * descent_direction)
        qp_t = compute_q_k(grid);

        //// check wolfe conditions
        wolfe_cond_ok = check_wolfe_cond(grid, q_0, q_t, qp_0, qp_t, td, tg, step_length);

        // update the model with the initial step size
        if (subdom_main && IP.get_verbose_output_level() && id_sim==0) {
            // store kernel only in the first src datafile
            io.change_group_name_for_model();

            // gradient
            io.write_Ks_update(grid, i_inv*100+subiter_count);
            io.write_Keta_update(grid, i_inv*100+subiter_count);
            io.write_Kxi_update(grid, i_inv*100+subiter_count);

            // output descent direction
            io.write_Ks_descent_dir(grid, i_inv*100+subiter_count);
            io.write_Keta_descent_dir(grid, i_inv*100+subiter_count);
            io.write_Kxi_descent_dir(grid, i_inv*100+subiter_count);

        }

        // writeout temporary xdmf file
        if (IP.get_verbose_output_level()&& id_sim==0)
            io.update_xdmf_file();

        if (wolfe_cond_ok) {
            // log out
            if(myrank == 0 && id_sim ==0)
                out_main \
                         << std::setw(5)  << i_inv \
                  << "," << std::setw(5)  << subiter_count \
                  << "," << std::setw(15) << step_length \
                  << "," << std::setw(15) << q_0 \
                  << "," << std::setw(15) << q_t \
                  << "," << std::setw(15) << v_obj_reg \
                  << "," << std::setw(15) << qp_0 \
                  << "," << std::setw(15) << qp_t \
                  << "," << std::setw(15) << td \
                  << "," << std::setw(15) << tg \
                  << "," << std::setw(15) << wolfe_c1*qp_0 \
                  << "," << std::setw(15) << wolfe_c2*qp_0 \
                  << "," << std::setw(15) << wolfe_cond_ok << std::endl;

            goto end_of_subiteration;
        } else if (subiter_count > max_sub_iterations){
            // reached max subiter
            // exit
            std::cout << "reached max subiterations" << std::endl;

            return false;

        } else {
            // wolfe conditions not satisfied
            subiter_count++;
            if(subdom_main) grid.restore_fun_xi_eta_bcf();
        }

    }

end_of_subiteration:
    //// store model and gradient,
    store_model_and_gradient(grid, i_inv+1);

    //// Q0=Qt
    q_0 = q_t;
    qp_0 = qp_t;

    step_length_init = step_length; // use current step size for the next iteration


    if (myrank == 0)
        std::cout << "Wolfe conditions satisfied at iteration " << subiter_count << std::endl;

    return true;

}



#endif // MODEL_OPTIMIZATION_ROUTINES_H