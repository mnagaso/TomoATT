#ifndef MODEL_UPDATE_H
#define MODEL_UPDATE_H

#include <cmath>
#include "config.h"
#include "grid.h"
#include "smooth.h"
#include "smooth_descent_dir.h"
#include "lbfgs.h"


// K*_loc -> K*_update_loc
void smooth_kernels(Grid& grid) {
    if (subdom_main){
        // initiaize update params
        for (int k = 0; k < loc_K; k++) {
            for (int j = 0; j < loc_J; j++) {
                for (int i = 0; i < loc_I; i++) {
                    grid.Ks_update_loc[I2V(i,j,k)]   = _0_CR;
                    grid.Keta_update_loc[I2V(i,j,k)] = _0_CR;
                    grid.Kxi_update_loc[I2V(i,j,k)]  = _0_CR;
                }
            }
        }

        if (id_sim == 0) { // calculation of the update model is only done in the main simultaneous run

            if (smooth_method == 0) {
                // grid based smoothing
                smooth_inv_kernels_orig(grid);
            } else if (smooth_method == 1) {
                // CG smoothing
                smooth_inv_kernels_CG(grid, smooth_lr, smooth_lt, smooth_lp);
            }

            // shared values on the boundary
            grid.send_recev_boundary_data(grid.Ks_update_loc);
            grid.send_recev_boundary_data(grid.Keta_update_loc);
            grid.send_recev_boundary_data(grid.Kxi_update_loc);

        } // end if id_sim == 0

        // send the updated model to all the simultaneous run
        broadcast_cr_inter_sim(grid.Ks_update_loc, loc_I*loc_J*loc_K, 0);
        broadcast_cr_inter_sim(grid.Kxi_update_loc, loc_I*loc_J*loc_K, 0);
        broadcast_cr_inter_sim(grid.Keta_update_loc, loc_I*loc_J*loc_K, 0);

    }
}


void smooth_descent_direction(Grid& grid){
    if (subdom_main){
       if (id_sim == 0) { // calculation of the update model is only done in the main simultaneous run
            if (smooth_method == 0) {
                // grid based smoothing
                smooth_descent_dir(grid);
            } else if (smooth_method == 1) {
                // CG smoothing not implemented yet
                smooth_descent_dir(grid);
            }

            // shared values on the boundary
            grid.send_recev_boundary_data(  grid.Ks_descent_dir_loc);
            grid.send_recev_boundary_data(grid.Keta_descent_dir_loc);
            grid.send_recev_boundary_data( grid.Kxi_descent_dir_loc);

        } // end if id_sim == 0

        // send the updated model to all the simultaneous run
        broadcast_cr_inter_sim(  grid.Ks_descent_dir_loc, loc_I*loc_J*loc_K, 0);
        broadcast_cr_inter_sim( grid.Kxi_descent_dir_loc, loc_I*loc_J*loc_K, 0);
        broadcast_cr_inter_sim(grid.Keta_descent_dir_loc, loc_I*loc_J*loc_K, 0);

    }
}


void calc_descent_direction(Grid& grid, int i_inv) {

    if (subdom_main) {

        // routines for LBFGS
        if (optim_method == LBFGS_MODE) {

            // calculate the descent direction
            if (i_inv > 0) {
                calculate_descent_direction_lbfgs(grid, i_inv);
            // use gradient for the first iteration
            } else {
                // sum kernels among all simultaneous runs
                sumup_kernels(grid);
                // smooth kernels
                smooth_kernels(grid);

                int n_grid = loc_I*loc_J*loc_K;
                std::memcpy(grid.Ks_descent_dir_loc,   grid.Ks_update_loc,   n_grid*sizeof(CUSTOMREAL));
                std::memcpy(grid.Keta_descent_dir_loc, grid.Keta_update_loc, n_grid*sizeof(CUSTOMREAL));
                std::memcpy(grid.Kxi_descent_dir_loc,  grid.Kxi_update_loc,  n_grid*sizeof(CUSTOMREAL));


                //// inverse the gradient to fit the update scheme for LBFGS
                for (int i = 0; i < n_grid; i++){
                    grid.Ks_descent_dir_loc[i]   = -grid.Ks_descent_dir_loc[i];
                    grid.Keta_descent_dir_loc[i] = -grid.Keta_descent_dir_loc[i];
                    grid.Kxi_descent_dir_loc[i]  = -grid.Kxi_descent_dir_loc[i];
                }
            }

        }


    } // end if subdom_main

    synchronize_all_world();

}


void set_new_model(Grid& grid, CUSTOMREAL step_size_new) {

    if (subdom_main) {

        // for LBFGS mode. K*_update_loc is not directly the descent direction but smoothed gradient
        // so for non LBFGS_mode, nextstep will be calculated with K*_update_loc
        if (optim_method != LBFGS_MODE) {
            // get the scaling factor
            CUSTOMREAL Linf_Ks = _0_CR, Linf_Keta = _0_CR, Linf_Kxi = _0_CR;
            for (int k = 1; k < loc_K-1; k++) {
                for (int j = 1; j < loc_J-1; j++) {
                    for (int i = 1; i < loc_I-1; i++) {
                        Linf_Ks   = std::max(Linf_Ks,   std::abs(grid.Ks_update_loc[I2V(i,j,k)]));
                        Linf_Keta = std::max(Linf_Keta, std::abs(grid.Keta_update_loc[I2V(i,j,k)]));
                        Linf_Kxi  = std::max(Linf_Kxi,  std::abs(grid.Kxi_update_loc[I2V(i,j,k)]));
                    }
                }
            }

            // get the maximum scaling factor among subdomains
            CUSTOMREAL Linf_tmp;
            allreduce_cr_single_max(Linf_Ks, Linf_tmp);   Linf_Ks = Linf_tmp;
            allreduce_cr_single_max(Linf_Keta, Linf_tmp); Linf_Keta = Linf_tmp;
            allreduce_cr_single_max(Linf_Kxi, Linf_tmp);  Linf_Kxi = Linf_tmp;

            CUSTOMREAL Linf_all = _0_CR;
            Linf_all = std::max(Linf_Ks, std::max(Linf_Keta, Linf_Kxi));
            Linf_Ks = Linf_all;
            Linf_Keta = Linf_all;
            Linf_Kxi = Linf_all;

            if (myrank == 0 && id_sim == 0)
                //std::cout << "Scaring factor for all kernels: " << Linf_all << std::endl;
                std::cout << "Scaring factor for model update for Ks, Keta, Kx, stepsize: " << Linf_Ks << ", " << Linf_Keta << ", " << Linf_Kxi << ", " << step_size_new << std::endl;

            // update the model
            for (int k = 0; k < loc_K; k++) {
                for (int j = 0; j < loc_J; j++) {
                    for (int i = 0; i < loc_I; i++) {

                        // update
                        grid.fun_loc[I2V(i,j,k)] *= (_1_CR - grid.Ks_update_loc[I2V(i,j,k)  ] / (Linf_Ks   /step_size_new) );
                        grid.xi_loc[I2V(i,j,k)]  -=          grid.Kxi_update_loc[I2V(i,j,k) ] / (Linf_Kxi  /step_size_new)  ;
                        grid.eta_loc[I2V(i,j,k)] -=          grid.Keta_update_loc[I2V(i,j,k)] / (Linf_Keta /step_size_new)  ;

                        grid.fac_b_loc[I2V(i,j,k)] = _1_CR - _2_CR * grid.xi_loc[I2V(i,j,k)];
                        grid.fac_c_loc[I2V(i,j,k)] = _1_CR + _2_CR * grid.xi_loc[I2V(i,j,k)];
                        grid.fac_f_loc[I2V(i,j,k)] =       - _2_CR * grid.eta_loc[I2V(i,j,k)];
                    }
                }
            }


        } else { // for LBFGS routine


//            // get the scaling factor
//            CUSTOMREAL Linf_Ks = _0_CR, Linf_Keta = _0_CR, Linf_Kxi = _0_CR;
//            CUSTOMREAL Linf_all = _0_CR;
//            for (int k = 0; k < loc_K; k++) {
//                for (int j = 0; j < loc_J; j++) {
//                    for (int i = 0; i < loc_I; i++) {
//                        Linf_Ks   = std::max(Linf_Ks,   std::abs(grid.Ks_descent_dir_loc[I2V(i,j,k)]));
//                        Linf_Keta = std::max(Linf_Keta, std::abs(grid.Keta_descent_dir_loc[I2V(i,j,k)]));
//                        Linf_Kxi  = std::max(Linf_Kxi,  std::abs(grid.Kxi_descent_dir_loc[I2V(i,j,k)]));
//                    }
//                }
//            }
//
//            // get the maximum scaling factor among subdomains
//            CUSTOMREAL Linf_tmp;
//            allreduce_cr_single_max(Linf_Ks, Linf_tmp);   Linf_Ks   = Linf_tmp;
//            allreduce_cr_single_max(Linf_Keta, Linf_tmp); Linf_Keta = Linf_tmp;
//            allreduce_cr_single_max(Linf_Kxi, Linf_tmp);  Linf_Kxi  = Linf_tmp;
//
//            Linf_all = std::max(Linf_Ks, std::max(Linf_Keta, Linf_Kxi));
//            Linf_Ks = Linf_all;
//            Linf_Keta = Linf_all;
//            Linf_Kxi = Linf_all;
//
//
//            if (myrank == 0 && id_sim == 0)
//                //std::cout << "Scaring factor for all kernels: " << Linf_all << std::endl;
//                std::cout << "Scaring factor for model update for Ks, Keta, Kx, stepsize: " << Linf_Ks << ", " << Linf_Keta << ", " << Linf_Kxi << ", " << step_size_new << std::endl;

            // update the model
            for (int k = 0; k < loc_K; k++) {
                for (int j = 0; j < loc_J; j++) {
                    for (int i = 0; i < loc_I; i++) {
                        // update
                        //grid.fun_loc[I2V(i,j,k)] *= (_1_CR -  grid.Ks_descent_dir_loc[I2V(i,j,k)]   * step_size_new);
                        //grid.xi_loc[I2V(i,j,k)]  -=           grid.Kxi_descent_dir_loc[I2V(i,j,k) ] * step_size_new;
                        //grid.eta_loc[I2V(i,j,k)] -=           grid.Keta_descent_dir_loc[I2V(i,j,k)] * step_size_new;
                        grid.fun_loc[I2V(i,j,k)] += grid.Ks_descent_dir_loc[I2V(i,j,k)]  *step_size_new;
                        grid.xi_loc[I2V(i,j,k)]  += grid.Kxi_descent_dir_loc[I2V(i,j,k)] *step_size_new;
                        grid.eta_loc[I2V(i,j,k)] += grid.Keta_descent_dir_loc[I2V(i,j,k)]*step_size_new;

                        grid.fac_b_loc[I2V(i,j,k)] = _1_CR - _2_CR * grid.xi_loc[I2V(i,j,k)];
                        grid.fac_c_loc[I2V(i,j,k)] = _1_CR + _2_CR * grid.xi_loc[I2V(i,j,k)];
                        grid.fac_f_loc[I2V(i,j,k)] =       - _2_CR * grid.eta_loc[I2V(i,j,k)];
                    }
                }
            }

        } // end LBFGS model update

        // shared values on the boundary
        grid.send_recev_boundary_data(grid.fun_loc);
        grid.send_recev_boundary_data(grid.xi_loc);
        grid.send_recev_boundary_data(grid.eta_loc);
        grid.send_recev_boundary_data(grid.fac_b_loc);
        grid.send_recev_boundary_data(grid.fac_c_loc);
        grid.send_recev_boundary_data(grid.fac_f_loc);

    }
}



// compute p_k * grad(f_k)
CUSTOMREAL compute_q_k(Grid& grid) {

    CUSTOMREAL tmp_qk = _0_CR;

    if (subdom_main) {

        // grad * descent_direction
        for (int k = 0; k < loc_K; k++) {
            for (int j = 0; j < loc_J; j++) {
                for (int i = 0; i < loc_I; i++) {
                    // we can add here a precondition (here use *_update_loc for storing descent direction)
                    tmp_qk += grid.Ks_update_loc[I2V(i,j,k)]   * grid.Ks_descent_dir_loc[I2V(i,j,k)]
                            + grid.Keta_update_loc[I2V(i,j,k)] * grid.Keta_descent_dir_loc[I2V(i,j,k)]
                            + grid.Kxi_update_loc[I2V(i,j,k)]  * grid.Kxi_descent_dir_loc[I2V(i,j,k)];
                    //tmp_qk += grid.Ks_loc[I2V(i,j,k)]
                    //        + grid.Keta_loc[I2V(i,j,k)]
                    //        + grid.Kxi_loc[I2V(i,j,k)];

                }
            }
        }

        // add tmp_qk among all subdomain
        allreduce_cr_single(tmp_qk, tmp_qk);
    }

    return tmp_qk;

}


// check if the wolfe conditions are satisfied
bool check_wolfe_cond(Grid& grid, \
                    CUSTOMREAL f_k, CUSTOMREAL f_new, \
                    CUSTOMREAL q_k, CUSTOMREAL q_new, \
                    CUSTOMREAL& td, CUSTOMREAL& tg, CUSTOMREAL& step_size_sub) {

    bool step_accepted = false;

    // check  if the wolfe conditions are satisfied and update the step_size_sub
    CUSTOMREAL qpt = (f_new - f_k) / step_size_sub;

    if (subdom_main) {
        // good step size
        if (qpt   <= wolfe_c1*q_k \
         && q_new >= wolfe_c2*q_k ){
            if (myrank==0)
                std::cout << "Wolfe rules:  step accepted" << std::endl;
            step_accepted = true;
        } else {
            // modify the stepsize

            if (wolfe_c1*q_k < qpt) {
                td = step_size_sub;
                if (myrank==0)
                    std::cout << "Wolfe rules:  right step size updated." << std::endl;
            }
            //if (q_new < wolfe_c2*q_k) {
            if (qpt <= wolfe_c1*q_k && q_new < wolfe_c2*q_k) {
                tg = step_size_sub;
                if (myrank==0)
                    std::cout << "Wolfe rules:  left step size updated." << std::endl;
            }
            if (isZero(td)) {
                step_size_sub = _2_CR * step_size_sub;
                if (myrank==0)
                    std::cout << "Wolfe rules:  step size too small. Increased to " << step_size_sub << std::endl;
            } else {
                step_size_sub = _0_5_CR * (td+tg);
                if (myrank==0)
                    std::cout << "Wolfe rules:  step size too large. Decreased to " << step_size_sub << std::endl;
            }
        }

        // share step_accepted among all subdomains
        broadcast_bool_single(step_accepted,0);
        broadcast_cr_single(step_size_sub,0);
    } // end if subdom_main

    broadcast_bool_single_sub(step_accepted,0);
    broadcast_cr_single_sub(step_size_sub,0);

    return step_accepted;
}



#endif // MODEL_UPDATE_H