#ifndef LBFGS_H
#define LBFGS_H

#include "config.h"
#include "utils.h"
#include "grid.h"
#include "mpi_funcs.h"

void store_model_and_gradient(Grid& grid, int i_inv) {

    if (subdom_main) {

        if (i_inv < Mbfgs){
            for (int k = 0; k < loc_K; k++) {
                for (int j = 0; j < loc_J; j++) {
                    for (int i = 0; i < loc_I; i++) {
                        grid.Ks_model_store_loc[  I2VLBFGS(i_inv,i,j,k)] = grid.fun_loc[I2V(i,j,k)];
                        grid.Ks_grad_store_loc[   I2VLBFGS(i_inv,i,j,k)] = grid.Ks_update_loc[  I2V(i,j,k)];
                    }
                }
            }

        } else {
            // replace the current stored models and gradients
            for (int istore = 0; istore < Mbfgs-1; istore++){
                for (int k = 0; k < loc_K; k++) {
                    for (int j = 0; j < loc_J; j++) {
                        for (int i = 0; i < loc_I; i++) {
                            grid.Ks_model_store_loc[  I2VLBFGS(istore,i,j,k)] = grid.Ks_model_store_loc[  I2VLBFGS(istore+1,i,j,k)];
                            grid.Ks_grad_store_loc[   I2VLBFGS(istore,i,j,k)] = grid.Ks_grad_store_loc[   I2VLBFGS(istore+1,i,j,k)];
                        }
                    }
                }
            }

            // store new model and gradient
            for (int k = 0; k < loc_K; k++) {
                for (int j = 0; j < loc_J; j++) {
                    for (int i = 0; i < loc_I; i++) {
                        grid.Ks_model_store_loc[  I2VLBFGS(Mbfgs-1,i,j,k)] = grid.fun_loc[I2V(i,j,k)];
                        grid.Ks_grad_store_loc[   I2VLBFGS(Mbfgs-1,i,j,k)] = grid.Ks_update_loc[  I2V(i,j,k)];
                    }
                }
            }
        }

   } // end subdom main

}


void calculate_descent_direction_lbfgs(Grid& grid, int i_inv) {

    if (subdom_main) {
        int imin = 0;
        int imax = 0;
        if (i_inv >= Mbfgs)
            imax = Mbfgs-2;
        else
            imax = i_inv-1;

        CUSTOMREAL* ak_store = new CUSTOMREAL[imax-imin+1];
        CUSTOMREAL* pk_store = new CUSTOMREAL[imax-imin+1];

        CUSTOMREAL* desc_wks_Ks   = new CUSTOMREAL[loc_I*loc_J*loc_K];

        CUSTOMREAL* wks_1_Ks   = new CUSTOMREAL[loc_I*loc_J*loc_K];
        CUSTOMREAL* wks_2_Ks   = new CUSTOMREAL[loc_I*loc_J*loc_K];

        std::cout << "imax, imin: " << imax << ", " << imin << std::endl;

        // initialize
        for (int i = 0; i < imax-imin+1; i++) {
            ak_store[i] = 0.0; // alpha
            pk_store[i] = 0.0; //rho
        }

        for (int k = 0; k < loc_K; k++) {
            for (int j = 0; j < loc_J; j++) {
                for (int i = 0; i < loc_I; i++) {
                    desc_wks_Ks[I2V(i,j,k)]   = grid.Ks_grad_store_loc[  I2VLBFGS(imax+1,i,j,k)]; // d
                    wks_1_Ks[I2V(i,j,k)]   = grid.Ks_grad_store_loc[  I2VLBFGS(imax+1,i,j,k)] - grid.Ks_grad_store_loc[  I2VLBFGS(imax,i,j,k)];
                }
            }
        }

        //calculate l2 norms
        CUSTOMREAL norm_yiter = 0.0;

        norm_yiter += calc_l2norm(wks_1_Ks,loc_I*loc_J*loc_K);

        CUSTOMREAL tmp = norm_yiter;
        allreduce_cr_single(tmp, norm_yiter);

        for (int iinv = imax; iinv >= imin; iinv--){

            for (int k = 0; k < loc_K; k++) {
                for (int j = 0; j < loc_J; j++) {
                    for (int i = 0; i < loc_I; i++) {
                        wks_1_Ks[I2V(i,j,k)]   = grid.Ks_grad_store_loc[  I2VLBFGS(iinv+1,i,j,k)]   - grid.Ks_grad_store_loc[I2VLBFGS(iinv,i,j,k)];
                        wks_2_Ks[I2V(i,j,k)]   = grid.Ks_model_store_loc[  I2VLBFGS(iinv+1,i,j,k)]   - grid.Ks_model_store_loc[I2VLBFGS(iinv,i,j,k)];
                    }
                }
            }

            CUSTOMREAL pk = _0_CR;
            pk += dot_product(wks_1_Ks,wks_2_Ks,loc_I*loc_J*loc_K);
            CUSTOMREAL tmp = pk;
            allreduce_cr_single(tmp, pk);
            pk_store[iinv] = _1_CR / pk;

            CUSTOMREAL ak = _0_CR;
            ak += dot_product(wks_2_Ks, desc_wks_Ks,loc_I*loc_J*loc_K);
            // print ak
            if(myrank== 0) std::cout << "ak: " << ak << std::endl;
            tmp = ak;
            allreduce_cr_single(tmp, ak);
            // print ak
            if(myrank== 0) std::cout << "ak gathered: " << ak << std::endl;
            ak_store[iinv] = pk_store[iinv] * ak;

            for (int k = 0; k < loc_K; k++) {
                for (int j = 0; j < loc_J; j++) {
                    for (int i = 0; i < loc_I; i++) {
                        desc_wks_Ks[I2V(i,j,k)]   -= ak_store[iinv] * wks_1_Ks[I2V(i,j,k)];
                    }
                }
            }

        } // end loop iinv

        // Nocedal's default preconditionning
        CUSTOMREAL pk = _1_CR / (pk_store[imax] * norm_yiter);
        for (int k = 0; k < loc_K; k++) {
            for (int j = 0; j < loc_J; j++) {
                for (int i = 0; i < loc_I; i++) {
                    desc_wks_Ks[I2V(i,j,k)]   *= pk;
                }
            }
        }

        // print ak_store and pk_store
        if (myrank == 0) {
            std::cout << "ak_store: ";
            for (int i = 0; i < imax+1; i++) {
                std::cout << ak_store[i] << " ";
            }
            std::cout << std::endl;
            std::cout << "pk_store: ";
            for (int i = 0; i < imax+1; i++) {
                std::cout << pk_store[i] << " ";
            }
            std::cout << std::endl;
        }


        // custom preconditionning
        // diagonal preconditionner

        for (int iinv = imin; iinv <= imax; iinv++){

            for (int k = 0; k < loc_K; k++) {
                for (int j = 0; j < loc_J; j++) {
                    for (int i = 0; i < loc_I; i++) {
                        wks_1_Ks[I2V(i,j,k)]   = grid.Ks_grad_store_loc[  I2VLBFGS(iinv+1,i,j,k)]   - grid.Ks_grad_store_loc[I2VLBFGS(iinv,i,j,k)];
                        wks_2_Ks[I2V(i,j,k)]   = grid.Ks_model_store_loc[  I2VLBFGS(iinv+1,i,j,k)]   - grid.Ks_model_store_loc[I2VLBFGS(iinv,i,j,k)];
                    }
                }
            }

            CUSTOMREAL beta = _0_CR;
            beta += dot_product(wks_1_Ks,desc_wks_Ks,loc_I*loc_J*loc_K);
            CUSTOMREAL tmp = beta;
            allreduce_cr_single(tmp, beta);
            beta = pk_store[iinv] * beta;

            for (int k = 0; k < loc_K; k++) {
                for (int j = 0; j < loc_J; j++) {
                    for (int i = 0; i < loc_I; i++) {
                        desc_wks_Ks[I2V(i,j,k)]   += (ak_store[iinv]-beta) * wks_2_Ks[I2V(i,j,k)];
                    }
                }
            }
        }

        // descent directions
        for (int k = 0; k < loc_K; k++) {
            for (int j = 0; j < loc_J; j++) {
                for (int i = 0; i < loc_I; i++) {
                    grid.Ks_descent_dir_loc[I2V(i,j,k)]   = desc_wks_Ks[I2V(i,j,k)];
                }
            }
        }

        delete [] ak_store;
        delete [] pk_store;
        delete [] desc_wks_Ks;
        delete [] wks_1_Ks;
        delete [] wks_2_Ks;

    } // end of subdom_main

}


inline void calc_laplacian_field(Grid& grid, CUSTOMREAL* arr_in, CUSTOMREAL* arr_res) {
    if (subdom_main) {

        CUSTOMREAL lr = 1.0;
        CUSTOMREAL lt = 1.0;
        CUSTOMREAL lp = 1.0;

        // calculate L(m)
        for (int k = 1; k < loc_K-1; k++) {
            for (int j = 1; j < loc_J-1; j++) {
                for (int i = 1; i < loc_I-1; i++) {
                    arr_res[I2V(i,j,k)] = \
                     - _2_CR * (lr+lt+lp) * arr_in[I2V(i,j,k)] \
                        + lr * (arr_in[I2V(i,j,k+1)] + arr_in[I2V(i,j,k-1)]) \
                        + lt * (arr_in[I2V(i,j+1,k)] + arr_in[I2V(i,j-1,k)]) \
                        + lp * (arr_in[I2V(i+1,j,k)] + arr_in[I2V(i-1,j,k)]);
                }
            }
        }

        // communicate with adjacent subdomains
        grid.send_recev_boundary_data(arr_res);
    }
}

// add 1/2 * L(m)^2 to the objective function
inline CUSTOMREAL add_regularization_obj(Grid& grid) {

    CUSTOMREAL regularization_term = _0_CR;

    if (subdom_main) {
        CUSTOMREAL* tmp_arr = new CUSTOMREAL[loc_I*loc_J*loc_K];
        int ngrid = loc_I*loc_J*loc_K;
        for (int i = 0; i < ngrid; i++) tmp_arr[i] = _0_CR;

        regularization_term += calc_l2norm(grid.fun_regularization_penalty_loc, ngrid);

        // gather from all subdomain
        CUSTOMREAL tmp = regularization_term;
        allreduce_cr_single(tmp, regularization_term);

        delete [] tmp_arr;
    }

    synchronize_all_world();

    return _0_5_CR * regularization_term * regularization_weight;

}


// add grad(L(m)) to gradient
inline void add_regularization_grad(Grid& grid) {
    if (subdom_main) {
        CUSTOMREAL* tmp_arr = new CUSTOMREAL[loc_I*loc_J*loc_K];
        int ngrid = loc_I*loc_J*loc_K;

        // calculate LL(m) on fun (Ks)
        for (int i = 0; i < ngrid; i++) tmp_arr[i] = _0_CR;
        calc_laplacian_field(grid, grid.fun_loc, grid.fun_regularization_penalty_loc);
        calc_laplacian_field(grid, grid.fun_regularization_penalty_loc, grid.Ks_regularization_penalty_loc);

        // add LL(m) to gradient
        for (int k = 0; k < loc_K; k++) {
            for (int j = 0; j < loc_J; j++) {
                for (int i = 0; i < loc_I; i++) {

                    //grid.Ks_regularization_penalty_loc[I2V(i,j,k)]   /= Linf_Ks;
                    //grid.Keta_regularization_penalty_loc[I2V(i,j,k)] /= Linf_Keta;
                    //grid.Kxi_regularization_penalty_loc[I2V(i,j,k)]  /= Linf_Kxi;

                    grid.Ks_update_loc[I2V(i,j,k)]   += regularization_weight * grid.Ks_regularization_penalty_loc[I2V(i,j,k)];
                }
            }
        }

        delete [] tmp_arr;

    }

}


// compute initial guess for step length to try for line search
void initial_guess_step(Grid& grid, CUSTOMREAL& step_size, CUSTOMREAL SC_VAL) {
    // find the max value in descent_direction
    CUSTOMREAL max_val = _0_CR;
    for (int k = 0; k < loc_K; k++) {
        for (int j = 0; j < loc_J; j++) {
            for (int i = 0; i < loc_I; i++) {
                if (std::abs(grid.Ks_descent_dir_loc[I2V(i,j,k)]  ) > max_val) max_val = std::abs(grid.Ks_descent_dir_loc[I2V(i,j,k)]);
            }
        }
    }

    // reduce
    allreduce_cr_single_max(max_val, max_val);

    // debug print max_val
    if (subdom_main) {
        std::cout << "max_val: " << max_val << std::endl;
    }

    // step size
    step_size = SC_VAL / max_val;

}


#endif