#ifndef SMOOTH_H
#define SMOOTH_H

#include <iostream>
#include "../../include/config.h"
#include "../../include/utils.h"


void CG_smooth(CUSTOMREAL* arr_in, CUSTOMREAL* arr_out, CUSTOMREAL lr, CUSTOMREAL lt, CUSTOMREAL lp) {
    // arr: array to be smoothed
    // lr: smooth length on r
    // lt: smooth length on theta
    // lp: smooth length on phi

    // arrays :
    // x_array  model
    // g_array  gradient
    // d_array  descent direction
    // Ap  stiffness scalar (laplacian)
    // pAp = d_array*Ap
    // rr = g_array * g_array (dot product)

    const int max_iter_cg = 300;
    const CUSTOMREAL xtol = 0.001;
    const bool use_scaling = true;

    // allocate memory
    CUSTOMREAL* x_array = new CUSTOMREAL[loc_I*loc_J*loc_K];
    CUSTOMREAL* g_array = new CUSTOMREAL[loc_I*loc_J*loc_K];
    CUSTOMREAL* d_array = new CUSTOMREAL[loc_I*loc_J*loc_K];
    CUSTOMREAL* Ap = new CUSTOMREAL[loc_I*loc_J*loc_K];
    CUSTOMREAL pAp=_0_CR, rr_0=_0_CR, rr=_0_CR, rr_new=_0_CR, aa=_0_CR, bb=_0_CR, tmp=_0_CR;

    CUSTOMREAL scaling_A=_1_CR, scaling_coeff = _1_CR;

    if (use_scaling) {
        // calculate scaling factor
        scaling_A = std::sqrt(_1_CR / (_8_CR * PI * lr * lt * lp));
        // scaling coefficient for gradient
        scaling_coeff = find_absmax(arr_in, loc_I*loc_J*loc_K);
        //tmp = scaling_coeff;
        //allreduce_cr_single_max(tmp, scaling_coeff);
        //if (scaling_coeff == _0_CR)
        if (isZero(scaling_coeff))
            scaling_coeff = _1_CR;
    }
    // std out scaling factors
    //if (myrank == 0) {
        std::cout << "scaling_A = " << scaling_A << std::endl;
        std::cout << "scaling_coeff = " << scaling_coeff << std::endl;
    //}


    // array initialization
    for (int i=0; i<loc_I*loc_J*loc_K; i++) {
        x_array[i] = _0_CR;
        g_array[i] = arr_in[i]/scaling_coeff;
        d_array[i] = arr_in[i]/scaling_coeff;
        Ap[i] = _0_CR;
    }

    // initial rr
    rr = dot_product(g_array, g_array, loc_I*loc_J*loc_K);
    //tmp = rr;
    // sum all rr among all processors
    //allreduce_cr_single(tmp, rr);
    rr_0 = rr; // record initial rr

    // CG loop
    for (int iter=0; iter<max_iter_cg; iter++) {

        // calculate laplacian
        for (int k = 1; k < loc_K-1; k++) {
            for (int j = 1; j < loc_J-1; j++) {
                for (int i = 1; i < loc_I-1; i++) {
                    //scaling_coeff = std::max(scaling_coeff, std::abs(arr[I2V(i,j,k)]));

                    Ap[I2V(i,j,k)] = (
                        - _2_CR * (lr+lt+lp) * d_array[I2V(i,j,k)] \
                        + lr * (d_array[I2V(i,j,k+1)] + d_array[I2V(i,j,k-1)]) \
                        + lt * (d_array[I2V(i,j+1,k)] + d_array[I2V(i,j-1,k)]) \
                        + lp * (d_array[I2V(i+1,j,k)] + d_array[I2V(i-1,j,k)]) \
                    );

                    // scaling
                    Ap[I2V(i,j,k)] = Ap[I2V(i,j,k)]*scaling_A;
                    //Ap[I2V(i,j,k)] = (d_array[I2V(i,j,k)]+Ap[I2V(i,j,k)])*scaling_A;
                }
            }
        }

        // get the values on the boundaries
        //grid.send_recev_boundary_data(Ap);

        // calculate pAp
        pAp = dot_product(d_array, Ap, loc_I*loc_J*loc_K);
        //tmp = pAp;
        // sum all pAp among all processors
        //allreduce_cr_single(tmp, pAp);

        // compute step length
        aa = rr / pAp;

        // update x_array (model)
        for (int i=0; i<loc_I*loc_J*loc_K; i++) {
            x_array[i] += aa * d_array[i];
        }

        // update g_array (gradient)
        for (int i=0; i<loc_I*loc_J*loc_K; i++) {
            g_array[i] -= aa * Ap[i];
        }

        // update rr
        rr_new = dot_product(g_array, g_array, loc_I*loc_J*loc_K);
        //tmp = rr_new;
        // sum all rr among all processors
        //allreduce_cr_single(tmp, rr_new);

        // update d_array (descent direction)
        bb = rr_new / rr;
        for (int i=0; i<loc_I*loc_J*loc_K; i++) {
            d_array[i] = g_array[i] + bb * d_array[i];
        }

        //if (myrank == 0 && iter%100==0){//} && if_verbose) {
        if (iter%100==0){//} && if_verbose) {
            std::cout << "iter: " << iter << " rr: " << rr << " rr_new: " << rr_new << " rr/rr_0: " << rr/rr_0 << " pAp: " << pAp << " aa: " << aa << " bb: " << bb << std::endl;
        }

        // update rr
        rr = rr_new;

        // check convergence
        if (rr / rr_0 < xtol) {
            std::cout << "CG converged in " << iter << " iterations." << std::endl;
            break;
        }

    } // end of CG loop


    // copy x_array to arr_out
    for (int i=0; i<loc_I*loc_J*loc_K; i++) {
        arr_out[i] = x_array[i]*scaling_coeff;
    }

    // deallocate
    delete[] x_array;
    delete[] g_array;
    delete[] d_array;
    delete[] Ap;

}



#endif