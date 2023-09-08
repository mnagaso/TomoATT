#ifndef SMOOTH_H
#define SMOOTH_H

#include <iostream>
#include "grid.h"
#include "config.h"


void CG_smooth(Grid& grid, CUSTOMREAL* arr_in, CUSTOMREAL* arr_out, CUSTOMREAL lr, CUSTOMREAL lt, CUSTOMREAL lp) {
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

    const int max_iter_cg = 1000;
    const CUSTOMREAL xtol = 0.01;
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
        tmp = scaling_coeff;
        allreduce_cr_single_max(tmp, scaling_coeff);
        if (scaling_coeff == _0_CR)
            scaling_coeff = _1_CR;
    }
    // std out scaling factors
    if (myrank == 0) {
        std::cout << "scaling_A = " << scaling_A << std::endl;
        std::cout << "scaling_coeff = " << scaling_coeff << std::endl;
    }


    // array initialization
    for (int i=0; i<loc_I*loc_J*loc_K; i++) {
        x_array[i] = _0_CR;
        g_array[i] = arr_in[i]/scaling_coeff;
        d_array[i] = arr_in[i]/scaling_coeff;
        Ap[i] = _0_CR;
    }

    // initial rr
    rr = dot_product(g_array, g_array, loc_I*loc_J*loc_K);
    tmp = rr;
    // sum all rr among all processors
    allreduce_cr_single(tmp, rr);
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
                    Ap[I2V(i,j,k)] = -Ap[I2V(i,j,k)]*scaling_A;

                }
            }
        }

        // get the values on the boundaries
        grid.send_recev_boundary_data(Ap);

        // calculate pAp
        pAp = dot_product(d_array, Ap, loc_I*loc_J*loc_K);
        tmp = pAp;
        // sum all pAp among all processors
        allreduce_cr_single(tmp, pAp);

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
        tmp = rr_new;
        // sum all rr among all processors
        allreduce_cr_single(tmp, rr_new);

        // update d_array (descent direction)
        bb = rr_new / rr;
        for (int i=0; i<loc_I*loc_J*loc_K; i++) {
            d_array[i] = g_array[i] + bb * d_array[i];
        }

        if (myrank == 0 && iter%100==0){//} && if_verbose) {
            std::cout << "iter: " << iter << " rr: " << rr << " rr_new: " << rr_new << " rr/rr_0: " << rr/rr_0 << " pAp: " << pAp << " aa: " << aa << " bb: " << bb << std::endl;
        }

        // update rr
        rr = rr_new;

        // check convergence
        if (rr / rr_0 < xtol) {
            if (myrank == 0) {
                std::cout << "CG converged in " << iter << " iterations." << std::endl;
            }
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


void smooth_inv_kernels_CG(Grid& grid, CUSTOMREAL lr, CUSTOMREAL lt, CUSTOMREAL lp) {
    // smoothing kernels with Conjugate Gradient method
    // lr,lt,lp: smoothing length in r, theta, phi direction

    // smooth Ks
    CG_smooth(grid, grid.Ks_loc, grid.Ks_update_loc, lr, lt, lp);
    // smooth Keta
    CG_smooth(grid, grid.Keta_loc, grid.Keta_update_loc, lr, lt, lp);
    // smooth Kxi
    CG_smooth(grid, grid.Kxi_loc, grid.Kxi_update_loc, lr, lt, lp);

}





// original method for smoothing kernels
inline void smooth_inv_kernels_orig(Grid& grid, InputParams& IP){
    // necessary params
    CUSTOMREAL r_r, r_t, r_p;
    CUSTOMREAL r_r_ani, r_t_ani, r_p_ani;
    
    int kdr = 0, jdt = 0, idp = 0;
    int kdr_ani = 0, jdt_ani = 0, idp_ani = 0;

    int k_start = 0; //grid.get_k_start_loc();
    int j_start = 0; //grid.get_j_start_loc();
    int i_start = 0; //grid.get_i_start_loc();
    int k_end   = ngrid_k; //grid.get_k_end_loc();
    int j_end   = ngrid_j; //grid.get_j_end_loc();
    int i_end   = ngrid_i; //grid.get_i_end_loc();

    CUSTOMREAL weight   = _1_CR;
    CUSTOMREAL * taper  = IP.get_depth_taper();

    //
    // sum the kernel values on the inversion grid
    //

    for (int i_grid = 0; i_grid < n_inv_grids; i_grid++) {

        // init coarser kernel arrays
            // velocity
        for (int k = 0; k < n_inv_K_loc; k++) {
            for (int j = 0; j < n_inv_J_loc; j++) {
                for (int i = 0; i < n_inv_I_loc; i++) {
                    grid.Ks_inv_loc[  I2V_INV_KNL(i,j,k)] = _0_CR;
                }
            }
        }
            // anisotropy
        for (int k = 0; k < n_inv_K_loc_ani; k++) {
            for (int j = 0; j < n_inv_J_loc_ani; j++) {
                for (int i = 0; i < n_inv_I_loc_ani; i++) {
                    grid.Keta_inv_loc[I2V_INV_ANI_KNL(i,j,k)] = _0_CR;
                    grid.Kxi_inv_loc[ I2V_INV_ANI_KNL(i,j,k)] = _0_CR;
                }
            }
        }


        for (int k = k_start; k < k_end; k++) {
            CUSTOMREAL r_glob = grid.get_r_min() + k*grid.get_delta_r(); // global coordinate of r
            r_r = -_1_CR;
            for (int ii_invr = 0; ii_invr < n_inv_K_loc-1; ii_invr++){
                // increasing or decreasing order
                if ((r_glob - grid.r_loc_inv[I2V_INV_GRIDS_1DK(ii_invr,i_grid)]) * (r_glob - grid.r_loc_inv[I2V_INV_GRIDS_1DK(ii_invr+1,i_grid)]) <= 0) {
                    kdr = ii_invr;
                    r_r = (r_glob - grid.r_loc_inv[I2V_INV_GRIDS_1DK(ii_invr,i_grid)]) / (grid.r_loc_inv[I2V_INV_GRIDS_1DK(ii_invr+1,i_grid)] - grid.r_loc_inv[I2V_INV_GRIDS_1DK(ii_invr,i_grid)]);
                    break;
                }
            }
            r_r_ani = -_1_CR;
            for (int ii_invr = 0; ii_invr < n_inv_K_loc_ani-1; ii_invr++){
                // increasing or decreasing order
                if ((r_glob - grid.r_loc_inv_ani[I2V_INV_ANI_GRIDS_1DK(ii_invr,i_grid)]) * (r_glob - grid.r_loc_inv_ani[I2V_INV_ANI_GRIDS_1DK(ii_invr+1,i_grid)]) <= 0) {
                    kdr_ani = ii_invr;
                    r_r_ani = (r_glob - grid.r_loc_inv_ani[I2V_INV_ANI_GRIDS_1DK(ii_invr,i_grid)]) / (grid.r_loc_inv_ani[I2V_INV_ANI_GRIDS_1DK(ii_invr+1,i_grid)] - grid.r_loc_inv_ani[I2V_INV_ANI_GRIDS_1DK(ii_invr,i_grid)]);
                    break;
                }
            }
            // continue if r is out of the inversion grid
            if (r_r < 0 || r_r_ani < 0) continue;

            for (int j = j_start; j < j_end; j++) {
                CUSTOMREAL t_glob = grid.get_lat_min() + j*grid.get_delta_lat();
                r_t = -_1_CR;
                for (int ii_invt = 0; ii_invt < n_inv_J_loc-1; ii_invt++){
                    if ((t_glob - grid.t_loc_inv[I2V_INV_GRIDS_1DJ(ii_invt,i_grid)]) * (t_glob - grid.t_loc_inv[I2V_INV_GRIDS_1DJ(ii_invt+1,i_grid)]) <= 0) {
                        jdt = ii_invt;
                        r_t = (t_glob - grid.t_loc_inv[I2V_INV_GRIDS_1DJ(ii_invt,i_grid)]) / (grid.t_loc_inv[I2V_INV_GRIDS_1DJ(ii_invt+1,i_grid)] - grid.t_loc_inv[I2V_INV_GRIDS_1DJ(ii_invt,i_grid)]);
                        break;
                    }
                }

                r_t_ani = -_1_CR;
                for (int ii_invt = 0; ii_invt < n_inv_J_loc_ani-1; ii_invt++){
                    if ((t_glob - grid.t_loc_inv_ani[I2V_INV_ANI_GRIDS_1DJ(ii_invt,i_grid)]) * (t_glob - grid.t_loc_inv_ani[I2V_INV_ANI_GRIDS_1DJ(ii_invt+1,i_grid)]) <= 0) {
                        jdt_ani = ii_invt;
                        r_t_ani = (t_glob - grid.t_loc_inv_ani[I2V_INV_ANI_GRIDS_1DJ(ii_invt,i_grid)]) / (grid.t_loc_inv_ani[I2V_INV_ANI_GRIDS_1DJ(ii_invt+1,i_grid)] - grid.t_loc_inv_ani[I2V_INV_ANI_GRIDS_1DJ(ii_invt,i_grid)]);
                        break;
                    }
                }

                // continue if t is out of the inversion grid
                if (r_t < 0 || r_t_ani < 0) continue;

                for (int i = i_start; i < i_end; i++) {
                    CUSTOMREAL p_glob = grid.get_lon_min() + i*grid.get_delta_lon();
                    r_p = -_1_CR;
                    for (int ii_invp = 0; ii_invp < n_inv_I_loc-1; ii_invp++){
                        if ((p_glob - grid.p_loc_inv[I2V_INV_GRIDS_1DI(ii_invp,i_grid)]) * (p_glob - grid.p_loc_inv[I2V_INV_GRIDS_1DI(ii_invp+1,i_grid)]) <= 0) {
                            idp = ii_invp;
                            r_p = (p_glob - grid.p_loc_inv[I2V_INV_GRIDS_1DI(ii_invp,i_grid)]) / (grid.p_loc_inv[I2V_INV_GRIDS_1DI(ii_invp+1,i_grid)] - grid.p_loc_inv[I2V_INV_GRIDS_1DI(ii_invp,i_grid)]);
                            break;
                        }
                    }

                    r_p_ani = -_1_CR;
                    for (int ii_invp = 0; ii_invp < n_inv_I_loc_ani-1; ii_invp++){
                        if ((p_glob - grid.p_loc_inv_ani[I2V_INV_ANI_GRIDS_1DI(ii_invp,i_grid)]) * (p_glob - grid.p_loc_inv_ani[I2V_INV_ANI_GRIDS_1DI(ii_invp+1,i_grid)]) <= 0) {
                            idp_ani = ii_invp;
                            r_p_ani = (p_glob - grid.p_loc_inv_ani[I2V_INV_ANI_GRIDS_1DI(ii_invp,i_grid)]) / (grid.p_loc_inv_ani[I2V_INV_ANI_GRIDS_1DI(ii_invp+1,i_grid)] - grid.p_loc_inv_ani[I2V_INV_ANI_GRIDS_1DI(ii_invp,i_grid)]);
                            break;
                        }
                    }

                    // continue if p is out of the inversion grid
                    if (r_p < 0 || r_p_ani < 0) continue;

                    // global grid id to local id
                    int k_loc = k - grid.get_offset_k();
                    int j_loc = j - grid.get_offset_j();
                    int i_loc = i - grid.get_offset_i();

                    // check if *_loc are inside the local subdomain
                    if (k_loc < 0 || k_loc > loc_K-1) continue;
                    if (j_loc < 0 || j_loc > loc_J-1) continue;
                    if (i_loc < 0 || i_loc > loc_I-1) continue;

                    // update Ks_inv Keta_inv Kxi_inv
                    grid.Ks_inv_loc[  I2V_INV_KNL(idp,jdt,kdr)] += (_1_CR-r_r)*(_1_CR-r_t)*(_1_CR-r_p)*grid.Ks_loc[  I2V(i_loc,j_loc,k_loc)];
                    grid.Keta_inv_loc[I2V_INV_ANI_KNL(idp_ani,jdt_ani,kdr_ani)] += (_1_CR-r_r_ani)*(_1_CR-r_t_ani)*(_1_CR-r_p_ani)*grid.Keta_loc[I2V(i_loc,j_loc,k_loc)];
                    grid.Kxi_inv_loc[ I2V_INV_ANI_KNL(idp_ani,jdt_ani,kdr_ani)] += (_1_CR-r_r_ani)*(_1_CR-r_t_ani)*(_1_CR-r_p_ani)*grid.Kxi_loc[ I2V(i_loc,j_loc,k_loc)];

                    grid.Ks_inv_loc[  I2V_INV_KNL(idp,jdt,kdr+1)] += r_r*(_1_CR-r_t)*(_1_CR-r_p)*grid.Ks_loc[  I2V(i_loc,j_loc,k_loc)];
                    grid.Keta_inv_loc[I2V_INV_ANI_KNL(idp_ani,jdt_ani,kdr_ani+1)] += r_r_ani*(_1_CR-r_t_ani)*(_1_CR-r_p_ani)*grid.Keta_loc[I2V(i_loc,j_loc,k_loc)];
                    grid.Kxi_inv_loc[ I2V_INV_ANI_KNL(idp_ani,jdt_ani,kdr_ani+1)] += r_r_ani*(_1_CR-r_t_ani)*(_1_CR-r_p_ani)*grid.Kxi_loc[ I2V(i_loc,j_loc,k_loc)];

                    grid.Ks_inv_loc[  I2V_INV_KNL(idp,jdt+1,kdr)] += (_1_CR-r_r)*r_t*(_1_CR-r_p)*grid.Ks_loc[  I2V(i_loc,j_loc,k_loc)];
                    grid.Keta_inv_loc[I2V_INV_ANI_KNL(idp_ani,jdt_ani+1,kdr_ani)] += (_1_CR-r_r_ani)*r_t_ani*(_1_CR-r_p_ani)*grid.Keta_loc[I2V(i_loc,j_loc,k_loc)];
                    grid.Kxi_inv_loc[ I2V_INV_ANI_KNL(idp_ani,jdt_ani+1,kdr_ani)] += (_1_CR-r_r_ani)*r_t_ani*(_1_CR-r_p_ani)*grid.Kxi_loc[ I2V(i_loc,j_loc,k_loc)];

                    grid.Ks_inv_loc[  I2V_INV_KNL(idp+1,jdt,kdr)] += (_1_CR-r_r)*(_1_CR-r_t)*r_p*grid.Ks_loc[  I2V(i_loc,j_loc,k_loc)];
                    grid.Keta_inv_loc[I2V_INV_ANI_KNL(idp_ani+1,jdt_ani,kdr_ani)] += (_1_CR-r_r_ani)*(_1_CR-r_t_ani)*r_p_ani*grid.Keta_loc[I2V(i_loc,j_loc,k_loc)];
                    grid.Kxi_inv_loc[ I2V_INV_ANI_KNL(idp_ani+1,jdt_ani,kdr_ani)] += (_1_CR-r_r_ani)*(_1_CR-r_t_ani)*r_p_ani*grid.Kxi_loc[ I2V(i_loc,j_loc,k_loc)];

                    grid.Ks_inv_loc[  I2V_INV_KNL(idp,jdt+1,kdr+1)] += r_r*r_t*(_1_CR-r_p)*grid.Ks_loc[  I2V(i_loc,j_loc,k_loc)];
                    grid.Keta_inv_loc[I2V_INV_ANI_KNL(idp_ani,jdt_ani+1,kdr_ani+1)] += r_r_ani*r_t_ani*(_1_CR-r_p_ani)*grid.Keta_loc[I2V(i_loc,j_loc,k_loc)];
                    grid.Kxi_inv_loc[ I2V_INV_ANI_KNL(idp_ani,jdt_ani+1,kdr_ani+1)] += r_r_ani*r_t_ani*(_1_CR-r_p_ani)*grid.Kxi_loc[ I2V(i_loc,j_loc,k_loc)];

                    grid.Ks_inv_loc[  I2V_INV_KNL(idp+1,jdt,kdr+1)] += r_r*(_1_CR-r_t)*r_p*grid.Ks_loc[  I2V(i_loc,j_loc,k_loc)];
                    grid.Keta_inv_loc[I2V_INV_ANI_KNL(idp_ani+1,jdt_ani,kdr_ani+1)] += r_r_ani*(_1_CR-r_t_ani)*r_p_ani*grid.Keta_loc[I2V(i_loc,j_loc,k_loc)];
                    grid.Kxi_inv_loc[ I2V_INV_ANI_KNL(idp_ani+1,jdt_ani,kdr_ani+1)] += r_r_ani*(_1_CR-r_t_ani)*r_p_ani*grid.Kxi_loc[ I2V(i_loc,j_loc,k_loc)];

                    grid.Ks_inv_loc[  I2V_INV_KNL(idp+1,jdt+1,kdr)] += (_1_CR-r_r)*r_t*r_p*grid.Ks_loc[  I2V(i_loc,j_loc,k_loc)];
                    grid.Keta_inv_loc[I2V_INV_ANI_KNL(idp_ani+1,jdt_ani+1,kdr_ani)] += (_1_CR-r_r_ani)*r_t_ani*r_p_ani*grid.Keta_loc[I2V(i_loc,j_loc,k_loc)];
                    grid.Kxi_inv_loc[ I2V_INV_ANI_KNL(idp_ani+1,jdt_ani+1,kdr_ani)] += (_1_CR-r_r_ani)*r_t_ani*r_p_ani*grid.Kxi_loc[ I2V(i_loc,j_loc,k_loc)];

                    grid.Ks_inv_loc[  I2V_INV_KNL(idp+1,jdt+1,kdr+1)] += r_r*r_t*r_p*grid.Ks_loc[  I2V(i_loc,j_loc,k_loc)];
                    grid.Keta_inv_loc[I2V_INV_ANI_KNL(idp_ani+1,jdt_ani+1,kdr_ani+1)] += r_r_ani*r_t_ani*r_p_ani*grid.Keta_loc[I2V(i_loc,j_loc,k_loc)];
                    grid.Kxi_inv_loc[ I2V_INV_ANI_KNL(idp_ani+1,jdt_ani+1,kdr_ani+1)] += r_r_ani*r_t_ani*r_p_ani*grid.Kxi_loc[ I2V(i_loc,j_loc,k_loc)];

                } // end for i
            } // end for j
        } // end for k

        // rescale the kernel:  kernel -> kernel / (volume of the inversion grid)
        if (IP.get_invgrid_volume_rescale()){
            CUSTOMREAL volume_r = 1.0;
            CUSTOMREAL volume_t = 1.0;
            CUSTOMREAL volume_p = 1.0;
            CUSTOMREAL volume   = 1.0;
            for (int ii_invr = 0; ii_invr < n_inv_K_loc; ii_invr++){
                if (ii_invr == 0)
                    volume_r = abs(grid.r_loc_inv[I2V_INV_GRIDS_1DK(1,i_grid)] - grid.r_loc_inv[I2V_INV_GRIDS_1DK(0,i_grid)]);
                else if (ii_invr == n_inv_K_loc - 1)
                    volume_r = abs(grid.r_loc_inv[I2V_INV_GRIDS_1DK(n_inv_K_loc-1,i_grid)] - grid.r_loc_inv[I2V_INV_GRIDS_1DK(n_inv_K_loc-2,i_grid)]);
                else
                    volume_r = abs(grid.r_loc_inv[I2V_INV_GRIDS_1DK(ii_invr+1,i_grid)] - grid.r_loc_inv[I2V_INV_GRIDS_1DK(ii_invr-1,i_grid)]);

                for (int ii_invt = 0; ii_invt < n_inv_J_loc; ii_invt++){
                    if (ii_invt == 0)
                        volume_t = abs(grid.t_loc_inv[I2V_INV_GRIDS_1DJ(1,i_grid)] - grid.t_loc_inv[I2V_INV_GRIDS_1DJ(0,i_grid)]);
                    else if (ii_invt == n_inv_J_loc - 1)
                        volume_t = abs(grid.t_loc_inv[I2V_INV_GRIDS_1DJ(n_inv_J_loc-1,i_grid)] - grid.t_loc_inv[I2V_INV_GRIDS_1DJ(n_inv_J_loc-2,i_grid)]);
                    else
                        volume_t = abs(grid.t_loc_inv[I2V_INV_GRIDS_1DJ(ii_invt+1,i_grid)] - grid.t_loc_inv[I2V_INV_GRIDS_1DJ(ii_invt-1,i_grid)]);

                    for (int ii_invp = 0; ii_invp < n_inv_I_loc; ii_invp++){
                        if (ii_invp == 0)
                            volume_p = abs(grid.p_loc_inv[I2V_INV_GRIDS_1DI(1,i_grid)] - grid.p_loc_inv[I2V_INV_GRIDS_1DI(0,i_grid)]);  
                        else if (ii_invp == n_inv_I_loc - 1)
                            volume_p = abs(grid.p_loc_inv[I2V_INV_GRIDS_1DI(n_inv_I_loc-1,i_grid)] - grid.p_loc_inv[I2V_INV_GRIDS_1DI(n_inv_I_loc-2,i_grid)]);
                        else
                            volume_p = abs(grid.p_loc_inv[I2V_INV_GRIDS_1DI(ii_invp+1,i_grid)] - grid.p_loc_inv[I2V_INV_GRIDS_1DI(ii_invp-1,i_grid)]);
                        
                        volume = volume_r * volume_t * volume_p;
                        grid.Ks_inv_loc[    I2V_INV_KNL(ii_invp,ii_invt,ii_invr)]   *= volume;
                    }
                }
            }

            for (int ii_invr = 0; ii_invr < n_inv_K_loc_ani; ii_invr++){
                if (ii_invr == 0)
                    volume_r = abs(grid.r_loc_inv_ani[I2V_INV_ANI_GRIDS_1DK(1,i_grid)] - grid.r_loc_inv_ani[I2V_INV_ANI_GRIDS_1DK(0,i_grid)]);
                else if (ii_invr == n_inv_K_loc_ani - 1)
                    volume_r = abs(grid.r_loc_inv_ani[I2V_INV_ANI_GRIDS_1DK(n_inv_K_loc_ani-1,i_grid)] - grid.r_loc_inv_ani[I2V_INV_ANI_GRIDS_1DK(n_inv_K_loc_ani-2,i_grid)]);
                else
                    volume_r = abs(grid.r_loc_inv_ani[I2V_INV_ANI_GRIDS_1DK(ii_invr+1,i_grid)] - grid.r_loc_inv_ani[I2V_INV_ANI_GRIDS_1DK(ii_invr-1,i_grid)]);

                for (int ii_invt = 0; ii_invt < n_inv_J_loc_ani; ii_invt++){
                    if (ii_invt == 0)
                        volume_t = abs(grid.t_loc_inv_ani[I2V_INV_ANI_GRIDS_1DJ(1,i_grid)] - grid.t_loc_inv_ani[I2V_INV_ANI_GRIDS_1DJ(0,i_grid)]);
                    else if (ii_invt == n_inv_J_loc_ani - 1)
                        volume_t = abs(grid.t_loc_inv_ani[I2V_INV_ANI_GRIDS_1DJ(n_inv_J_loc_ani-1,i_grid)] - grid.t_loc_inv_ani[I2V_INV_ANI_GRIDS_1DJ(n_inv_J_loc_ani-2,i_grid)]);
                    else
                        volume_t = abs(grid.t_loc_inv_ani[I2V_INV_ANI_GRIDS_1DJ(ii_invt+1,i_grid)] - grid.t_loc_inv_ani[I2V_INV_ANI_GRIDS_1DJ(ii_invt-1,i_grid)]);

                    for (int ii_invp = 0; ii_invp < n_inv_I_loc_ani; ii_invp++){
                        if (ii_invp == 0)
                            volume_p = abs(grid.p_loc_inv_ani[I2V_INV_ANI_GRIDS_1DI(1,i_grid)] - grid.p_loc_inv_ani[I2V_INV_ANI_GRIDS_1DI(0,i_grid)]);  
                        else if (ii_invp == n_inv_I_loc_ani - 1)
                            volume_p = abs(grid.p_loc_inv_ani[I2V_INV_ANI_GRIDS_1DI(n_inv_I_loc_ani-1,i_grid)] - grid.p_loc_inv_ani[I2V_INV_ANI_GRIDS_1DI(n_inv_I_loc_ani-2,i_grid)]);
                        else
                            volume_p = abs(grid.p_loc_inv_ani[I2V_INV_ANI_GRIDS_1DI(ii_invp+1,i_grid)] - grid.p_loc_inv_ani[I2V_INV_ANI_GRIDS_1DI(ii_invp-1,i_grid)]);
                        
                        volume = volume_r * volume_t * volume_p;
                        grid.Keta_inv_loc[  I2V_INV_KNL(ii_invp,ii_invt,ii_invr)]   *= volume;
                        grid.Kxi_inv_loc[   I2V_INV_KNL(ii_invp,ii_invt,ii_invr)]   *= volume;
                    }
                }
            }
        }


        // sum over all sub-domains
        allreduce_cr_inplace(grid.Ks_inv_loc,   n_inv_I_loc*n_inv_J_loc*n_inv_K_loc);
        allreduce_cr_inplace(grid.Keta_inv_loc, n_inv_I_loc_ani*n_inv_J_loc_ani*n_inv_K_loc_ani);
        allreduce_cr_inplace(grid.Kxi_inv_loc,  n_inv_I_loc_ani*n_inv_J_loc_ani*n_inv_K_loc_ani);

        //
        // update factors
        //
        for (int k = k_start; k < k_end; k++) {
            CUSTOMREAL r_glob = grid.get_r_min() + k*grid.get_delta_r(); // global coordinate of r
            r_r = -_1_CR;
            for (int ii_invr = 0; ii_invr < n_inv_K_loc-1; ii_invr++){
                if ((r_glob - grid.r_loc_inv[I2V_INV_GRIDS_1DK(ii_invr,i_grid)]) * (r_glob - grid.r_loc_inv[I2V_INV_GRIDS_1DK(ii_invr+1,i_grid)]) <= 0) {
                    kdr = ii_invr;
                    r_r = (r_glob - grid.r_loc_inv[I2V_INV_GRIDS_1DK(ii_invr,i_grid)]) / (grid.r_loc_inv[I2V_INV_GRIDS_1DK(ii_invr+1,i_grid)] - grid.r_loc_inv[I2V_INV_GRIDS_1DK(ii_invr,i_grid)]);
                    break;
                }
            }

            r_r_ani = -_1_CR;
            for (int ii_invr = 0; ii_invr < n_inv_K_loc_ani-1; ii_invr++){
                // increasing or decreasing order
                if ((r_glob - grid.r_loc_inv_ani[I2V_INV_ANI_GRIDS_1DK(ii_invr,i_grid)]) * (r_glob - grid.r_loc_inv_ani[I2V_INV_ANI_GRIDS_1DK(ii_invr+1,i_grid)]) <= 0) {
                    kdr_ani = ii_invr;
                    r_r_ani = (r_glob - grid.r_loc_inv_ani[I2V_INV_ANI_GRIDS_1DK(ii_invr,i_grid)]) / (grid.r_loc_inv_ani[I2V_INV_ANI_GRIDS_1DK(ii_invr+1,i_grid)] - grid.r_loc_inv_ani[I2V_INV_ANI_GRIDS_1DK(ii_invr,i_grid)]);
                    break;
                }
            }
            // depth model update taper
            CUSTOMREAL depth = radius2depth(r_glob);
            if (depth < taper[0]) {     // weight = 0;
                weight = _0_CR;
            } else if (depth < taper[1]) {
                weight = (_1_CR - std::cos(PI*(depth - taper[0])/(taper[1] - taper[0]))) / _2_CR;
            } else {
                weight = _1_CR;
            }

            // continue if r is out of the inversion grid
            if (r_r < 0 || r_r_ani < 0) continue;

            for (int j = j_start; j < j_end; j++) {
                CUSTOMREAL t_glob = grid.get_lat_min() + j*grid.get_delta_lat();
                r_t = -_1_CR;
                for (int ii_invt = 0; ii_invt < n_inv_J_loc-1; ii_invt++){
                    if ((t_glob - grid.t_loc_inv[I2V_INV_GRIDS_1DJ(ii_invt,i_grid)]) * (t_glob - grid.t_loc_inv[I2V_INV_GRIDS_1DJ(ii_invt+1,i_grid)]) <= 0) {
                        jdt = ii_invt;
                        r_t = (t_glob - grid.t_loc_inv[I2V_INV_GRIDS_1DJ(ii_invt,i_grid)]) / (grid.t_loc_inv[I2V_INV_GRIDS_1DJ(ii_invt+1,i_grid)] - grid.t_loc_inv[I2V_INV_GRIDS_1DJ(ii_invt,i_grid)]);
                        break;
                    }
                }

                r_t_ani = -_1_CR;
                for (int ii_invt = 0; ii_invt < n_inv_J_loc_ani-1; ii_invt++){
                    if ((t_glob - grid.t_loc_inv_ani[I2V_INV_ANI_GRIDS_1DJ(ii_invt,i_grid)]) * (t_glob - grid.t_loc_inv_ani[I2V_INV_ANI_GRIDS_1DJ(ii_invt+1,i_grid)]) <= 0) {
                        jdt_ani = ii_invt;
                        r_t_ani = (t_glob - grid.t_loc_inv_ani[I2V_INV_ANI_GRIDS_1DJ(ii_invt,i_grid)]) / (grid.t_loc_inv_ani[I2V_INV_ANI_GRIDS_1DJ(ii_invt+1,i_grid)] - grid.t_loc_inv_ani[I2V_INV_ANI_GRIDS_1DJ(ii_invt,i_grid)]);
                        break;
                    }
                }

                // continue if t is out of the inversion grid
                if (r_t < 0 || r_t_ani < 0) continue;

                for (int i = i_start; i < i_end; i++) {
                    CUSTOMREAL p_glob = grid.get_lon_min() + i*grid.get_delta_lon();
                    r_p = -_1_CR;
                    for (int ii_invp = 0; ii_invp < n_inv_I_loc-1; ii_invp++){
                        if ((p_glob - grid.p_loc_inv[I2V_INV_GRIDS_1DI(ii_invp,i_grid)]) * (p_glob - grid.p_loc_inv[I2V_INV_GRIDS_1DI(ii_invp+1,i_grid)]) <= 0) {
                            idp = ii_invp;
                            r_p = (p_glob - grid.p_loc_inv[I2V_INV_GRIDS_1DI(ii_invp,i_grid)]) / (grid.p_loc_inv[I2V_INV_GRIDS_1DI(ii_invp+1,i_grid)] - grid.p_loc_inv[I2V_INV_GRIDS_1DI(ii_invp,i_grid)]);
                            break;
                        }
                    }

                    r_p_ani = -_1_CR;
                    for (int ii_invp = 0; ii_invp < n_inv_I_loc_ani-1; ii_invp++){
                        if ((p_glob - grid.p_loc_inv_ani[I2V_INV_ANI_GRIDS_1DI(ii_invp,i_grid)]) * (p_glob - grid.p_loc_inv_ani[I2V_INV_ANI_GRIDS_1DI(ii_invp+1,i_grid)]) <= 0) {
                            idp_ani = ii_invp;
                            r_p_ani = (p_glob - grid.p_loc_inv_ani[I2V_INV_ANI_GRIDS_1DI(ii_invp,i_grid)]) / (grid.p_loc_inv_ani[I2V_INV_ANI_GRIDS_1DI(ii_invp+1,i_grid)] - grid.p_loc_inv_ani[I2V_INV_ANI_GRIDS_1DI(ii_invp,i_grid)]);
                            break;
                        }
                    }

                    // continue if p is out of the inversion grid
                    if (r_p < 0 || r_p_ani < 0) continue;

                    // global grid id to local id
                    int k_loc = k - grid.get_offset_k();
                    int j_loc = j - grid.get_offset_j();
                    int i_loc = i - grid.get_offset_i();

                    // check if *_loc are inside the local subdomain
                    if (k_loc < 0 || k_loc > loc_K-1) continue;
                    if (j_loc < 0 || j_loc > loc_J-1) continue;
                    if (i_loc < 0 || i_loc > loc_I-1) continue;

                    CUSTOMREAL pert_Ks = 0.0;
                    CUSTOMREAL pert_Keta = 0.0;
                    CUSTOMREAL pert_Kxi = 0.0;

                    pert_Ks   += (_1_CR-r_r)*(_1_CR-r_t)*(_1_CR-r_p)*grid.Ks_inv_loc[  I2V_INV_KNL(idp,jdt,kdr)];
                    pert_Keta += (_1_CR-r_r_ani)*(_1_CR-r_t_ani)*(_1_CR-r_p_ani)*grid.Keta_inv_loc[I2V_INV_ANI_KNL(idp_ani,jdt_ani,kdr_ani)];
                    pert_Kxi  += (_1_CR-r_r_ani)*(_1_CR-r_t_ani)*(_1_CR-r_p_ani)*grid.Kxi_inv_loc[ I2V_INV_ANI_KNL(idp_ani,jdt_ani,kdr_ani)];

                    pert_Ks    += r_r*(_1_CR-r_t)*(_1_CR-r_p)*grid.Ks_inv_loc[  I2V_INV_KNL(idp,jdt,kdr+1)];
                    pert_Keta  += r_r_ani*(_1_CR-r_t_ani)*(_1_CR-r_p_ani)*grid.Keta_inv_loc[I2V_INV_ANI_KNL(idp_ani,jdt_ani,kdr_ani+1)];
                    pert_Kxi   += r_r_ani*(_1_CR-r_t_ani)*(_1_CR-r_p_ani)*grid.Kxi_inv_loc[ I2V_INV_ANI_KNL(idp_ani,jdt_ani,kdr_ani+1)];

                    pert_Ks    += (_1_CR-r_r)*r_t*(_1_CR-r_p)*grid.Ks_inv_loc[  I2V_INV_KNL(idp,jdt+1,kdr)];
                    pert_Keta  += (_1_CR-r_r_ani)*r_t_ani*(_1_CR-r_p_ani)*grid.Keta_inv_loc[I2V_INV_ANI_KNL(idp_ani,jdt_ani+1,kdr_ani)];
                    pert_Kxi   += (_1_CR-r_r_ani)*r_t_ani*(_1_CR-r_p_ani)*grid.Kxi_inv_loc[ I2V_INV_ANI_KNL(idp_ani,jdt_ani+1,kdr_ani)];

                    pert_Ks    += (_1_CR-r_r)*(_1_CR-r_t)*r_p*grid.Ks_inv_loc[  I2V_INV_KNL(idp+1,jdt,kdr)];
                    pert_Keta  += (_1_CR-r_r_ani)*(_1_CR-r_t_ani)*r_p_ani*grid.Keta_inv_loc[I2V_INV_ANI_KNL(idp_ani+1,jdt_ani,kdr_ani)];
                    pert_Kxi   += (_1_CR-r_r_ani)*(_1_CR-r_t_ani)*r_p_ani*grid.Kxi_inv_loc[ I2V_INV_ANI_KNL(idp_ani+1,jdt_ani,kdr_ani)];

                    pert_Ks    += r_r*r_t*(_1_CR-r_p)*grid.Ks_inv_loc[  I2V_INV_KNL(idp,jdt+1,kdr+1)];
                    pert_Keta  += r_r_ani*r_t_ani*(_1_CR-r_p_ani)*grid.Keta_inv_loc[I2V_INV_ANI_KNL(idp_ani,jdt_ani+1,kdr_ani+1)];
                    pert_Kxi   += r_r_ani*r_t_ani*(_1_CR-r_p_ani)*grid.Kxi_inv_loc[ I2V_INV_ANI_KNL(idp_ani,jdt_ani+1,kdr_ani+1)];

                    pert_Ks    += r_r*(_1_CR-r_t)*r_p*grid.Ks_inv_loc[  I2V_INV_KNL(idp+1,jdt,kdr+1)];
                    pert_Keta  += r_r_ani*(_1_CR-r_t_ani)*r_p_ani*grid.Keta_inv_loc[I2V_INV_ANI_KNL(idp_ani+1,jdt_ani,kdr_ani+1)];
                    pert_Kxi   += r_r_ani*(_1_CR-r_t_ani)*r_p_ani*grid.Kxi_inv_loc[ I2V_INV_ANI_KNL(idp_ani+1,jdt_ani,kdr_ani+1)];

                    pert_Ks    += (_1_CR-r_r)*r_t*r_p*grid.Ks_inv_loc[  I2V_INV_KNL(idp+1,jdt+1,kdr)];
                    pert_Keta  += (_1_CR-r_r_ani)*r_t_ani*r_p_ani*grid.Keta_inv_loc[I2V_INV_ANI_KNL(idp_ani+1,jdt_ani+1,kdr_ani)];
                    pert_Kxi   += (_1_CR-r_r_ani)*r_t_ani*r_p_ani*grid.Kxi_inv_loc[ I2V_INV_ANI_KNL(idp_ani+1,jdt_ani+1,kdr_ani)];

                    pert_Ks    += r_r*r_t*r_p*grid.Ks_inv_loc[  I2V_INV_KNL(idp+1,jdt+1,kdr+1)];
                    pert_Keta  += r_r_ani*r_t_ani*r_p_ani*grid.Keta_inv_loc[I2V_INV_ANI_KNL(idp_ani+1,jdt_ani+1,kdr_ani+1)];
                    pert_Kxi   += r_r_ani*r_t_ani*r_p_ani*grid.Kxi_inv_loc[ I2V_INV_ANI_KNL(idp_ani+1,jdt_ani+1,kdr_ani+1)];


                    // update para
                    grid.Ks_update_loc[  I2V(i_loc,j_loc,k_loc)] += weight * pert_Ks;
                    grid.Keta_update_loc[I2V(i_loc,j_loc,k_loc)] += weight * pert_Keta;
                    grid.Kxi_update_loc[ I2V(i_loc,j_loc,k_loc)] += weight * pert_Kxi;


                } // end for i
            } // end for j
        } // end for k

    } // end i_grid




}


#endif