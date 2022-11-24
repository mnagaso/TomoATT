#ifndef KERNEL_H
#define KERNEL_H

#include "config.h"
#include "grid.h"
#include "input_params.h"

void calculate_sensitivity_kernel(Grid& grid, InputParams& IP){
    // calculate sensitivity kernel

    // kernel calculation will be done only by the subdom_main
    if (subdom_main) {
        // get the necessary parameters
        int np             = loc_I;
        int nt             = loc_J;
        int nr             = loc_K;
        CUSTOMREAL dr      = grid.dr;
        CUSTOMREAL dt      = grid.dt;
        CUSTOMREAL dp      = grid.dp;
        CUSTOMREAL src_lon = IP.get_src_lon();
        CUSTOMREAL src_lat = IP.get_src_lat();
        CUSTOMREAL src_r   = IP.get_src_radius();

        CUSTOMREAL weight   = _1_CR;
        // CUSTOMREAL * taper  = IP.get_kernel_taper();

        // inner points
        for (int kkr = 1; kkr < nr-1; kkr++) {
            for (int jjt = 1; jjt < nt-1; jjt++) {
                for (int iip = 1; iip < np-1; iip++) {

                    // mask within one grid around the source
                    if (std::abs(grid.r_loc_1d[kkr]-src_r)   >= dr \
                     || std::abs(grid.t_loc_1d[jjt]-src_lat) >= dt \
                     || std::abs(grid.p_loc_1d[iip]-src_lon) >= dp) {
                        // Kernel w r t slowness s
                        grid.Ks_loc[I2V(iip,jjt,kkr)] += weight * grid.Tadj_loc[I2V(iip,jjt,kkr)] * my_square(grid.fun_loc[I2V(iip,jjt,kkr)]);
                    }
                }
            }
        }

        // boundary
        for (int kkr = 0; kkr < nr; kkr++) {
            for (int jjt = 0; jjt < nt; jjt++) {
                // set Ks Kxi Keta to zero
                if (grid.i_first()){
                    grid.Ks_loc[I2V(0,jjt,kkr)]      = _0_CR;
                }
                if (grid.i_last()){
                    grid.Ks_loc[I2V(np-1,jjt,kkr)]   = _0_CR;
                }
           }
        }
        for (int kkr = 0; kkr < nr; kkr++) {
            for (int iip = 0; iip < np; iip++) {
                // set Ks Kxi Keta to zero
                if (grid.j_first()){
                    grid.Ks_loc[I2V(iip,0,kkr)]      = _0_CR;
                }
                if (grid.j_last()){
                    grid.Ks_loc[I2V(iip,nt-1,kkr)]   = _0_CR;
                }
            }
        }
        for (int jjt = 0; jjt < nt; jjt++) {
            for (int iip = 0; iip < np; iip++) {
                // set Ks Kxi Keta to zero
                if (grid.k_first()){
                    grid.Ks_loc[I2V(iip,jjt,0)]      = _0_CR;
                }
                if (grid.k_last()){
                    grid.Ks_loc[I2V(iip,jjt,nr-1)]   = _0_CR;
                }
            }
        }

    } // end if subdom_main
}


void sumup_kernels(Grid& grid) {
    if(subdom_main){
        int n_grids = loc_I*loc_J*loc_K;

        allreduce_cr_sim_inplace(grid.Ks_loc, n_grids);

        // share the values on boundary
        grid.send_recev_boundary_data(grid.Ks_loc);
    }
}

#endif