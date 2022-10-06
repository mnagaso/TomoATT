#include "iterator_level.h"


Iterator_level::Iterator_level(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                : Iterator(IP, grid, src, io, first_init, is_teleseismic_in, is_second_run_in) {
    // do nothing
}


void Iterator_level::do_sweep_adj(int iswp, Grid& grid, InputParams& IP){

    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;
    for (int i_level = st_level; i_level <= ed_level; i_level++) {
        for (auto& ijk : ijk_for_this_subproc[i_level-st_level]) {

            if (r_dirc < 0) kkr = nr-ijk.at(2); //kk-1;
            else            kkr = ijk.at(2)-1;  //nr-kk;
            if (t_dirc < 0) jjt = nt-ijk.at(1); //jj-1;
            else            jjt = ijk.at(1)-1;  //nt-jj;
            if (p_dirc < 0) iip = np-ijk.at(0); //ii-1;
            else            iip = ijk.at(0)-1;  //np-ii;

            //
            // calculate stencils
            //
            calculate_stencil_adj(grid, iip, jjt, kkr);
        } // end ijk

        // mpi synchronization
        synchronize_all_sub();

    } // end loop i_level

}


Iterator_level_tele::Iterator_level_tele(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                : Iterator(IP, grid, src, io, first_init, is_teleseismic_in, is_second_run_in) {
    // do nothing
}


void Iterator_level_tele::do_sweep_adj(int iswp, Grid& grid, InputParams& IP){
    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;

    for (int i_level = st_level; i_level <= ed_level; i_level++) {
        for (auto& ijk : ijk_for_this_subproc[i_level-st_level]) {

            if (r_dirc < 0) kkr = nr-ijk.at(2); //kk-1;
            else            kkr = ijk.at(2)-1;  //nr-kk;
            if (t_dirc < 0) jjt = nt-ijk.at(1); //jj-1;
            else            jjt = ijk.at(1)-1;  //nt-jj;
            if (p_dirc < 0) iip = np-ijk.at(0); //ii-1;
            else            iip = ijk.at(0)-1;  //np-ii;

            //
            // calculate stencils
            //
            if (iip != 0    && jjt != 0    && kkr != 0 \
             && iip != np-1 && jjt != nt-1 && kkr != nr-1) {
                // calculate stencils
                calculate_stencil_adj(grid, iip, jjt, kkr);
            } else {
                calculate_boundary_nodes_tele_adj(grid, iip, jjt, kkr);
            }
        } // end ijk

        // mpi synchronization
        synchronize_all_sub();

    } // end loop i_level

}


Iterator_level_1st_order::Iterator_level_1st_order(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                         : Iterator_level(IP, grid, src, io, first_init, is_teleseismic_in, is_second_run_in) {
    // initialization is done in the base class
}


void Iterator_level_1st_order::do_sweep(int iswp, Grid& grid, InputParams& IP){
    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;

    for (int i_level = st_level; i_level <= ed_level; i_level++) {
        for (auto& ijk : ijk_for_this_subproc[i_level-st_level]) {

            if (r_dirc < 0) kkr = nr-ijk.at(2); //kk-1;
            else            kkr = ijk.at(2)-1;  //nr-kk;
            if (t_dirc < 0) jjt = nt-ijk.at(1); //jj-1;
            else            jjt = ijk.at(1)-1;  //nt-jj;
            if (p_dirc < 0) iip = np-ijk.at(0); //ii-1;
            else            iip = ijk.at(0)-1;  //np-ii;

            //
            // calculate stencils
            //
            if (grid.is_changed[I2V(iip, jjt, kkr)]) {
                calculate_stencil_1st_order(grid, iip, jjt, kkr);
            } // is_changed == true
        } // end ijk

        // mpi synchronization
        synchronize_all_sub();

    } // end loop i_level

    // update boundary
    if (subdom_main) {
        calculate_boundary_nodes(grid);
    }
}


Iterator_level_3rd_order::Iterator_level_3rd_order(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                         : Iterator_level(IP, grid, src, io, first_init, is_teleseismic_in, is_second_run_in) {
    // initialization is done in the base class
}


void Iterator_level_3rd_order::do_sweep(int iswp, Grid& grid, InputParams& IP){
    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;

    for (int i_level = st_level; i_level <= ed_level; i_level++) {
        for (auto& ijk : ijk_for_this_subproc[i_level-st_level]) {

            if (r_dirc < 0) kkr = nr-ijk.at(2); //kk-1;
            else            kkr = ijk.at(2)-1;  //nr-kk;
            if (t_dirc < 0) jjt = nt-ijk.at(1); //jj-1;
            else            jjt = ijk.at(1)-1;  //nt-jj;
            if (p_dirc < 0) iip = np-ijk.at(0); //ii-1;
            else            iip = ijk.at(0)-1;  //np-ii;

            //
            // calculate stencils
            //
            if (grid.is_changed[I2V(iip, jjt, kkr)]) {
                calculate_stencil_3rd_order(grid, iip, jjt, kkr);
            } // is_changed == true
        } // end ijk

        // mpi synchronization
        synchronize_all_sub();

    } // end loop i_level

    // update boundary
    if (subdom_main) {
        calculate_boundary_nodes(grid);
    }

}


Iterator_level_1st_order_tele::Iterator_level_1st_order_tele(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                         : Iterator_level_tele(IP, grid, src, io, first_init, is_teleseismic_in, is_second_run_in) {
    // initialization is done in the base class
}


void Iterator_level_1st_order_tele::do_sweep(int iswp, Grid& grid, InputParams& IP){
    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;

    for (int i_level = st_level; i_level <= ed_level; i_level++) {
        for (auto& ijk : ijk_for_this_subproc[i_level-st_level]) {

            if (r_dirc < 0) kkr = nr-ijk.at(2); //kk-1;
            else            kkr = ijk.at(2)-1;  //nr-kk;
            if (t_dirc < 0) jjt = nt-ijk.at(1); //jj-1;
            else            jjt = ijk.at(1)-1;  //nt-jj;
            if (p_dirc < 0) iip = np-ijk.at(0); //ii-1;
            else            iip = ijk.at(0)-1;  //np-ii;

            //
            // calculate stencils
            //
            if (iip != 0 && iip != np-1 && jjt != 0 && jjt != nt-1 && kkr != 0 && kkr != nr-1) {
                // calculate stencils
                calculate_stencil_1st_order_tele(grid, iip, jjt, kkr);
            } else {
                // update boundary
                calculate_boundary_nodes_tele(grid, iip, jjt, kkr);
            }
        } // end ijk

        // mpi synchronization
        synchronize_all_sub();

    } // end loop i_level
}


Iterator_level_3rd_order_tele::Iterator_level_3rd_order_tele(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in, bool is_second_run_in) \
                         : Iterator_level_tele(IP, grid, src, io, first_init, is_teleseismic_in, is_second_run_in) {
    // initialization is done in the base class
}


void Iterator_level_3rd_order_tele::do_sweep(int iswp, Grid& grid, InputParams& IP){
    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;

    for (int i_level = st_level; i_level <= ed_level; i_level++) {
        for (auto& ijk : ijk_for_this_subproc[i_level-st_level]) {

            if (r_dirc < 0) kkr = nr-ijk.at(2); //kk-1;
            else            kkr = ijk.at(2)-1;  //nr-kk;
            if (t_dirc < 0) jjt = nt-ijk.at(1); //jj-1;
            else            jjt = ijk.at(1)-1;  //nt-jj;
            if (p_dirc < 0) iip = np-ijk.at(0); //ii-1;
            else            iip = ijk.at(0)-1;  //np-ii;

            //
            // calculate stencils
            //
            if (iip != 0 && iip != np-1 && jjt != 0 && jjt != nt-1 && kkr != 0 && kkr != nr-1) {
                // calculate stencils
                calculate_stencil_3rd_order_tele(grid, iip, jjt, kkr);
            } else {
                // update boundary
                calculate_boundary_nodes_tele(grid, iip, jjt, kkr);
            }
        } // end ijk

        // mpi synchronization
        synchronize_all_sub();

    } // end loop i_level

}