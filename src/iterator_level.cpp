#include "iterator_level.h"


Iterator_level::Iterator_level(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in) \
                : Iterator(IP, grid, src, io, first_init, is_teleseismic_in) {
    // do nothing
}


void Iterator_level::do_sweep_adj(int iswp, Grid& grid, InputParams& IP){
    if (subdom_main) {

        // set sweep direction
        set_sweep_direction(iswp);

        int r_start, t_start, p_start, r_end, t_end, p_end;

        // set loop range
        if (r_dirc < 0) {
            r_start = nr-2;
            r_end   = 0;
        } else {
            r_start = 1;
            r_end   = nr-1;
        }
        if (t_dirc < 0) {
            t_start = nt-2;
            t_end   = 0;
        } else {
            t_start = 1;
            t_end   = nt-1;
        }
        if (p_dirc < 0) {
            p_start = np-2;
            p_end   = 0;
        } else {
            p_start = 1;
            p_end   = np-1;
        }

        for (int kkr = r_start; kkr != r_end; kkr+=r_dirc) {
            for (int jjt = t_start; jjt != t_end; jjt+=t_dirc) {
                for (int iip = p_start; iip != p_end; iip+=p_dirc) {
                    // calculate stencils
                    calculate_stencil_adj(grid, iip, jjt, kkr);
                }
            }
        }

    } // end if subdom_main
}


Iterator_level_tele::Iterator_level_tele(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in) \
                : Iterator(IP, grid, src, io, first_init, is_teleseismic_in) {
    // do nothing
}


void Iterator_level_tele::do_sweep_adj(int iswp, Grid& grid, InputParams& IP){
    if (subdom_main) {

        // set sweep direction
        set_sweep_direction(iswp);

        int r_start, t_start, p_start, r_end, t_end, p_end;

        // set loop range
        if (r_dirc < 0) {
            r_start = nr-1;
            r_end   = -1;
        } else {
            r_start = 0;
            r_end   = nr;
        }
        if (t_dirc < 0) {
            t_start = nt-1;
            t_end   = -1;
        } else {
            t_start = 0;
            t_end   = nt;
        }
        if (p_dirc < 0) {
            p_start = np-1;
            p_end   = -1;
        } else {
            p_start = 0;
            p_end   = np;
        }

        for (int kkr = r_start; kkr != r_end; kkr+=r_dirc) {
            for (int jjt = t_start; jjt != t_end; jjt+=t_dirc) {
                for (int iip = p_start; iip != p_end; iip+=p_dirc) {
                    if (iip != 0    && jjt != 0    && kkr != 0 \
                     && iip != np-1 && jjt != nt-1 && kkr != nr-1) {
                        // calculate stencils
                        calculate_stencil_adj(grid, iip, jjt, kkr);
                    } else {
                        calculate_boundary_nodes_tele_adj(grid, iip, jjt, kkr);
                    }
                }
            }
        }
    } // end if subdom_main
}


Iterator_level_1st_order::Iterator_level_1st_order(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in) \
                         : Iterator_level(IP, grid, src, io, first_init, is_teleseismic_in) {
    // initialization is done in the base class
}


void Iterator_level_1st_order::do_sweep(int iswp, Grid& grid, InputParams& IP){
    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;

    for (int i_level = st_level; i_level <= ed_level; i_level++) {
        for (auto& ijk : ijk_for_this_subproc[i_level-st_level]) {
            int kk = ijk.at(2);
            int jj = ijk.at(1);
            int ii = ijk.at(0);

            if (r_dirc < 0) kkr = nr-kk; //kk-1;
            else            kkr = kk-1;  //nr-kk;
            if (t_dirc < 0) jjt = nt-jj; //jj-1;
            else            jjt = jj-1;  //nt-jj;
            if (p_dirc < 0) iip = np-ii; //ii-1;
            else            iip = ii-1;  //np-ii;

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


Iterator_level_3rd_order::Iterator_level_3rd_order(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in) \
                         : Iterator_level(IP, grid, src, io, first_init, is_teleseismic_in) {
    // initialization is done in the base class
}


void Iterator_level_3rd_order::do_sweep(int iswp, Grid& grid, InputParams& IP){
    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;

    for (int i_level = st_level; i_level <= ed_level; i_level++) {
        for (auto& ijk : ijk_for_this_subproc[i_level-st_level]) {
            int kk = ijk.at(2);
            int jj = ijk.at(1);
            int ii = ijk.at(0);

            if (r_dirc < 0) kkr = nr-kk; //kk-1;
            else            kkr = kk-1;  //nr-kk;
            if (t_dirc < 0) jjt = nt-jj; //jj-1;
            else            jjt = jj-1;  //nt-jj;
            if (p_dirc < 0) iip = np-ii; //ii-1;
            else            iip = ii-1;  //np-ii;

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


Iterator_level_1st_order_tele::Iterator_level_1st_order_tele(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in) \
                         : Iterator_level_tele(IP, grid, src, io, first_init, is_teleseismic_in) {
    // initialization is done in the base class
}


void Iterator_level_1st_order_tele::do_sweep(int iswp, Grid& grid, InputParams& IP){
    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;

    for (int i_level = st_level; i_level <= ed_level; i_level++) {
        for (auto& ijk : ijk_for_this_subproc[i_level-st_level]) {
            int kk = ijk.at(2);
            int jj = ijk.at(1);
            int ii = ijk.at(0);

            //std::cout << sub_rank << ": " << i_level << " " << ii << " " << jj << " " << kk << std::endl;

            if (r_dirc < 0) kkr = nr-kk; //kk-1;
            else            kkr = kk-1;  //nr-kk;
            if (t_dirc < 0) jjt = nt-jj; //jj-1;
            else            jjt = jj-1;  //nt-jj;
            if (p_dirc < 0) iip = np-ii; //ii-1;
            else            iip = ii-1; //np-ii;

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


Iterator_level_3rd_order_tele::Iterator_level_3rd_order_tele(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic_in) \
                         : Iterator_level_tele(IP, grid, src, io, first_init, is_teleseismic_in) {
    // initialization is done in the base class
}


void Iterator_level_3rd_order_tele::do_sweep(int iswp, Grid& grid, InputParams& IP){
    // set sweep direction
    set_sweep_direction(iswp);

    int iip, jjt, kkr;

    for (int i_level = st_level; i_level <= ed_level; i_level++) {
        for (auto& ijk : ijk_for_this_subproc[i_level-st_level]) {
            int kk = ijk.at(2);
            int jj = ijk.at(1);
            int ii = ijk.at(0);

            //std::cout << sub_rank << ": " << i_level << " " << ii << " " << jj << " " << kk << std::endl;

            if (r_dirc < 0) kkr = nr-kk; //kk-1;
            else            kkr = kk-1;  //nr-kk;
            if (t_dirc < 0) jjt = nt-jj; //jj-1;
            else            jjt = jj-1;  //nt-jj;
            if (p_dirc < 0) iip = np-ii; //ii-1;
            else            iip = ii-1; //np-ii;

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