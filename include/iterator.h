#ifndef ITERATOR_H
#define ITERATOR_H

#include <algorithm>
#include <initializer_list>
#include <vector>
#include "grid.h"
#include "input_params.h"
#include "utils.h"
#include "source.h"
#include "io.h"
#include "timer.h"

#ifdef USE_CUDA
#include "iterator_wrapper.cuh"
#endif

class Iterator {
public:
    Iterator(InputParams&, Grid&, Source&, IO_utils&, bool, bool);
    virtual ~Iterator();
    // regional source
    void run_iteration_forward(InputParams&, Grid&, IO_utils&, bool&); // run forward iteratiom till convergence
    void run_iteration_adjoint(InputParams&, Grid&, IO_utils&);        // run adjoint iteratiom till convergence

    // teleseismic source (teleseismic adjoint use the same function with reginal source)
    //void run_iteration_forward_teleseismic(InputParams&, Grid&, IO_utils&, bool&); // run forward iteratiom till convergence

    void initialize_arrays(InputParams&, Grid&, Source&); // initialize factors etc.

protected:
    void assign_processes_for_levels(); // assign intra-node processes for each sweeping level
    void set_sweep_direction(int);      // set sweep direction
    // regional source
    virtual void do_sweep(int, Grid&, InputParams&){};               // do sweeping with ordinal method
    void calculate_stencil_1st_order(Grid&, int&, int&, int&);     // calculate stencil for 1st order
    void calculate_stencil_3rd_order(Grid&, int&, int&, int&);     // calculate stencil for 3rd order
    void calculate_boundary_nodes(Grid&);                          // calculate boundary values
//    // teleseismic source
    void calculate_stencil_1st_order_tele(Grid&, int&, int&, int&); // calculate stencil for 1st order
    void calculate_stencil_3rd_order_tele(Grid&, int&, int&, int&); // calculate stencil for 3rd order
    void calculate_boundary_nodes_tele(Grid&, int&, int&, int&);    // calculate boundary values for teleseismic source
    void calculate_boundary_nodes_tele_adj(Grid&, int&, int&, int&);// calculate boundary values for teleseismic adjoint source

    // Hamiltonian calculation
    inline CUSTOMREAL calc_LF_Hamiltonian(Grid&, CUSTOMREAL& ,CUSTOMREAL& , \
                                                 CUSTOMREAL& ,CUSTOMREAL& , \
                                                 CUSTOMREAL& ,CUSTOMREAL&, \
                                                 int&, int&, int& );
    inline CUSTOMREAL calc_LF_Hamiltonian_tele(Grid&, CUSTOMREAL& ,CUSTOMREAL& , \
                                                 CUSTOMREAL& ,CUSTOMREAL& , \
                                                 CUSTOMREAL& ,CUSTOMREAL&, \
                                                 int&, int&, int& );


    // methods for adjoint field calculation
    void init_delta_and_Tadj(Grid&, InputParams&);                     // initialize delta and Tadj
    void fix_boundary_Tadj(Grid&);                                     // fix boundary values for Tadj
    virtual void do_sweep_adj(int, Grid&, InputParams&){};                       // do sweeping with ordinal method for adjoint field
    void calculate_stencil_adj(Grid&, int&, int&, int&);               // calculate stencil for 1st order for adjoint field

    // grid point information
    int* _nr, *_nt, *_np;           // number of grid points on the direction r, theta, phi
    CUSTOMREAL* _dr, *_dt, *_dp;    // grid spacing on the direction r, theta, phi
    int nr, nt, np;                 // number of grid points on the direction r, theta, phi
    CUSTOMREAL dr, dt, dp;          // grid spacing on the direction r, theta, phi
    int st_level;                   // start level for sweeping
    int ed_level;                   // end level for sweeping
    MPI_Win win_nr, win_nt, win_np; // windows for grid point information
    MPI_Win win_dr, win_dt, win_dp; // windows for grid point information

    std::vector< std::vector< std::vector<int> > > ijk_for_this_subproc; // ii,jj,kk for this process (level, kk)
    const int nswp = 8;          // number of sweeping directions
    int r_dirc, t_dirc, p_dirc;  // sweeping directions
    CUSTOMREAL sigr, sigt, sigp; //
    CUSTOMREAL coe;
    CUSTOMREAL wp1, pp1, wp2, pp2;
    CUSTOMREAL wt1, pt1, wt2, pt2;
    CUSTOMREAL wr1, pr1, wr2, pr2;
    CUSTOMREAL Htau, tpT;

    // iteration control
    int iter_count = 0;
    CUSTOMREAL ini_diff_L1 = HUGE_VAL, ini_diff_Linf = HUGE_VAL;
    CUSTOMREAL ini_err_L1  = HUGE_VAL, ini_err_Linf  = HUGE_VAL;
    CUSTOMREAL cur_diff_L1 = HUGE_VAL, cur_diff_Linf = HUGE_VAL;
    CUSTOMREAL cur_err_L1  = HUGE_VAL, cur_err_Linf  = HUGE_VAL;

    // teleseismic flag
    bool is_teleseismic = false;

};

// define derived classes for each iteration scheme


#endif // ITERATOR_H