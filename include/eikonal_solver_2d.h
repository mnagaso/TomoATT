#ifndef EIKONAL_SOLVER_2D_H
#define EIKONAL_SOLVER_2D_H

#include <fstream>
#include <string>
#include <iostream>

#include "utils.h"
#include "input_params.h"
#include "grid.h"
#include "1d_models.h"
#include "mpi_funcs.h"
#include "io.h"


class PlainGrid {

public:
    PlainGrid(Source& src, InputParams& IP);
    ~PlainGrid();

    void run_iteration(InputParams& IP);
    void run_iteration_upwind(InputParams& IP);

//private:
    // grid parameters
    CUSTOMREAL src_r, src_t, src_p;
    CUSTOMREAL src_t_dummy;
    CUSTOMREAL rmin_3d, rmax_3d, tmin_3d, tmax_3d, pmin_3d, pmax_3d;
    CUSTOMREAL max_degree;
    int nr_2d, nt_2d;

    // for eikonal solver
    int count_cand;
    int ii,ii_nr,ii_n2r,ii_pr,ii_p2r,ii_nt,ii_n2t,ii_pt,ii_p2t;
    CUSTOMREAL ar,br,at,bt,ar1,ar2,at1,at2,br1,br2,bt1,bt2;
    CUSTOMREAL eqn_a, eqn_b, eqn_c, eqn_Delta;
    CUSTOMREAL tmp_tau, T_r, T_t, charact_r, charact_t;
    bool is_causality;
    std::vector<CUSTOMREAL> canditate = std::vector<CUSTOMREAL>(60);

    // arrays
    // azimuth
    CUSTOMREAL *azimuth_2d;
    // epicentral distance
    CUSTOMREAL *epicentral_distance_2d;
    // activated boundaries
    //bool *activated_boundaries;
    // grid
    CUSTOMREAL *r_2d, *t_2d; // coordinates
    CUSTOMREAL *fun_2d, *fac_a_2d, *fac_b_2d, *u_2d, *T_2d; // functions
    CUSTOMREAL *T0v_2d, *T0r_2d, *T0t_2d, *tau_2d, *tau_old_2d;
    bool *is_changed_2d;

    // member functions
    void allocate_arrays();
    void deallocate_arrays();
    void load_1d_model(std::string);
    void select_1d_model(std::string);
    void calculate_stencil_2d(int&, int&);
    void calculate_stencil_upwind_2d(int&, int&);
    CUSTOMREAL eps_2d = 1e-12;

private:

    // # TODO: make id model type selectable by input file
    // 1D model used for 2D solver
    std::vector<std::vector<CUSTOMREAL>> model_1d_ref = model_1d_prem;


};

void prepare_teleseismic_boundary_conditions(InputParams&, Grid&, IO_utils&);
void run_2d_solver(InputParams&, Source&, IO_utils&);
void interp2d(PlainGrid&, CUSTOMREAL, CUSTOMREAL, CUSTOMREAL&);
void load_2d_traveltime(InputParams&, Source&, Grid&, IO_utils&);
std::string get_2d_tt_filename(const std::string&, Source&);

#endif // EIKONAL_SOLVER_2D_H