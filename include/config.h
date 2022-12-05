#ifndef CONFIG_H
#define CONFIG_H

#include <math.h>
#include <mpi.h>
#include <string>
#include <vector>
#include <limits>

// custom floating point accuracy
//#define CUSTOMREAL float
//#define MPI_CR MPI_FLOAT
#define CUSTOMREAL double
#define MPI_CR MPI_DOUBLE


#define MPI_DUMMY_TAG 1000

inline int loc_I, loc_J, loc_K;
inline int loc_I_vis, loc_J_vis, loc_K_vis;
inline int loc_I_excl_ghost, loc_J_excl_ghost, loc_K_excl_ghost;
inline int n_inv_grids;
inline int n_inv_I_loc, n_inv_J_loc, n_inv_K_loc;

// 3d indices to 1d index
#define I2V(A,B,C) ((C)*loc_I*loc_J + (B)*loc_I + A)
inline void V2I(const int& ijk, int& i, int& j, int& k) {
    k = ijk / (loc_I * loc_J);
    j = (ijk - k * loc_I * loc_J) / loc_I;
    i = ijk - k * loc_I * loc_J - j * loc_I;
}

#define I2V_INV_GRIDS(A,B,C,D) ((D)*n_inv_I_loc*n_inv_J_loc*n_inv_K_loc + (C)*n_inv_I_loc*n_inv_J_loc + (B)*n_inv_I_loc + A)
#define I2V_INV_KNL(A,B,C)     ((C)*n_inv_I_loc*n_inv_J_loc + (B)*n_inv_I_loc + A)
#define I2V_INV_GRIDS_1DK(A,B)  ((B)*n_inv_K_loc + (A))
#define I2V_INV_GRIDS_1DJ(A,B)  ((B)*n_inv_J_loc + (A))
#define I2V_INV_GRIDS_1DI(A,B)  ((B)*n_inv_I_loc + (A))

#define I2V_EXCL_GHOST(A,B,C)  ((C)* loc_I_excl_ghost   * loc_J_excl_ghost +    (B)* loc_I_excl_ghost    + A)
//#define I2V_ELM_CONN(A,B,C)   ((C)*(loc_I_excl_ghost-1)*(loc_J_excl_ghost-1) + (B)*(loc_I_excl_ghost-1) + A)
#define I2V_VIS(A,B,C)         ((C)* loc_I_vis   * loc_J_vis +    (B)* loc_I_vis    + A)
#define I2V_3D(A,B,C)          ((C)* loc_I_excl_ghost*loc_J_excl_ghost + (B)* loc_I_excl_ghost + A)
#define I2V_ELM_CONN(A,B,C)    ((C)*(loc_I_vis-1)*(loc_J_vis-1) + (B)*(loc_I_vis-1) + A)
#define I2VLBFGS(A,B,C,D)      ((A)*loc_I*loc_J*loc_K + (D)*loc_I*loc_J + (C)*loc_I + B)

// 2d indices to 1d index
#define IJ2V(A,B,C) (C*loc_I*loc_J + (B)*loc_I + A)
#define JK2V(A,B,C) (C*loc_J*loc_K + (B)*loc_J + A)
#define IK2V(A,B,C) (C*loc_I*loc_K + (B)*loc_I + A)

inline const CUSTOMREAL eps    = 1e-12;
inline const CUSTOMREAL epsAdj = 1e-6;

inline bool isZero(CUSTOMREAL x) {
    return fabs(x) < eps;
}

inline bool isZeroAdj(CUSTOMREAL x) {
    return fabs(x) < epsAdj;
}

// constants
inline const CUSTOMREAL PI        = 3.14159265358979323846264338327950288;
inline const CUSTOMREAL DEG2RAD   = PI/180.0;
inline const CUSTOMREAL RAD2DEG   = 180.0/PI;
inline const CUSTOMREAL M2KM      = 0.001;
inline const CUSTOMREAL _0_CR     = 0.0;
inline const CUSTOMREAL _0_5_CR   = 0.5;
inline const CUSTOMREAL _1_CR     = 1.0;
inline const CUSTOMREAL _1_5_CR   = 1.5;
inline const CUSTOMREAL _2_CR     = 2.0;
inline const CUSTOMREAL _3_CR     = 3.0;
inline const CUSTOMREAL _4_CR     = 4.0;
inline const CUSTOMREAL _8_CR     = 8.0;
inline const CUSTOMREAL TAU_INITIAL_VAL = 1.0;
inline const CUSTOMREAL TAU_INF_VAL     = 20.0;
inline const int        SWEEPING_COEFF  = 1.0; // coefficient for calculationg sigr/sigt/sigp for cpu

// ASCII output precision
typedef std::numeric_limits< CUSTOMREAL > dbl;
inline const int ASCII_OUTPUT_PRECISION = dbl::max_digits10;

// radious of the earth in km
inline CUSTOMREAL       R_earth        = 6371.0; // for compatibility with fortran code
inline const CUSTOMREAL GAMMA          = 0.0;
inline const CUSTOMREAL r_kermel_mask  = 40.0;
inline CUSTOMREAL       step_size_init = 0.01; // update step size limit
inline CUSTOMREAL       step_size_lbfgs;

// RUN MODE TYPE FLAG
inline const int ONLY_FORWARD        = 0;
inline const int DO_INVERSION        = 1;
inline const int SRC_RELOCATION      = 2;

// SWEEPING TYPE FLAG
inline const int SWEEP_TYPE_LEGACY = 0;
inline const int SWEEP_TYPE_LEVEL  = 1;
inline      bool hybrid_stencil_order = false; // if true, code at first run 1st order, then change to 3rd order (Dong Cui 2020)

// convert depth <-> radius
inline CUSTOMREAL depth2radius(CUSTOMREAL depth) {
    return R_earth - depth;
}
inline CUSTOMREAL radius2depth(CUSTOMREAL radius) {
    return R_earth - radius;
}

// parallel calculation parameters
// here we use so called "inline variables" for defining global variables
inline int ngrid_i     = 0; // number of grid points in i direction
inline int ngrid_j     = 0; // number of grid points in j direction
inline int ngrid_k     = 0; // number of grid points in k direction
inline int ngrid_i_inv = 0; // number of inversion grid points in i direction
inline int ngrid_j_inv = 0; // number of inversion grid points in j direction
inline int ngrid_k_inv = 0; // number of inversion grid points in k direction

// mpi parameters
inline int      world_nprocs;     // total number of processes (all groups)
inline int      world_rank;       // mpi rank of this process  (all groups)
inline int      sim_nprocs;       // number of processes in this simulation group
inline int      sim_rank;         // mpi rank of this process in this simulation group
inline int      inter_sim_rank;
inline int      sub_nprocs;       // total number of processes
inline int      sub_rank;         // mpi rank of this process
inline int      inter_sub_rank;   // mpi rank of this process in the inter-subdomain communicator
inline int      inter_sub_nprocs; // number of processes in the inter-subdomain communicator
inline int      nprocs;           // = n subdomains
inline int      myrank;           // = id subdomain
inline MPI_Comm sim_comm, inter_sim_comm, sub_comm, inter_sub_comm; // mpi communicator for simulation, inter-simulation, subdomain, and inter subdomains
inline int      n_sims           = 1; // number of mpi groups for simultaneous runs
inline int      n_procs_each_sim = 1; // number of processes in each simulation group
inline int      n_subdomains     = 1; // number of subdomains
inline int      n_subprocs       = 1; // number of sub processes in each subdomain
inline int      ndiv_i           = 1; // number of divisions in x direction
inline int      ndiv_j           = 1; // number of divisions in y direction
inline int      ndiv_k           = 1; // number of divisions in z direction
inline int      id_sim           = 0; // simultaneous run id  (not equal to src id)
inline int      id_sim_src       = 0; // id of current target source
inline int      id_subdomain     = 0; // subdomain id
inline bool     subdom_main      = false; // true if this process is main process in subdomain

// read input yaml file
inline std::string input_file="input_params.yaml";
// output directory name
inline std::string output_dir="./OUTPUT_FILES/";
// output file format
#define OUTPUT_FORMAT_HDF5 11111
#define OUTPUT_FORMAT_ASCII 22222
//#define OUTPUT_FORMAT_BINARY 33333 //#TODO: add
inline int output_format = OUTPUT_FORMAT_HDF5; // 0 - ascii, 1 - hdf5, 2 - binary

// smooth parameters
inline int smooth_method = 0; // 0: multi grid parametrization, 1: laplacian smoothing
inline CUSTOMREAL smooth_lp = 1.0;
inline CUSTOMREAL smooth_lt = 1.0;
inline CUSTOMREAL smooth_lr = 1.0;
inline const int GRADIENT_DESCENT    = 0;
inline const int HALVE_STEPPING_MODE = 1;
inline const int LBFGS_MODE          = 2;
inline int optim_method              = 0; // 0: gradient descent, 1: halve_stepping, 2: LBFGS
inline const CUSTOMREAL wolfe_c1     = 1e-4;
inline const CUSTOMREAL wolfe_c2     = 0.9;
inline const int        Mbfgs        = 5;            // number of gradients/models stored in memory
inline CUSTOMREAL       regularization_weight = 0.5; // regularization weight
inline int              max_sub_iterations    = 20;  // maximum number of sub-iterations
inline const CUSTOMREAL LBFGS_RELATIVE_STEP_SIZE = 0.3; // relative step size for the second and later iteration

// variables for test
inline bool if_test = false;

// verboose mode
inline bool if_verbose = false;
// inline bool if_verbose = true;

// if use gpu
inline int use_gpu = 0; // 0: no, 1: yes

// earthquake relocation
inline CUSTOMREAL step_length_src_reloc = 0.0001; // step length for source relocation
inline const int  N_ITER_MAX_SRC_RELOC  = 1000;
inline const CUSTOMREAL TOL_SRC_RELOC   = 1e-2;

#endif // CONFIG_H