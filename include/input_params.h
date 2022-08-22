#ifndef INPUT_PARAMS_H
#define INPUT_PARAMS_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <unordered_set>

#include "config.h"
#include "yaml-cpp/yaml.h"
#include "utils.h"
#include "mpi_funcs.h"


// strucutre for storing source or receiver information
class SrcRec {
public:

    int id_src = 0;
    int id_rec = -9999;
    int n_rec  = 0;

    CUSTOMREAL dep; // stored as depth (km), so need to depth2radious function when using it in other part of the code
    CUSTOMREAL lat; // stored in degree, but convarted to radian when passed through get_src_* function
    CUSTOMREAL lon; // stored in degree, but convarted to radian when passed through get_src_* function

    CUSTOMREAL arr_time;     // calculated arrival time will be stored and updated during the simulation
    CUSTOMREAL arr_time_ori; // recorded/original arrival time (written in the input file)
    CUSTOMREAL t_adj;        // adjoint source time = calculated (arr_time) - recorded (arr_time_ori)
    CUSTOMREAL weight=1.0;   // weight of the receiver

    int year             = 9999;
    int month            = 99;
    int day              = 99;
    int hour             = 99;
    int min              = 99;
    CUSTOMREAL sec       = 99;
    CUSTOMREAL mag       = 99;
    std::string id_event = "src_name_dummy";
    std::string name_rec = "rec_name_dummy";
    std::string phase    = "P";
    CUSTOMREAL epi_dist  = 0;

    // original ids of source/receiver before swapping
    std::vector<int> id_srcs_ori;

    //
    // variables for teleseismic source
    //
    bool is_teleseismic = false;
    // arrays for storing arrival times on boundary surfaces, calculated by 2D Eikonal solver
    CUSTOMREAL* arr_times_bound_N; // arrival time of the receiver at the north boundary of the subdomain
    CUSTOMREAL* arr_times_bound_E; // arrival time of the receiver at the east boundary of the subdomain
    CUSTOMREAL* arr_times_bound_W; // arrival time of the receiver at the west boundary of the subdomain
    CUSTOMREAL* arr_times_bound_S; // arrival time of the receiver at the south boundary of the subdomain
    CUSTOMREAL* arr_times_bound_Bot; // arrival time of the receiver at the bottom boundary of the subdomain
    bool*       is_bound_src; // true if the source is on the boundary surface
};


// InputParams class for reading/storing and outputing input parameters and src_rec information
class InputParams {

public:
    InputParams(std::string&);
    ~InputParams();

    // get parameters
    CUSTOMREAL get_min_dep(){return min_dep;};
    CUSTOMREAL get_max_dep(){return max_dep;};
    CUSTOMREAL get_min_lat(){return min_lat*DEG2RAD;};
    CUSTOMREAL get_max_lat(){return max_lat*DEG2RAD;};
    CUSTOMREAL get_min_lon(){return min_lon*DEG2RAD;};
    CUSTOMREAL get_max_lon(){return max_lon*DEG2RAD;};

    CUSTOMREAL get_min_dep_inv(){return min_dep_inv;};
    CUSTOMREAL get_max_dep_inv(){return max_dep_inv;};
    CUSTOMREAL get_min_lat_inv(){return min_lat_inv*DEG2RAD;};
    CUSTOMREAL get_max_lat_inv(){return max_lat_inv*DEG2RAD;};
    CUSTOMREAL get_min_lon_inv(){return min_lon_inv*DEG2RAD;};
    CUSTOMREAL get_max_lon_inv(){return max_lon_inv*DEG2RAD;};

    CUSTOMREAL           get_src_radius();
    CUSTOMREAL           get_src_lat();
    CUSTOMREAL           get_src_lon();
    std::string          get_src_rec_file(){return src_rec_file;};
    bool                 get_src_rec_file_exist(){return src_rec_file_exist;};
    SrcRec&              get_src_point(int);  // return SrcRec object
    std::vector<SrcRec>& get_rec_points(int); // return receivers for the current source

    CUSTOMREAL get_conv_tol(){return conv_tol;};
    CUSTOMREAL get_max_iter(){return max_iter;};

    int get_stencil_order(){return stencil_order;};
    int get_sweep_type()   {return sweep_type;};

    std::string get_init_model_path(){return init_model_path;};

    int get_do_inversion()    {return do_inversion;};
    int get_n_inversion_grid(){return n_inversion_grid;};
    int get_n_inv_r()         {return n_inv_r;};
    int get_n_inv_t()         {return n_inv_t;};
    int get_n_inv_p()         {return n_inv_p;};
    int get_max_iter_inv()    {return max_iter_inv;};

private:
    // boundary information
    CUSTOMREAL min_dep; // minimum depth in km
    CUSTOMREAL max_dep; // maximum depth in km
    CUSTOMREAL min_lat; // minimum latitude in degree
    CUSTOMREAL max_lat; // maximum latitude in degree
    CUSTOMREAL min_lon; // minimum longitude in degree
    CUSTOMREAL max_lon; // maximum longitude in degree
    // inversion grid
    CUSTOMREAL min_dep_inv=-99999; // minimum depth in km
    CUSTOMREAL max_dep_inv=-99999; // maximum depth in km
    CUSTOMREAL min_lat_inv=-99999; // minimum latitude
    CUSTOMREAL max_lat_inv=-99999; // maximum latitude
    CUSTOMREAL min_lon_inv=-99999; // minimum longitude
    CUSTOMREAL max_lon_inv=-99999; // maximum longitude


    // source information
    CUSTOMREAL src_dep;              // source depth in km
    CUSTOMREAL src_lat;              // source latitude in degrees
    CUSTOMREAL src_lon;              // source longitude in degrees
    std::string src_rec_file;        // source receiver file
    std::string src_rec_file_out;    // source receiver file to be output
    bool src_rec_file_exist = false; // source receiver file exist

    // model input files
    std::string init_model_type; // model type
    std::string init_model_path; // model file path init

    // inversion
    int do_inversion=0;                  // do inversion or not (0: no, 1: yes)
    int n_inversion_grid=1;              // number of inversion grid
    int n_inv_r=1, n_inv_t=1, n_inv_p=1; // number of inversion grid in r, t, p direction

    // convergence setting
    CUSTOMREAL conv_tol;       // convergence tolerance
    int        max_iter;       // maximum number of iterations
    int        max_iter_inv=1; // maximum number of iterations for inversion

    // calculation method
    int stencil_order; // stencil order
    int sweep_type; // sweep type (0: legacy, 1: cuthill-mckee, 2: cuthil-mckee with shm parallelization)

    // parse src_rec_file
    void parse_src_rec_file();
    // list for all of src point data
    std::vector<SrcRec> src_points;
    std::vector<SrcRec> src_points_back; // for backup purposes
    std::vector<SrcRec> src_points_out;  // temporal storage for output

    // list for all of rec point data
    std::vector< std::vector<SrcRec> > rec_points;
    std::vector< std::vector<SrcRec> > rec_points_back; // for backup purposes
    std::vector< std::vector<SrcRec> > rec_points_out;  // temporal storage for output
    // gather all arrival times to a main process
    void gather_all_arrival_times_to_main();
    // swap the sources and receivers
    void do_swap_src_rec();
    bool swap_src_rec = false;     // whether the src/rec are swapped or not
    void reverse_src_rec_points(); // reverse the swapped src/rec

    // tele seismic source management
    void separate_region_and_tele_src();               // check if the source is tele seismic or not
    void merge_region_and_tele_src();                  // merge tele seismic source to the region source
    std::vector<SrcRec> tele_src_points;               // tele seismic source points
    std::vector<std::vector<SrcRec> > tele_rec_points; // tele seismic receiver points
    bool i_first=false, i_last=false, j_first=false, j_last=false, k_first=false; // store info if this subdomain has outer boundary

public:
    void allocate_memory_tele_boundaries(int, int, int, int,
        bool, bool, bool, bool, bool); // allocate memory for tele boundaries
private:
    // check contradictions in input parameters
    void check_contradictions();

public:
    // prepare source list for this simulation group
    void prepare_src_list();
    // src points for this sim group
    std::vector<int> src_ids_this_sim;
    // write out src_rec_file
    void write_src_rec_file();

};



#endif