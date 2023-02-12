#ifndef SRC_REC_H
#define SRC_REC_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>

#include "config.h"
#include "utils.h"
#include "timer.h"
#include "mpi_funcs.h"


// strucutre for storing source or receiver information
class SrcRec {
public:
    //
    // default values, for the case of single source - receiver
    //
    int id_src = 0;
    int id_rec = -9999;
    int n_rec  = 0;
    int n_data = 0;  // = n_rec + n_rec_pair

    CUSTOMREAL dep = 0.0; // stored as depth (km), so need to depth2radious function when using it in other part of the code
    CUSTOMREAL lat = 0.0; // stored in degree, but convarted to radian when passed through get_src_* function
    CUSTOMREAL lon = 0.0; // stored in degree, but convarted to radian when passed through get_src_* function

    CUSTOMREAL arr_time     = 0.0; // calculated arrival time will be stored and updated during the simulation
    CUSTOMREAL arr_time_ori = 0.0; // recorded/original arrival time (written in the input file)
    CUSTOMREAL t_adj        = 0.0; // adjoint source time = calculated (arr_time) - recorded (arr_time_ori)
    CUSTOMREAL weight       = 1.0; // weight

    //
    // another case of receiver pair (now only for teleseismicity)
    //
    bool is_rec_pair = false;
    int n_rec_pair = 0;

    std::vector<int> id_rec_pair{std::vector<int>(2)};
    std::vector<CUSTOMREAL> dep_pair{std::vector<CUSTOMREAL>(2)};
    std::vector<CUSTOMREAL> lat_pair{std::vector<CUSTOMREAL>(2)};
    std::vector<CUSTOMREAL> lon_pair{std::vector<CUSTOMREAL>(2)};
    std::vector<CUSTOMREAL> ddt_adj_pair{std::vector<CUSTOMREAL>(2)}; // adjoint source time = [calculated (dif_arr_time) - recorded (dif_arr_time_ori)] * 1 (for 1) or * -1 (for 2)
    std::vector<std::string> name_rec_pair = std::vector<std::string>(2);
    std::vector<CUSTOMREAL> station_correction_pair = {0.0, 0.0}; // station correction value (only use for teleseismic data)

    CUSTOMREAL dif_arr_time     = 0.0; // calculated differential arrival time will be stored and updated during the simulation,  arr_time1 - arr_time2
    CUSTOMREAL dif_arr_time_ori = 0.0; // recorded/original differential arrival time (written in the input file)


    // common parameters for both cases

    int year             = 9999;
    int month            = 99;
    int day              = 99;
    int hour             = 99;
    int min              = 99;
    CUSTOMREAL sec       = 99;
    CUSTOMREAL mag       = 99;
    std::string name_src = "src_name_dummy";
    std::string name_rec = "rec_name_dummy";
    std::string phase    = "P";

    // original ids of source/receiver before swapping
    std::vector<int> id_srcs_ori;
    // used for swapping case at writing process
    int id_src_ori, id_rec_ori;

    bool is_teleseismic = false;
    // arrays for storing arrival times on boundary surfaces, calculated by 2D Eikonal solver
    CUSTOMREAL* arr_times_bound_N  = nullptr; // arrival time of the receiver at the north boundary of the subdomain
    CUSTOMREAL* arr_times_bound_E  = nullptr; // arrival time of the receiver at the east boundary of the subdomain
    CUSTOMREAL* arr_times_bound_W  = nullptr; // arrival time of the receiver at the west boundary of the subdomain
    CUSTOMREAL* arr_times_bound_S  = nullptr; // arrival time of the receiver at the south boundary of the subdomain
    CUSTOMREAL* arr_times_bound_Bot= nullptr; // arrival time of the receiver at the bottom boundary of the subdomain
    bool*       is_bound_src       = nullptr; // true if the source is on the boundary surface

    // params for source relocation
    CUSTOMREAL DTk, DTj, DTi;                      // gradient of traveltime
    CUSTOMREAL tau_opt;                            // optimal origin time
    CUSTOMREAL sum_weight;                         // sum of weights of all sources
    CUSTOMREAL grad_chi_k, grad_chi_j, grad_chi_i; // gradient of objective function
    int        id_unique_list;                     // id of the source in the unique list
    CUSTOMREAL vobj_src_reloc;                     // value of objective function
    CUSTOMREAL vobj_grad_norm_src_reloc;           // norm of gradient of objective function

};

// methods for managing SrcRec objects/lists

// parse src_rec_file
void parse_src_rec_file(std::string&                      , \
                        std::vector<SrcRec>&              , \
                        std::vector<std::vector<SrcRec>>& , \
                        std::map<std::string, SrcRec>&    , \
                        std::map<std::string, CUSTOMREAL>&, \
                        std::map<std::string, CUSTOMREAL>&);


// swap the sources and receivers
void do_swap_src_rec(std::vector<SrcRec>& src_points,      std::vector<std::vector<SrcRec>>& rec_points, \
                     std::vector<SrcRec>& src_points_back, std::vector<std::vector<SrcRec>>& rec_points_back);

// reverse the swapped src/rec
void reverse_src_rec_points(std::vector<SrcRec>& src_points,      std::vector<std::vector<SrcRec>>& rec_points, \
                            std::vector<SrcRec>& src_points_back, std::vector<std::vector<SrcRec>>& rec_points_back, \
                            std::vector<SrcRec>& src_points_out,  std::vector<std::vector<SrcRec>>& rec_points_out, \
                            const bool& swap_src_rec, const int& run_mode); // reverse the swapped src/rec

// writeout src_rec_file
void writeout_src_rec_file(std::string&                      src_rec_file, \
                           std::vector<SrcRec>&              src_points_out, \
                           std::vector<std::vector<SrcRec>>& rec_points_out);

#endif //SRC_REC_H