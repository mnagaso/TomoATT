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


// information of source or receiver
class SrcRecInfo {
public:
    int id = -1;
    std::string name = "unknown";
    CUSTOMREAL dep; // stored as depth (km), so need to depth2radious function when using it in other part of the code
    CUSTOMREAL lat; // stored in degree, but convarted to radian when passed through get_src_* function
    CUSTOMREAL lon; // stored in degree, but convarted to radian when passed through get_src_* function

    int year             = -1;
    int month            = -1;
    int day              = -1;
    int hour             = -1;
    int min              = -1;
    CUSTOMREAL sec       = -1.0;
    CUSTOMREAL mag       = -1.0;

    CUSTOMREAL adjoint_source = 0.0;

    int n_data = 0;

    // arrays for storing arrival times on boundary surfaces, calculated by 2D Eikonal solver
    bool is_out_of_region = false;   // is the source or receiver in the region; false: in the refion; true: teleseismic
    CUSTOMREAL* arr_times_bound_N;   // arrival time of the receiver at the north boundary of the subdomain
    CUSTOMREAL* arr_times_bound_E;   // arrival time of the receiver at the east boundary of the subdomain
    CUSTOMREAL* arr_times_bound_W;   // arrival time of the receiver at the west boundary of the subdomain
    CUSTOMREAL* arr_times_bound_S;   // arrival time of the receiver at the south boundary of the subdomain
    CUSTOMREAL* arr_times_bound_Bot; // arrival time of the receiver at the bottom boundary of the subdomain
    bool*       is_bound_src;        // true if the source is on the boundary surface

    // kernel
    CUSTOMREAL sta_correct = 0.0;
    CUSTOMREAL sta_correct_kernel = 0.0;

    // relocation
    bool       is_stop      = false;
    CUSTOMREAL tau_opt      = 0.0;
    CUSTOMREAL grad_chi_i   = 0.0;
    CUSTOMREAL grad_chi_j   = 0.0;
    CUSTOMREAL grad_chi_k   = 0.0;
    CUSTOMREAL sum_weight   = 0.0;
    CUSTOMREAL vobj_src_reloc_old       = 99999999.9;
    CUSTOMREAL vobj_src_reloc           = 0.0;
    CUSTOMREAL vobj_grad_norm_src_reloc = 0.0;
    CUSTOMREAL DTi          = 0.0;
    CUSTOMREAL DTj          = 0.0;
    CUSTOMREAL DTk          = 0.0;
    CUSTOMREAL step_length_max  = step_length_src_reloc;  // 2 km default, step length for relocation
};

class DataInfo {
public:

    CUSTOMREAL data_weight = 1.0;   // the weight in the src_rec file
    CUSTOMREAL weight      = 1.0;   // the actual weight in the inversion, equal   data_weight * weight about the data type;

    std::string phase = "unknown";

    bool is_src_rec = false; // absolute traveltime, single source - receiver

    // source information
    int         id_src   = -1;
    std::string name_src = "unknown";

    // receiver information
    int         id_rec   = -1;
    std::string name_rec = "unknown";

    // traveltime
    CUSTOMREAL travel_time     = -999.0;
    CUSTOMREAL travel_time_obs = -999.0;

    bool is_rec_pair = false;   // common source differential traveltime

    // source infomation
    int         id_src_single   = -1;
    std::string name_src_single = "unknown";

    // receiver pair infomation
    std::vector<int>         id_rec_pair   = {-1,-1};
    std::vector<std::string> name_rec_pair = {"unknown","unknown"};

    // common source differential travel time
    CUSTOMREAL cs_dif_travel_time     = -999.0;
    CUSTOMREAL cs_dif_travel_time_obs = -999.0;

    // source pair infomation
    bool                     is_src_pair   = false;   // common receiver differential traveltime (for future)
    std::vector<int>         id_src_pair   = {-1,-1};
    std::vector<std::string> name_src_pair = {"unknown","unknown"};

    // receiver infomation
    int id_rec_single           = -1;
    std::string name_rec_single = "unknown";

    // common receiver differential travel time
    CUSTOMREAL cr_dif_travel_time     = -999.0;
    CUSTOMREAL cr_dif_travel_time_obs = -999.0;

};



// methods for managing SrcRec objects/lists

// parse src_rec_file
void parse_src_rec_file(std::string&                      , \
                        std::map<std::string, SrcRecInfo>&, \
                        std::map<std::string, SrcRecInfo>&, \
                        std::vector<DataInfo>&            , \
                        std::vector<std::string>&);

// parse sta_correctoion_file
void parse_sta_correction_file(std::string&, \
                               std::map<std::string, SrcRecInfo>&);

// swap the sources and receivers
void do_swap_src_rec(std::map<std::string, SrcRecInfo> &, \
                     std::map<std::string, SrcRecInfo> &, \
                     std::vector<DataInfo>             &);

#endif //SRC_REC_H