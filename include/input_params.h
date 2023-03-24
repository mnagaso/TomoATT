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
#include <map>

#include "config.h"
#include "yaml-cpp/yaml.h"
#include "utils.h"
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

    CUSTOMREAL dep; // stored as depth (km), so need to depth2radious function when using it in other part of the code
    CUSTOMREAL lat; // stored in degree, but convarted to radian when passed through get_src_* function
    CUSTOMREAL lon; // stored in degree, but convarted to radian when passed through get_src_* function

    CUSTOMREAL arr_time = 0.0;     // calculated arrival time will be stored and updated during the simulation
    CUSTOMREAL arr_time_ori = 0.0; // recorded/original arrival time (written in the input file)
    CUSTOMREAL t_adj;        // adjoint source time = calculated (arr_time) - recorded (arr_time_ori)
    CUSTOMREAL weight=1.0;   // weight

    //
    // another case of receiver pair (now only for teleseismicity)
    //
    bool is_src_rec = false;    // absolute traveltime, single source - receiver
    bool is_rec_pair = false;   // common source differential traveltime 
    bool is_src_pair = false;   // common receiver differential traveltime (for future)

    int n_rec_pair = 0;

    std::vector<int> id_rec_pair         = std::vector<int>(2);
    std::vector<CUSTOMREAL> dep_pair     = std::vector<CUSTOMREAL>(2);
    std::vector<CUSTOMREAL> lat_pair     = std::vector<CUSTOMREAL>(2);
    std::vector<CUSTOMREAL> lon_pair     = std::vector<CUSTOMREAL>(2);
    std::vector<CUSTOMREAL> ddt_adj_pair = std::vector<CUSTOMREAL>(2);    // adjoint source time = [calculated (dif_arr_time) - recorded (dif_arr_time_ori)] * 1 (for 1) or * -1 (for 2)
    std::vector<std::string> name_rec_pair    = std::vector<std::string>(2);
    std::vector<CUSTOMREAL> station_correction_pair = {0.0, 0.0}; // station correction value (only use for teleseismic data)

    // name_rec_pair[0] = "rec1_name_dummy";
    // name_rec_pair[1] = "rec2_name_dummy";

    CUSTOMREAL dif_arr_time;     // calculated differential arrival time will be stored and updated during the simulation,  arr_time1 - arr_time2
    CUSTOMREAL dif_arr_time_ori; // recorded/original differential arrival time (written in the input file)


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
    int id_rec_ori; // used for swapping case at writing process

    bool is_teleseismic = false;
    // arrays for storing arrival times on boundary surfaces, calculated by 2D Eikonal solver
    CUSTOMREAL* arr_times_bound_N; // arrival time of the receiver at the north boundary of the subdomain
    CUSTOMREAL* arr_times_bound_E; // arrival time of the receiver at the east boundary of the subdomain
    CUSTOMREAL* arr_times_bound_W; // arrival time of the receiver at the west boundary of the subdomain
    CUSTOMREAL* arr_times_bound_S; // arrival time of the receiver at the south boundary of the subdomain
    CUSTOMREAL* arr_times_bound_Bot; // arrival time of the receiver at the bottom boundary of the subdomain
    bool*       is_bound_src; // true if the source is on the boundary surface

    // params for source relocation
    CUSTOMREAL DTk, DTj, DTi;  // gradient of traveltime
    CUSTOMREAL tau_opt;        // optimal origin time
    CUSTOMREAL sum_weight;     // sum of weights of all sources
    CUSTOMREAL grad_chi_k, grad_chi_j, grad_chi_i; // gradient of objective function
    int        id_unique_list; // id of the source in the unique list
    CUSTOMREAL vobj_src_reloc; // value of objective function
    CUSTOMREAL vobj_grad_norm_src_reloc; // norm of gradient of objective function
};

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
    bool is_out_of_region = false;     // is the source or receiver in the region; false: in the refion; true: teleseismic
    CUSTOMREAL* arr_times_bound_N; // arrival time of the receiver at the north boundary of the subdomain
    CUSTOMREAL* arr_times_bound_E; // arrival time of the receiver at the east boundary of the subdomain
    CUSTOMREAL* arr_times_bound_W; // arrival time of the receiver at the west boundary of the subdomain
    CUSTOMREAL* arr_times_bound_S; // arrival time of the receiver at the south boundary of the subdomain
    CUSTOMREAL* arr_times_bound_Bot; // arrival time of the receiver at the bottom boundary of the subdomain
    bool*       is_bound_src; // true if the source is on the boundary surface

    // kernel
    CUSTOMREAL sta_correct = 0.0;
    CUSTOMREAL sta_correct_kernel = 0.0;

    // relocation
    bool       is_stop      = false;
    CUSTOMREAL tau_opt      = 0.0;
    CUSTOMREAL grad_tau     = 0.0;
    CUSTOMREAL change_ortime = 0.0;

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
    CUSTOMREAL change_dep = 0.0;
    CUSTOMREAL change_lat = 0.0;
    CUSTOMREAL change_lon = 0.0;
};

class DataInfo {
public:

    CUSTOMREAL data_weight = 1.0;   // the weight in the src_rec file
    CUSTOMREAL weight      = 1.0;   // the actual weight in the inversion, equal   data_weight * weight about the data type;

    std::string phase = "unknown";

    bool is_src_rec = false;    // absolute traveltime, single source - receiver
    // source information
    int id_src = -1;
    std::string name_src = "unknown";
    // receiver information
    int id_rec = -1;
    std::string name_rec = "unknown";
    // traveltime
    CUSTOMREAL travel_time = -999.0;
    CUSTOMREAL travel_time_obs = -999.0;

    bool is_rec_pair = false;   // common source differential traveltime 
    // source infomation
    int id_src_single = -1; 
    std::string name_src_single = "unknown";
    // receiver pair infomation
    std::vector<int> id_rec_pair = {-1,-1};
    std::vector<std::string> name_rec_pair = {"unknown","unknown"};
    // common source differential travel time
    CUSTOMREAL cs_dif_travel_time = -999.0; 
    CUSTOMREAL cs_dif_travel_time_obs = -999.0; 

    bool is_src_pair = false;   // common receiver differential traveltime (for future)
    // source pair infomation
    std::vector<int> id_src_pair = {-1,-1};
    std::vector<std::string> name_src_pair = {"unknown","unknown"};
    // receiver infomation
    int id_rec_single = -1; 
    std::string name_rec_single = "unknown";
    // common receiver differential travel time
    CUSTOMREAL cr_dif_travel_time = -999.0;
    CUSTOMREAL cr_dif_travel_time_obs = -999.0;

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

    // CUSTOMREAL           get_src_radius();
    // CUSTOMREAL           get_src_lat();
    // CUSTOMREAL           get_src_lon();
    CUSTOMREAL           get_src_radius_nv();
    CUSTOMREAL           get_src_lat_nv();
    CUSTOMREAL           get_src_lon_nv();
    std::string          get_src_rec_file()      {return src_rec_file;};
    bool                 get_src_rec_file_exist(){return src_rec_file_exist;};
    // SrcRec&              get_src_point(int);  // return SrcRec object
    // std::vector<SrcRec>& get_rec_points(int); // return receivers for the current source
    SrcRecInfo&                 get_src_point_nv(std::string);  // return SrcRec object
    SrcRecInfo&                 get_rec_point_nv(std::string); // return receivers for the current source
    std::vector<std::string>    get_rec_points_nv(std::string);  // return SrcRec object
    

    CUSTOMREAL get_conv_tol(){return conv_tol;};
    void set_conv_tol(CUSTOMREAL conv_tol_){conv_tol = conv_tol_;};
    CUSTOMREAL get_max_iter(){return max_iter;};

    int get_stencil_order(){return stencil_order;};
    void set_stencil_order(int stencil_order_){stencil_order = stencil_order_;};
    int get_stencil_type(){return stencil_type;};
    int get_sweep_type()   {return sweep_type;};

    std::string get_init_model_path(){return init_model_path;};
    std::string get_model_1d_name()  {return model_1d_name;};

    int get_run_mode()        {return run_mode;};
    int get_n_inversion_grid(){return n_inversion_grid;};

    int get_type_dep_inv()         {return type_dep_inv;};
    int get_type_lat_inv()         {return type_lat_inv;};
    int get_type_lon_inv()         {return type_lon_inv;};
    // type = 0:
    int get_n_inv_r()         {return n_inv_r;};
    int get_n_inv_t()         {return n_inv_t;};
    int get_n_inv_p()         {return n_inv_p;};
    CUSTOMREAL get_min_dep_inv(){return min_dep_inv;};
    CUSTOMREAL get_max_dep_inv(){return max_dep_inv;};
    CUSTOMREAL get_min_lat_inv(){return min_lat_inv*DEG2RAD;};
    CUSTOMREAL get_max_lat_inv(){return max_lat_inv*DEG2RAD;};
    CUSTOMREAL get_min_lon_inv(){return min_lon_inv*DEG2RAD;};
    CUSTOMREAL get_max_lon_inv(){return max_lon_inv*DEG2RAD;};

    // type = 1:
    int get_n_inv_r_flex()         {return n_inv_r_flex;};
    int get_n_inv_t_flex()         {return n_inv_t_flex;};
    int get_n_inv_p_flex()         {return n_inv_p_flex;};
    CUSTOMREAL * get_dep_inv(){return dep_inv;};
    CUSTOMREAL * get_lat_inv(){return lat_inv;};
    CUSTOMREAL * get_lon_inv(){return lon_inv;};

    int get_max_iter_inv()    {return max_iter_inv;};

    bool get_is_srcrec_swap() {return swap_src_rec;};

    bool get_is_output_source_field() {return is_output_source_field;};
    bool get_is_output_model_dat()    {return is_output_model_dat;};
    bool get_is_verbose_output()      {return is_verbose_output;};
    bool get_is_output_final_model()  {return is_output_final_model;};

    bool get_is_inv_slowness()        {return is_inv_slowness;};
    bool get_is_inv_azi_ani()         {return is_inv_azi_ani;};
    bool get_is_inv_rad_ani()         {return is_inv_rad_ani;};
    CUSTOMREAL * get_kernel_taper()   {return kernel_taper;};
private:
    // boundary information
    CUSTOMREAL min_dep; // minimum depth in km
    CUSTOMREAL max_dep; // maximum depth in km
    CUSTOMREAL min_lat; // minimum latitude in degree
    CUSTOMREAL max_lat; // maximum latitude in degree
    CUSTOMREAL min_lon; // minimum longitude in degree
    CUSTOMREAL max_lon; // maximum longitude in degree

    // source information
    CUSTOMREAL src_dep;              // source depth in km
    CUSTOMREAL src_lat;              // source latitude in degrees
    CUSTOMREAL src_lon;              // source longitude in degrees
    std::string src_rec_file;        // source receiver file
    std::string src_list_file;       // source list file
    std::string rec_list_file;       // source list file
    std::string src_rec_file_out;    // source receiver file to be output
    std::string sta_correction_file;         // station correction file to be input
    std::string station_correction_file_out; // station correction file to be output
    bool src_rec_file_exist = false; // source receiver file exist
    bool sta_correction_file_exist = false;   // station correction file exist

    // model input files
    std::string init_model_type; // model type
    std::string init_model_path; // model file path init
    std::string model_1d_name;   // name of 1d model for teleseismic tomography

    // inversion
    int run_mode=0;                  // do inversion or not (0: no, 1: yes)
    int n_inversion_grid=1;              // number of inversion grid
    int type_dep_inv=0, type_lat_inv=0, type_lon_inv=0; // uniform or flexible inversion grid (0: uniform, 1: flexible)
    // type = 0: uniform inversion grid
    int n_inv_r=1, n_inv_t=1, n_inv_p=1; // number of inversion grid in r, t, p direction
    // inversion grid
    CUSTOMREAL min_dep_inv=-99999; // minimum depth in km
    CUSTOMREAL max_dep_inv=-99999; // maximum depth in km
    CUSTOMREAL min_lat_inv=-99999; // minimum latitude
    CUSTOMREAL max_lat_inv=-99999; // maximum latitude
    CUSTOMREAL min_lon_inv=-99999; // minimum longitude
    CUSTOMREAL max_lon_inv=-99999; // maximum longitude
    // type = 1: flexible inversion grid
    CUSTOMREAL *dep_inv, *lat_inv, *lon_inv;   // flexibly designed inversion grid
    int n_inv_r_flex=1, n_inv_t_flex=1, n_inv_p_flex=1; // number of flexibly designed inversion grid in r, t, p direction
    bool n_inv_r_flex_read = false, n_inv_t_flex_read = false, n_inv_p_flex_read = false; // flag if n inv grid flex is read or not. if false, code allocate dummy memory

    // convergence setting
    CUSTOMREAL conv_tol;       // convergence tolerance
    int        max_iter;       // maximum number of iterations
    int        max_iter_inv=1; // maximum number of iterations for inversion

    // calculation method
    int stencil_order; // stencil order
    int stencil_type; // stencil order
    int sweep_type; // sweep type (0: legacy, 1: cuthil-mckee with shm parallelization)

    // parse src_rec_file
    void parse_src_rec_file();
    // parse sta_correction_file
    void parse_sta_correction_file();

public:
    // new version for reading src rec data (by Chen Jing, 20230212)
    std::map<std::string, SrcRecInfo> src_list_nv;
    std::map<std::string, SrcRecInfo> rec_list_nv;
public:
    // functions for SrcRecInfo
    void set_adjoint_source(std::string, CUSTOMREAL);

    // new version for reading src rec data (by Chen Jing, 20230212)
    auto get_src_list_nv_begin()                        {return src_list_nv.begin();};
    auto get_src_list_nv_end()                          {return src_list_nv.end();};
    SrcRecInfo get_src_list_nv(std::string);
    
    std::map<std::string, SrcRecInfo> src_list_back_nv;     // for output purposes
    std::map<std::string, SrcRecInfo> src_list_prepare_nv;  // related to common receiver differential time, computed fist
    std::map<std::string, SrcRecInfo> src_list_tele_nv;    // source list for teleseismic
    std::vector<std::string> src_name_list;                 // name list for output

    auto get_rec_list_nv_begin()                        {return rec_list_nv.begin();};
    auto get_rec_list_nv_end()                          {return rec_list_nv.end();};
    SrcRecInfo get_rec_list_nv(std::string);    
    std::map<std::string, SrcRecInfo> rec_list_back_nv;    // for output purposes
    std::map<std::string, SrcRecInfo> rec_list_tele_nv;    // rec list for teleseismic

    std::vector<DataInfo> data_info_nv;
    std::vector<DataInfo> tele_data_info_nv;    // data list for teleseismic
    std::vector<DataInfo> data_info_back_nv;    // for backup purposes


    // std::map<std::vector<std::string>,CUSTOMREAL> syn_time_list;    // (evname, stname) -> syn_time
    std::map<std::string, std::map<std::string, CUSTOMREAL> > syn_time_list_sr;     // all used synthetic traveltime in forward modeling and inversion.  two level map, map1: source -> map2;  map 2: receiver -> time;
    std::map<std::string, std::vector<DataInfo> > data_info_smap;      // map source -> vector; vector: (related) Traveltime data
    std::map<std::string, std::vector<DataInfo> > data_info_smap_reloc;      // map source -> vector; vector: (related) Traveltime data
    
    // initialize_adjoint_source
    void initialize_adjoint_source();

    // weight of different types of data
    CUSTOMREAL abs_time_local_weight    = 1.0;    // weight of absolute traveltime data for local earthquake,                        default: 1.0
    CUSTOMREAL cr_dif_time_local_weight = 1.0;    // weight of common receiver differential traveltime data for local earthquake,    default: 1.0
    CUSTOMREAL cs_dif_time_local_weight = 1.0;    // weight of common source differential traveltime data for local earthquake,      default: 1.0    (not ready)
    CUSTOMREAL teleseismic_weight       = 1.0;    // weight of teleseismic data                                                      default: 1.0    (not ready)
    // misfit balance
    int is_balance_data_weight = 0;    // add the weight to normalize the initial objective function of different types of data. 1 for yes and 0 for no
    CUSTOMREAL total_abs_local_data_weight = 0.0;
    CUSTOMREAL total_cr_dif_local_data_weight = 0.0;
    CUSTOMREAL total_cs_dif_local_data_weight = 0.0;
    CUSTOMREAL total_teleseismic_data_weight = 0.0;
    // the number of data
    int N_abs_local_data = 0;
    int N_cr_dif_local_data = 0;
    int N_cs_dif_local_data = 0;
    int N_teleseismic_data = 0;
    // the number of the types of data
    int N_data_type = 0;
    std::map<std::string, int> data_type;     // element: "abs", "cs_dif", "cr_dif", "tele"

    CUSTOMREAL max_change_dep = 10.0; 
    CUSTOMREAL max_change_lat = 1.0; 
    CUSTOMREAL max_change_lon = 1.0; 
    int is_ortime_local_search = 0;
    CUSTOMREAL ref_ortime_change = 5.0;
    CUSTOMREAL max_change_ortime = 0.5;
    CUSTOMREAL step_length_ortime_rescale = 0.1;
    
    // std::map<std::string,std::string> name_for_reloc;
    std::vector<std::string> name_for_reloc;
private:
    // rec_map.  rec_map: rec_name -> rec_id;  rec_map_reverse: rec_id -> rec_name
    // std::map<std::string, SrcRec> rec_list;
    // std::map<std::string, CUSTOMREAL> station_correction;
    // std::map<std::string, CUSTOMREAL> station_correction_kernel;
    // CUSTOMREAL* station_correction_kernel;
    // CUSTOMREAL* station_correction_value;
    // int N_receiver = 0;

    // list for all of src point data
    // std::vector<SrcRec> src_points;
    // std::vector<SrcRec> src_points_back; // for backup purposes
    // std::vector<SrcRec> src_points_out;  // temporal storage for output

    
    // list for all of rec point data
    // std::vector< std::vector<SrcRec> > rec_points;
    // std::vector< std::vector<SrcRec> > rec_points_back; // for backup purposes
    // std::vector< std::vector<SrcRec> > rec_points_out;  // temporal storage for output
    // gather all arrival times to a main process
    // void gather_all_arrival_times_to_main();
    void gather_all_arrival_times_to_main_nv();

    // swap the sources and receivers
    void do_swap_src_rec();
    bool swap_src_rec = false;     // whether the src/rec are swapped or not
    // void reverse_src_rec_points(); // reverse the swapped src/rec

    // rearrange the data_info_nv to data_info_smap
    void rearrange_data_info();

    void generate_src_list_prepare();
public:
    // generate syn_time_list_nv based on data
    void generate_syn_time_list();
    void initialize_syn_time_list();
    void reduce_syn_time_list();
private:
    // tele seismic source management
    void separate_region_and_tele_src();               // check if the source is tele seismic or not
    void merge_region_and_tele_src();                  // merge tele seismic source to the region source
    // std::vector<SrcRec> tele_src_points;               // tele seismic source points
    // std::vector<std::vector<SrcRec> > tele_rec_points; // tele seismic receiver points
    bool i_first=false, i_last=false, j_first=false, j_last=false, k_first=false; // store info if this subdomain has outer boundary

public:
    void allocate_memory_tele_boundaries(int, int, int, std::string,
        bool, bool, bool, bool, bool); // allocate memory for tele boundaries
private:
    // check contradictions in input parameters
    void check_contradictions();

public:
    // prepare source list for this simulation group
    void prepare_src_list();
    // src points for this sim group
    std::vector<int> src_ids_this_sim;
    std::vector<std::string> src_names_this_sim;

    // traveltime of src should be prepared for this sim group (have common receivr differential traveltime)
    std::vector<int> src_ids_this_sim_prepare;
    std::vector<std::string> src_names_this_sim_prepare;

    // boundary of src should be computed for this sim group (teleseismic earthquake)
    std::vector<int> src_ids_this_sim_tele;
    std::vector<std::string> src_names_this_sim_tele;

    // (relocation) modify (swapped source) receiver's location and time
    void modift_swapped_source_location();

    // write out src_rec_file
    // void write_src_rec_file(int);
    void write_src_rec_file_nv(int);

    // write out station_correction_file
    void write_station_correction_file(int);

    // station correction
    // void station_correction_kernel();
    void station_correction_update( CUSTOMREAL );

private:
    // output setting
    bool is_output_source_field = false; // output out_data_sim_X.h or not.
    bool is_output_model_dat    = false; // output model_parameters_inv_0000.dat or not.
    bool is_verbose_output      = false; // output verbose information or not.
    bool is_output_final_model  = true;  // output merged final model or not.

    // inversion setting
    bool is_inv_slowness = true;  // update slowness (velocity) or not.
    bool is_inv_azi_ani  = false; // update azimuthal anisotropy (xi, eta) or not.
    bool is_inv_rad_ani  = false; // update radial anisotropy (in future) or not.

    CUSTOMREAL kernel_taper[2] = {-9999999, -9999998};   // kernel weight:  0: -inf ~ taper[0]; 0 ~ 1 : taper[0] ~ taper[1]; 1 : taper[1] ~ inf
    bool is_sta_correction = false; // apply station correction or not.
};



#endif