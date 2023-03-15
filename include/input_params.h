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
#include "timer.h"
#include "src_rec.h"


// InputParams class for reading/storing and outputing input parameters and src_rec information
class InputParams {

public:
    InputParams(std::string&);
    ~InputParams();

    // write parameters to output file
    void write_params_to_file();

    // get parameters
    CUSTOMREAL get_min_dep(){return min_dep;};
    CUSTOMREAL get_max_dep(){return max_dep;};
    CUSTOMREAL get_min_lat(){return min_lat*DEG2RAD;};
    CUSTOMREAL get_max_lat(){return max_lat*DEG2RAD;};
    CUSTOMREAL get_min_lon(){return min_lon*DEG2RAD;};
    CUSTOMREAL get_max_lon(){return max_lon*DEG2RAD;};

    CUSTOMREAL                  get_src_radius(const std::string&);
    CUSTOMREAL                  get_src_lat(   const std::string&);
    CUSTOMREAL                  get_src_lon(   const std::string&);
    std::string                 get_src_rec_file()      {return src_rec_file;};
    bool                        get_src_rec_file_exist(){return src_rec_file_exist;};
    bool                        get_if_src_teleseismic(const std::string&); // return true if the source is teleseismic
    SrcRecInfo&                 get_src_point(const std::string&);          // return SrcRec object
    SrcRecInfo&                 get_rec_point(const std::string&);          // return receivers for the current source
    std::vector<std::string>    get_rec_names(const std::string&);          // return SrcRec object

    CUSTOMREAL get_conv_tol()                    {return conv_tol;};
    void       set_conv_tol(CUSTOMREAL conv_tol_){conv_tol = conv_tol_;};
    CUSTOMREAL get_max_iter()                    {return max_iter;};

    int  get_stencil_order()                  {return stencil_order;};
    void set_stencil_order(int stencil_order_){stencil_order = stencil_order_;};
    int  get_stencil_type()                   {return stencil_type;};
    int  get_sweep_type()                     {return sweep_type;};

    std::string get_init_model_path(){return init_model_path;};
    std::string get_model_1d_name()  {return model_1d_name;};

    int get_run_mode()        {return run_mode;};
    int get_n_inversion_grid(){return n_inversion_grid;};

    int get_type_dep_inv(){return type_dep_inv;};
    int get_type_lat_inv(){return type_lat_inv;};
    int get_type_lon_inv(){return type_lon_inv;};

    // type = 0:
    int        get_n_inv_r()    {return n_inv_r;};
    int        get_n_inv_t()    {return n_inv_t;};
    int        get_n_inv_p()    {return n_inv_p;};
    CUSTOMREAL get_min_dep_inv(){return min_dep_inv;};
    CUSTOMREAL get_max_dep_inv(){return max_dep_inv;};
    CUSTOMREAL get_min_lat_inv(){return min_lat_inv*DEG2RAD;};
    CUSTOMREAL get_max_lat_inv(){return max_lat_inv*DEG2RAD;};
    CUSTOMREAL get_min_lon_inv(){return min_lon_inv*DEG2RAD;};
    CUSTOMREAL get_max_lon_inv(){return max_lon_inv*DEG2RAD;};

    // type = 1:
    int         get_n_inv_r_flex(){return n_inv_r_flex;};
    int         get_n_inv_t_flex(){return n_inv_t_flex;};
    int         get_n_inv_p_flex(){return n_inv_p_flex;};
    CUSTOMREAL* get_dep_inv()     {return dep_inv;};
    CUSTOMREAL* get_lat_inv()     {return lat_inv;};
    CUSTOMREAL* get_lon_inv()     {return lon_inv;};

    int  get_max_iter_inv() {return max_iter_inv;};

    bool get_is_srcrec_swap() {return swap_src_rec;};

    bool get_is_output_source_field() {return is_output_source_field;};
    bool get_is_output_model_dat()    {return is_output_model_dat;};
    bool get_is_verbose_output()      {return is_verbose_output;};
    bool get_is_output_final_model()  {return is_output_final_model;};
    bool get_is_output_in_process()   {return is_output_in_process;};

    bool get_is_inv_slowness()        {return is_inv_slowness;};
    bool get_is_inv_azi_ani()         {return is_inv_azi_ani;};
    bool get_is_inv_rad_ani()         {return is_inv_rad_ani;};
    CUSTOMREAL * get_kernel_taper()   {return kernel_taper;};

    bool get_is_single_precision_output() {return is_single_precision_output;};

    // prepare source list for this simulation group
    void prepare_src_map();

    // (relocation) modify (swapped source) receiver's location and time
    void modift_swapped_source_location();

    // write out src_rec_file
    void write_src_rec_file(int);

    // write out station_correction_file
    void write_station_correction_file(int);

    // station correction
    void station_correction_update( CUSTOMREAL );

    std::map<std::string, SrcRecInfo> src_map_all;          // map of all sources (full information is only stored by the main process)
    std::map<std::string, SrcRecInfo> src_map;              // map of sources belonging to this simultaneous group
    std::map<std::string, SrcRecInfo> src_map_comm_src;     // map of sources with common source
    std::map<std::string, SrcRecInfo> src_map_tele;         // source list for teleseismic

    std::map<std::string, SrcRecInfo> rec_map_all;     // map of all receivers (full information is only stored by the main process)
    std::map<std::string, SrcRecInfo> rec_map;         // map of receivers belonging to this simultaneous group
    std::map<std::string, SrcRecInfo> rec_map_tele;    // rec list for teleseismic

    // datainfo maps with key <id_src, rec_name>
    std::map< std::string, std::map<std::string, std::vector<DataInfo>>> data_map_all;     // data list for all data (full information is only stored by the main process)
    std::map< std::string, std::map<std::string, std::vector<DataInfo>>> data_map;         // data list for this simultaneous group
    std::map< std::string, std::map<std::string, std::vector<DataInfo>>> data_map_tele;    // data list for teleseismic

    std::vector<std::string> name_for_reloc;    // name list of receivers (swarpped sources) for location

    // src id <-> src name relations
    std::vector<std::string>          src_id2name;          // name list of sources belonging to this simultaneous group
    std::vector<std::string>          src_id2name_comm_src; // name list of sources with common source
    std::vector<std::string>          src_id2name_all;      // name list of all sources (store the order of sources in src_rec file)
    std::vector<std::string>          src_id2name_back;     // back up of name list of all sources (this will not be swapped)

    // backups used when outputing the data
    std::map<std::string, SrcRecInfo>                                   src_map_back;
    std::map<std::string, SrcRecInfo>                                   rec_map_back;
    std::map<std::string, std::map<std::string, std::vector<DataInfo>>> data_map_back;



    // the number of data
    int N_abs_local_data    = 0;
    int N_cr_dif_local_data = 0;
    int N_cs_dif_local_data = 0;
    int N_teleseismic_data  = 0;
    // the number of the types of data
    int N_data_type = 0;
    std::map<std::string, int> data_type;     // element: "abs", "cs_dif", "cr_dif", "tele"

    // initialize_adjoint_source
    void initialize_adjoint_source();
    // set adjoint source
    void set_adjoint_source(std::string, CUSTOMREAL);

    void allocate_memory_tele_boundaries(int, int, int, std::string,
        bool, bool, bool, bool, bool); // allocate memory for tele boundaries


private:
    // boundary information
    CUSTOMREAL min_dep; // minimum depth in km
    CUSTOMREAL max_dep; // maximum depth in km
    CUSTOMREAL min_lat; // minimum latitude in degree
    CUSTOMREAL max_lat; // maximum latitude in degree
    CUSTOMREAL min_lon; // minimum longitude in degree
    CUSTOMREAL max_lon; // maximum longitude in degree

    // source information
    CUSTOMREAL  src_dep;                     // source depth in km
    CUSTOMREAL  src_lat;                     // source latitude in degrees
    CUSTOMREAL  src_lon;                     // source longitude in degrees
    std::string src_rec_file;                // source receiver file
    std::string src_rec_file_out;            // source receiver file to be output
    std::string sta_correction_file;         // station correction file to be input
    std::string station_correction_file_out; // station correction file to be output
    bool src_rec_file_exist        = false;  // source receiver file exist
    bool sta_correction_file_exist = false;  // station correction file exist
    bool swap_src_rec              = false;  // whether the src/rec are swapped or not

    // model input files
    std::string init_model_path; // model file path init
    std::string model_1d_name;   // name of 1d model for teleseismic tomography

    // inversion
    int run_mode=0;                                     // do inversion or not (0: no, 1: yes)
    int n_inversion_grid=1;                             // number of inversion grid
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
    CUSTOMREAL *dep_inv, *lat_inv, *lon_inv;            // flexibly designed inversion grid
    int n_inv_r_flex=1, n_inv_t_flex=1, n_inv_p_flex=1; // number of flexibly designed inversion grid in r, t, p direction
    bool n_inv_r_flex_read = false, n_inv_t_flex_read = false, n_inv_p_flex_read = false; // flag if n inv grid flex is read or not. if false, code allocate dummy memory

    // convergence setting
    CUSTOMREAL conv_tol;       // convergence tolerance
    int        max_iter;       // maximum number of iterations
    int        max_iter_inv=1; // maximum number of iterations for inversion

    // calculation method
    int stencil_order = 3; // stencil order
    int stencil_type  = 0; // stencil type  (0: non upwind; 1: upwind)
    int sweep_type    = 0; // sweep type (0: legacy, 1: cuthil-mckee with shm parallelization)

    // gather all arrival times to a main process
    void gather_all_arrival_times_to_main();
    // geneerate a map of sources which include common source double difference data
    void generate_src_map_with_common_source();

    bool i_first=false, i_last=false, \
         j_first=false, j_last=false, \
         k_first=false; // store info if this subdomain has outer boundary

    // check contradictions in input parameters
    void check_contradictions();

    // output setting
    bool is_output_source_field = false; // output out_data_sim_X.h or not.
    bool is_output_model_dat    = false; // output model_parameters_inv_0000.dat or not.
    bool is_verbose_output      = false; // output verbose information or not.
    bool is_output_final_model  = true;  // output merged final model or not.
    bool is_output_in_process   = true;  // output merged model at each inv iteration or not.

    // inversion setting
    bool is_inv_slowness = true;  // update slowness (velocity) or not.
    bool is_inv_azi_ani  = true; // update azimuthal anisotropy (xi, eta) or not.
    bool is_inv_rad_ani  = false; // update radial anisotropy (in future) or not.

    CUSTOMREAL kernel_taper[2] = {-9999999, -9999998};   // kernel weight:  0: -inf ~ taper[0]; 0 ~ 1 : taper[0] ~ taper[1]; 1 : taper[1] ~ inf
    bool is_sta_correction = false; // apply station correction or not.

    // single precision (float) output mode
    bool is_single_precision_output = false;

};


// some utils

inline DataInfo& get_data_rec_pair(std::vector<DataInfo>& v){
    // return the first element in the vector with is_rec_pair = true
    for (auto it = v.begin(); it != v.end(); it++){
        if (it->is_rec_pair)
            return *it;
    }

    // error if no rec pair is found
    std::cerr << "Error: no rec pair is found in get_data_rec_pair" << std::endl;
    exit(1);

    // return the first element in the vector as a dummy
    return v[0];
}


#endif