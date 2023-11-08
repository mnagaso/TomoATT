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
#include <stdexcept>

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

    //
    // getter/setter
    //

    //
    // bondary information
    //
    CUSTOMREAL get_min_dep(){return min_dep;};
    CUSTOMREAL get_max_dep(){return max_dep;};
    CUSTOMREAL get_min_lat(){return min_lat*DEG2RAD;};
    CUSTOMREAL get_max_lat(){return max_lat*DEG2RAD;};
    CUSTOMREAL get_min_lon(){return min_lon*DEG2RAD;};
    CUSTOMREAL get_max_lon(){return max_lon*DEG2RAD;};

    //
    // source receiver data
    //

    // source receiver information for processes which stores source receiver information
    std::string                 get_src_rec_file()      {return src_rec_file;};
    bool                        get_src_rec_file_exist(){return src_rec_file_exist;};
    SrcRecInfo&                 get_src_point(const std::string&);          // return SrcRec object
    SrcRecInfo&                 get_rec_point(const std::string&);          // return receivers for the current source

    // source receiver information with broadcast to all subdom_main processes
    SrcRecInfo                  get_src_point_bcast(const std::string&);    // return SrcRec object
    SrcRecInfo                  get_rec_point_bcast(const std::string&);    // return receivers for the current source
    CUSTOMREAL                  get_src_radius(   const std::string&);
    CUSTOMREAL                  get_src_lat(      const std::string&);
    CUSTOMREAL                  get_src_lon(      const std::string&);
    SrcRecInfo                  get_src_point_bcast_2d(const std::string&);    // return SrcRec object
    CUSTOMREAL                  get_src_radius_2d(const std::string&);
    CUSTOMREAL                  get_src_lat_2d(   const std::string&);
    CUSTOMREAL                  get_src_lon_2d(   const std::string&);
    std::string                 get_src_name(const int&);                   // return source name from in-sim_group id
    std::string                 get_rec_name(const int&);                   // return receiver name from in-sim_group id
    std::string                 get_src_name_comm(const int&);              // return source name in common receiver list
    int                         get_src_id(const std::string&);             // return src global id from src name
    bool                        get_if_src_teleseismic(const std::string&); // return true if the source is teleseismic

    //
    // others
    //
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

    int get_type_invgrid_dep(){return type_invgrid_dep;};
    int get_type_invgrid_lat(){return type_invgrid_lat;};
    int get_type_invgrid_lon(){return type_invgrid_lon;};
    int get_type_invgrid_dep_ani(){return type_invgrid_dep_ani;};
    int get_type_invgrid_lat_ani(){return type_invgrid_lat_ani;};
    int get_type_invgrid_lon_ani(){return type_invgrid_lon_ani;};


    // type = 0:
    int get_n_inv_r()     {return n_inv_r;};
    int get_n_inv_t()     {return n_inv_t;};
    int get_n_inv_p()     {return n_inv_p;};
    int get_n_inv_r_ani() {return n_inv_r_ani;};
    int get_n_inv_t_ani() {return n_inv_t_ani;};
    int get_n_inv_p_ani() {return n_inv_p_ani;};

    CUSTOMREAL get_min_dep_inv(){return min_dep_inv;};
    CUSTOMREAL get_max_dep_inv(){return max_dep_inv;};
    CUSTOMREAL get_min_lat_inv(){return min_lat_inv*DEG2RAD;};
    CUSTOMREAL get_max_lat_inv(){return max_lat_inv*DEG2RAD;};
    CUSTOMREAL get_min_lon_inv(){return min_lon_inv*DEG2RAD;};
    CUSTOMREAL get_max_lon_inv(){return max_lon_inv*DEG2RAD;};
    CUSTOMREAL get_min_dep_inv_ani(){return min_dep_inv_ani;};
    CUSTOMREAL get_max_dep_inv_ani(){return max_dep_inv_ani;};
    CUSTOMREAL get_min_lat_inv_ani(){return min_lat_inv_ani*DEG2RAD;};
    CUSTOMREAL get_max_lat_inv_ani(){return max_lat_inv_ani*DEG2RAD;};
    CUSTOMREAL get_min_lon_inv_ani(){return min_lon_inv_ani*DEG2RAD;};
    CUSTOMREAL get_max_lon_inv_ani(){return max_lon_inv_ani*DEG2RAD;};

    // type = 1:
    int         get_n_inv_r_flex(){return n_inv_r_flex;};
    int         get_n_inv_t_flex(){return n_inv_t_flex;};
    int         get_n_inv_p_flex(){return n_inv_p_flex;};
    CUSTOMREAL* get_dep_inv()     {return dep_inv;};
    CUSTOMREAL* get_lat_inv()     {return lat_inv;};
    CUSTOMREAL* get_lon_inv()     {return lon_inv;};

    int         get_n_inv_r_flex_ani(){return n_inv_r_flex_ani;};
    int         get_n_inv_t_flex_ani(){return n_inv_t_flex_ani;};
    int         get_n_inv_p_flex_ani(){return n_inv_p_flex_ani;};
    CUSTOMREAL* get_dep_inv_ani()     {return dep_inv_ani;};
    CUSTOMREAL* get_lat_inv_ani()     {return lat_inv_ani;};
    CUSTOMREAL* get_lon_inv_ani()     {return lon_inv_ani;};

    // type = 2:
    int         get_n_inv_t_trape()     {return n_inv_t_trape;};
    int         get_n_inv_p_trape()     {return n_inv_p_trape;};
    CUSTOMREAL* get_lat_spacing_inv()   {return lat_spacing_inv;};
    CUSTOMREAL* get_lon_spacing_inv()   {return lon_spacing_inv;};

    int         get_n_inv_t_trape_ani()     {return n_inv_t_trape_ani;};
    int         get_n_inv_p_trape_ani()     {return n_inv_p_trape_ani;};
    CUSTOMREAL* get_lat_spacing_inv_ani()   {return lat_spacing_inv_ani;};
    CUSTOMREAL* get_lon_spacing_inv_ani()   {return lon_spacing_inv_ani;};

    // invgrid for ani
    bool        get_invgrid_ani()     {return invgrid_ani;};

    bool get_invgrid_volume_rescale(){return invgrid_volume_rescale;}

    int  get_max_iter_inv() {return max_iter_inv;};

    int  get_model_update_N_iter(){return model_update_N_iter;};
    int  get_relocation_N_iter()  {return relocation_N_iter;};
    int  get_max_loop()           {return max_loop;};

    bool get_is_srcrec_swap() {return swap_src_rec;};

    bool get_if_output_source_field()    {return output_source_field;};
    bool get_if_output_model_dat()       {return output_model_dat;};
    bool get_if_output_final_model()     {return output_final_model;};
    bool get_if_output_in_process()      {return output_in_process;};
    bool get_if_output_in_process_data() {return output_in_process_data;};
    bool get_if_single_precision_output(){return single_precision_output;};
    int  get_verbose_output_level()      {return verbose_output_level;}; // #TODO: modify codes for verbose_output_level > 1

    bool get_update_slowness()        {return update_slowness;};
    bool get_update_azi_ani()         {return update_azi_ani;};
    bool get_update_rad_ani()         {return update_rad_ani;};
    bool get_update_qp()              {return update_qp;};
    CUSTOMREAL * get_depth_taper()    {return depth_taper;};

    bool get_use_abs()                  {return use_abs;};
    bool get_use_cs()                   {return use_cs;};
    bool get_use_cr()                   {return use_cr;};

    bool get_use_abs_reloc()                  {return use_abs_reloc;};
    bool get_use_cr_reloc()                   {return use_cr_reloc;};

    CUSTOMREAL* get_residual_weight_abs()   {return residual_weight_abs;};
    CUSTOMREAL* get_distance_weight_abs()   {return distance_weight_abs;};
    CUSTOMREAL* get_residual_weight_cr()    {return residual_weight_cr;};
    CUSTOMREAL* get_azimuthal_weight_cr()   {return azimuthal_weight_cr;};
    CUSTOMREAL* get_residual_weight_cs()    {return residual_weight_cs;};
    CUSTOMREAL* get_azimuthal_weight_cs()   {return azimuthal_weight_cs;};

    CUSTOMREAL* get_residual_weight_abs_reloc()   {return residual_weight_abs_reloc;};
    CUSTOMREAL* get_distance_weight_abs_reloc()   {return distance_weight_abs_reloc;};
    CUSTOMREAL* get_residual_weight_cr_reloc()    {return residual_weight_cr_reloc;};
    CUSTOMREAL* get_azimuthal_weight_cr_reloc()   {return azimuthal_weight_cr_reloc;};

    // get if the T field is written into the file
    bool get_is_T_written_into_file(const std::string&);

    // prepare source list for this simulation group
    void prepare_src_map();

    // (relocation) modify (swapped source) receiver's location and time
    void modify_swapped_source_location();

    // write out src_rec_file
    void write_src_rec_file(int,int);

    // write out station_correction_file
    void write_station_correction_file(int);

    // station correction
    void station_correction_update( CUSTOMREAL );

    int n_src_this_sim_group = 0;          // number of sources in this simultaneous group
    int n_src_comm_rec_this_sim_group = 0; // number of sources with common receiver in this simultaneous group
    int n_src_2d_this_sim_group = 0;       // number of sources assigned for 2d solver in this simultaneous group
    int n_rec_this_sim_group = 0;          // number of receivers in this simultaneous group

    std::map<std::string, SrcRecInfo> src_map_all;          // map of all sources (full information is only stored by the main process)
    std::map<std::string, SrcRecInfo> src_map;              // map of sources belonging to this simultaneous group
    std::map<std::string, SrcRecInfo> src_map_comm_rec;     // map of sources with common receiver
    std::map<std::string, SrcRecInfo> src_map_2d;           // map of sources assigned for 2d solver
    std::map<std::string, SrcRecInfo> src_map_tele;         // source list for teleseismic

    std::map<std::string, SrcRecInfo> rec_map_all;     // map of all receivers (full information is only stored by the main process)
    std::map<std::string, SrcRecInfo> rec_map;         // map of receivers belonging to this simultaneous group
    std::map<std::string, SrcRecInfo> rec_map_tele;    // rec list for teleseismic

    // datainfo-vector maps <src_name, rec_name>
    std::map< std::string, std::map<std::string, std::vector<DataInfo>>> data_map_all;     // data list for all data (full information is only stored by the main process)
    std::map< std::string, std::map<std::string, std::vector<DataInfo>>> data_map;         // data list for this simultaneous group
    std::map< std::string, std::map<std::string, std::vector<DataInfo>>> data_map_tele;    // data list for teleseismic

    std::vector<std::string> name_for_reloc;    // name list of receivers (swarpped sources) for location

    // src id <-> src name relations
    std::vector<std::string>                           src_id2name;          // name list of sources belonging to this simultaneous group
    std::vector<std::string>                           rec_id2name;          // name list of receivers belongig to this simultaneous group
    std::vector<std::string>                           src_id2name_comm_rec; // name list of sources with common receiver
    std::vector<std::string>                           src_id2name_2d;       // name list of sources assigned for 2d solver.
    std::vector<std::string>                           src_id2name_all;      // name list of all sources (store the order of sources in src_rec file)
    std::vector<std::string>                           src_id2name_back;     // back up of name list of all sources (this will not be swapped)
    std::vector<std::vector<std::vector<std::string>>> rec_id2name_back;     // back up of the name list of all receivers for each source (this will not be swapped)

    // backups used when outputing the data
    std::map<std::string, SrcRecInfo>                                   src_map_back;
    std::map<std::string, SrcRecInfo>                                   rec_map_back;
    std::map<std::string, std::map<std::string, std::vector<DataInfo>>> data_map_back;

    // the number of data
    int N_abs_local_data    = 0;
    int N_cr_dif_local_data = 0;
    int N_cs_dif_local_data = 0;
    int N_teleseismic_data  = 0;
    int N_data              = 0;
    // the number of the types of data
    int N_data_type = 0;
    std::map<std::string, int> data_type;     // element: "abs", "cs_dif", "cr_dif", "tele"

    // initialize_adjoint_source
    void initialize_adjoint_source();
    // set adjoint source
    void set_adjoint_source(std::string, CUSTOMREAL);

    // gather traveltimes and calculate differences of synthetic data
    void gather_traveltimes_and_calc_syn_diff();

    // reduce necessary data in rec_map, which has differed elements in each sim group.
    template <typename T>
    void allreduce_rec_map_var(T&);

    void allreduce_rec_map_vobj_src_reloc();
    void allreduce_rec_map_grad_src();

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
    int run_mode=0;                                                 // do inversion or not (0: no, 1: yes)
    int n_inversion_grid=1;                                         // number of inversion grid
    // uniform or flexible inversion grid (0: uniform, 1: flexible)
    int type_invgrid_dep=0;
    int type_invgrid_lat=0;
    int type_invgrid_lon=0;
    // uniform or flexible anisotropic inversion grid (0: uniform, 1: flexible)
    int type_invgrid_dep_ani=0;
    int type_invgrid_lat_ani=0;
    int type_invgrid_lon_ani=0;

    //
    // variables for type = 0: uniform inversion grid
    //
    // number of uniform inversion grid nodes in r, t, p direction
    int n_inv_r=1;
    int n_inv_t=1;
    int n_inv_p=1;
    // number of uniform anisotropic inversion grid nodes in r, t, p direction
    int n_inv_r_ani=1;
    int n_inv_t_ani=1;
    int n_inv_p_ani=1;

    // min max values of inversion grid
    CUSTOMREAL min_dep_inv    =-99999; // minimum depth in km
    CUSTOMREAL max_dep_inv    =-99999; // maximum depth in km
    CUSTOMREAL min_lat_inv    =-99999; // minimum latitude
    CUSTOMREAL max_lat_inv    =-99999; // maximum latitude
    CUSTOMREAL min_lon_inv    =-99999; // minimum longitude
    CUSTOMREAL max_lon_inv    =-99999; // maximum longitude
    // min max values of anisotropic inversion grid
    CUSTOMREAL min_dep_inv_ani=-99999; // minimum depth in km
    CUSTOMREAL max_dep_inv_ani=-99999; // maximum depth in km
    CUSTOMREAL min_lat_inv_ani=-99999; // minimum latitude
    CUSTOMREAL max_lat_inv_ani=-99999; // maximum latitude
    CUSTOMREAL min_lon_inv_ani=-99999; // minimum longitude
    CUSTOMREAL max_lon_inv_ani=-99999; // maximum longitude

    //
    // variables fo type = 1: flexible inversion grid
    //
    CUSTOMREAL *dep_inv; // array for storing inversion grid points in depth
    CUSTOMREAL *lat_inv; // array for storing inversion grid points in latitude
    CUSTOMREAL *lon_inv; // array for storing inversion grid points in longitude
    CUSTOMREAL *dep_inv_ani; // array for storing inversion grid points in depth for anisotropy
    CUSTOMREAL *lat_inv_ani; // array for storing inversion grid points in latitude for anisotropy
    CUSTOMREAL *lon_inv_ani; // array for storing inversion grid points in longitude for anisotropy

    // number of flexibly designed inversion grid in r, t, p direction
    int n_inv_r_flex=1;
    int n_inv_t_flex=1;
    int n_inv_p_flex=1;
    // number of flexibly designed inversion grid in r, t, p direction
    int n_inv_r_flex_ani=1;
    int n_inv_t_flex_ani=1;
    int n_inv_p_flex_ani=1;

    // flag if n inv grid flex is read or not. if false, code allocate dummy memory
    bool n_inv_r_flex_read = false;
    bool n_inv_t_flex_read = false;
    bool n_inv_p_flex_read = false;
    // flag if n inv grid flex is read or not. if false, code allocate dummy memory
    bool n_inv_r_flex_ani_read = false;
    bool n_inv_t_flex_ani_read = false;
    bool n_inv_p_flex_ani_read = false;

    // if true, use defined inversion grid for anisotropy. Otherwise, use the same inversion grid of velocity for anisotropy
    bool invgrid_ani = false;

    //
    // variables fo type = 2: trapezoid inversion grid
    //
    int n_lat_lon_spacing_inv_trape=1; 
    int n_lat_lon_spacing_inv_trape_ani=1;
    CUSTOMREAL *lat_spacing_inv; // array for storing the spacing of inversion grid points in latitude
    CUSTOMREAL *lon_spacing_inv; // array for storing the spacing of inversion grid points in longitude
    CUSTOMREAL *lat_spacing_inv_ani; // array for storing the spacing of inversion grid points in latitude
    CUSTOMREAL *lon_spacing_inv_ani; // array for storing the spacing of inversion grid points in longitude

    // number of trapezoid designed inversion grid in t, p direction
    int n_inv_t_trape=1; 
    int n_inv_p_trape=1;
    // number of trapezoid designed inversion grid in t, p direction
    int n_inv_t_trape_ani=1; 
    int n_inv_p_trape_ani=1;

    // flag if n inv grid trapezoid is read or not. if false, code allocate dummy memory
    bool n_inv_t_trape_read = false;
    bool n_inv_p_trape_read = false;
    // flag if n inv grid trapezoid is read or not. if false, code allocate dummy memory
    bool n_inv_t_trape_ani_read = false;
    bool n_inv_p_trape_ani_read = false;

    // inversion grid volume rescale (kernel -> kernel / volume of inversion grid mesh)
    bool invgrid_volume_rescale = false;

    // date usage setting and weights
    bool use_abs = false; // use absolute travel time or not
    bool use_cs  = false; // use common source double difference or not
    bool use_cr  = false; // use common receiver double difference or not
    static const int n_weight = 4;
    CUSTOMREAL residual_weight_abs[n_weight];
    CUSTOMREAL residual_weight_cs[n_weight];
    CUSTOMREAL residual_weight_cr[n_weight];
    CUSTOMREAL distance_weight_abs[n_weight];
    CUSTOMREAL azimuthal_weight_cs[n_weight];
    CUSTOMREAL azimuthal_weight_cr[n_weight];

    // for relocation
    bool use_abs_reloc = false; // use absolute travel time or not
    bool use_cr_reloc  = false; // use common source double difference or not
    CUSTOMREAL residual_weight_abs_reloc[n_weight];
    CUSTOMREAL residual_weight_cr_reloc[n_weight];
    CUSTOMREAL distance_weight_abs_reloc[n_weight];
    CUSTOMREAL azimuthal_weight_cr_reloc[n_weight];

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
    // gather rec info to main process
    void gather_rec_info_to_main();

    // generate a map of sources which include common receiver double difference data
    void generate_src_map_with_common_receiver(std::map<std::string, std::map<std::string, std::vector<DataInfo>>>&,
                                             std::map<std::string, SrcRecInfo>&,
                                             std::vector<std::string>&);

    bool i_first=false, i_last=false, \
         j_first=false, j_last=false, \
         k_first=false; // store info if this subdomain has outer boundary

    // check contradictions in input parameters
    void check_contradictions();

    // output setting
    bool output_source_field    = false; // output out_data_sim_X.h or not.
    bool output_model_dat       = false; // output model_parameters_inv_0000.dat or not.
    bool output_final_model     = true;  // output merged final model or not.
    bool output_in_process      = true;  // output merged model at each inv iteration or not.
    bool output_in_process_data = true;  // output src_rec_file at each inv iteration or not.
    int  verbose_output_level   = 0;     // output verbose information or not.

    // inversion setting
    bool update_slowness = true;  // update slowness (velocity) or not.
    bool update_azi_ani  = false; // update azimuthal anisotropy (xi, eta) or not.
    bool update_rad_ani  = false; // update radial anisotropy (in future) or not.
    bool update_qp       = false; // update Qp or not.

    CUSTOMREAL depth_taper[2] = {-9999999, -9999998};   // kernel weight:  0: -inf ~ taper[0]; 0 ~ 1 : taper[0] ~ taper[1]; 1 : taper[1] ~ inf
    bool use_sta_correction = false; // apply station correction or not.

    // single precision (float) output mode
    bool single_precision_output = false;

};


//
// utils
//
inline DataInfo& get_data_src_rec(std::vector<DataInfo>& v){
    // return the first element in the vector with is_rec_pair = true
    for (auto it = v.begin(); it != v.end(); it++){
        if (it->is_src_rec)
            return *it;
    }

    // error if no rec pair is found
    std::cerr << "Error: no src/rec is found in get_data_src_rec" << std::endl;
    exit(1);

    // return the first element in the vector as a dummy
    return v[0];
}


inline DataInfo& get_data_rec_pair(std::map<std::string, std::map<std::string, std::vector<DataInfo>>>& v,
                                   const std::string& name_src,
                                   const std::string& name_rec1,
                                   const std::string& name_rec2){
    // return the first element in the vector with is_rec_pair = true
    auto& map1 = v[name_src][name_rec1];
    auto& map2 = v[name_src][name_rec2];

    for (auto it = map1.begin(); it != map1.end(); it++){
        if (it->is_rec_pair && it->name_rec_pair[0] == name_rec1 && it->name_rec_pair[1] == name_rec2)
            return *it;
    }

    for (auto it = map2.begin(); it != map2.end(); it++){
        if (it->is_rec_pair && it->name_rec_pair[0] == name_rec1 && it->name_rec_pair[1] == name_rec2)
            return *it;
    }

    // error if no rec pair is found
    std::cerr << "Error: no rec pair is found in get_data_rec_pair" << std::endl;
    exit(1);

    // return the first element in the vector as a dummy
    return v[name_src][name_rec1][0];
}


inline DataInfo& get_data_src_pair(std::map<std::string, std::map<std::string, std::vector<DataInfo>>>& v,
                                   const std::string& name_src1,
                                   const std::string& name_src2,
                                   const std::string& name_rec){

    // return the first element in the vector with is_rec_pair = true
    auto& map1 = v[name_src1][name_rec];
    auto& map2 = v[name_src2][name_rec];

    for (auto it = map1.begin(); it != map1.end(); it++){
        if (it->is_src_pair && it->name_src_pair[0] == name_src1 && it->name_src_pair[1] == name_src2)
            return *it;
    }

    for (auto it = map2.begin(); it != map2.end(); it++){
        if (it->is_src_pair && it->name_src_pair[0] == name_src2 && it->name_src_pair[1] == name_src1)
            return *it;
    }

    // error if no src pair is found
    std::cerr << "Error: no src pair is found in get_data_src_pair" << std::endl;
    exit(1);

    // return the first element in the vector as a dummy
    return v[name_src1][name_rec][0];
}


inline void set_cr_dif_to_src_pair(std::map<std::string, std::map< std::string, std::vector<DataInfo>>>& v,
                                   std::string& name_src1,
                                   std::string& name_src2,
                                   std::string& name_rec,
                                   CUSTOMREAL& cr_dif){
    // return the first element in the vector with is_rec_pair = true

    std::vector<DataInfo>& vdata = v[name_src1][name_rec];

    for (auto it = vdata.begin(); it != vdata.end(); it++){
        if (it->is_src_pair
        && ( (it->name_src_pair[0] == name_src1 && it->name_src_pair[1] == name_src2)
         ||  (it->name_src_pair[0] == name_src2 && it->name_src_pair[1] == name_src1) )) {
            it->cr_dif_travel_time = cr_dif;
        }
    }
}


inline bool get_if_any_src_pair(std::vector<DataInfo>& v){
    // return the first element in the vector with is_rec_pair = true
    for (auto it = v.begin(); it != v.end(); it++){
        if (it->is_src_pair)
            return true;
    }

    // or maybe this will be enough
    //return v.size() > 1;

    // return false if no rec pair is found
    return false;
}


#endif