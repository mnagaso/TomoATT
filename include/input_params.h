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

    CUSTOMREAL           get_src_radius();
    CUSTOMREAL           get_src_lat();
    CUSTOMREAL           get_src_lon();
    std::string          get_src_rec_file()      {return src_rec_file;};
    bool                 get_src_rec_file_exist(){return src_rec_file_exist;};
    SrcRec&              get_src_point(int);  // return SrcRec object
    std::vector<SrcRec>& get_rec_points(int); // return receivers for the current source
    bool get_if_src_teleseismic(int); // return true if the source is teleseismic

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
    bool get_is_output_in_process()   {return is_output_in_process;};

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
    CUSTOMREAL *dep_inv, *lat_inv, *lon_inv;   // flexibly designed inversion grid
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

    // rec_map.  rec_map: rec_name -> rec_id;  rec_map_reverse: rec_id -> rec_name
    std::map<std::string, SrcRec> rec_list;
    std::map<std::string, CUSTOMREAL> station_correction;
    std::map<std::string, CUSTOMREAL> station_correction_kernel;
    CUSTOMREAL* station_correction_value;
    int N_receiver = 0;

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

    bool swap_src_rec = false;     // whether the src/rec are swapped or not

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
    // prepare source list for a forward/inverse simulation group
    void prepare_src_list();
    // prepare source list for backward simulation (swap src/rec including teleseismic src)
    void prepare_src_list_for_backward_run();
    // src points for this sim group
    std::vector<int> src_ids_this_sim;
    // write out src_rec_file
    void write_src_rec_file(int);

    // write out station_correction_file
    void write_station_correction_file(int);

    // station correction
    // void station_correction_kernel();
    void station_correction_update( CUSTOMREAL );

private:
    // broadcast source list to all processes
    void broadcast_src_list();
    // free memory for teleseismic source data
    void free_memory_tele_src();

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

    // function to communicate src_rec info
    void send_src_inter_sim(SrcRec& src, int dest){
        // send source parameters of SrcRec class to dest
        // comm is inter_sub_comm or inter_comm

        // id_src id_rec n_rec n_data
        // lon lat dep weight is_teleseismic

        send_i_single_sim(&src.id_src, dest);
        send_i_single_sim(&src.id_rec, dest);
        send_i_single_sim(&src.n_rec, dest);
        send_i_single_sim(&src.n_data, dest);
        send_cr_single_sim(&src.lon, dest);
        send_cr_single_sim(&src.lat, dest);
        send_cr_single_sim(&src.dep, dest);
        send_cr_single_sim(&src.weight, dest);
        send_bool_single_sim(&src.is_teleseismic, dest);

        send_i_single_sim(&src.n_rec_pair, dest);

    }

    void send_rec_inter_sim(std::vector<SrcRec>& vrec, int dest){

        // send receiver parameters of SrcRec class to dest
        for (auto& rec: vrec) {
            // id_src id_rec n_rec n_data
            // lon lat dep weight arr_time arr_time_ori t_adj
            // is_rec_pair n_rec_pair
            send_i_single_sim(&rec.id_src, dest);
            send_i_single_sim(&rec.id_rec, dest);
            send_i_single_sim(&rec.n_rec, dest);
            send_i_single_sim(&rec.n_data, dest);
            send_cr_single_sim(&rec.lon, dest);
            send_cr_single_sim(&rec.lat, dest);
            send_cr_single_sim(&rec.dep, dest);
            send_cr_single_sim(&rec.weight, dest);
            send_cr_single_sim(&rec.arr_time, dest);
            send_cr_single_sim(&rec.arr_time_ori, dest);
            send_cr_single_sim(&rec.t_adj, dest);
            send_bool_single_sim(&rec.is_rec_pair, dest);

            // if (is_rec_pair) id_rec_pair dep_pair lat_pair lon_pair ddt_adj_pair
            if (rec.is_rec_pair){
                // TODO : check how to send double difference pairs
                for (int i = 0; i < 2; i++){
                    send_i_single_sim(&rec.id_rec_pair[i], dest);
                    send_cr_single_sim(&rec.dep_pair[i], dest);
                    send_cr_single_sim(&rec.lat_pair[i], dest);
                    send_cr_single_sim(&rec.lon_pair[i], dest);
                    send_cr_single_sim(&rec.ddt_adj_pair[i], dest);
                    send_str_sim(rec.name_rec_pair[i], dest);
                    send_cr_single_sim(&rec.station_correction_pair[i], dest);
                }

                send_cr_single_sim(&rec.dif_arr_time, dest);
                send_cr_single_sim(&rec.dif_arr_time_ori, dest);
           }
        }

    }

    void recv_src_inter_sim(SrcRec& src, int orig){
        // recv source parameters of SrcRec class from orig
        // comm is inter_sub_comm or inter_comm

        // id_src id_rec n_rec n_data
        // lon lat dep weight is_teleseismic

        recv_i_single_sim(&src.id_src, orig);
        recv_i_single_sim(&src.id_rec, orig);
        recv_i_single_sim(&src.n_rec, orig);
        recv_i_single_sim(&src.n_data, orig);
        recv_cr_single_sim(&src.lon, orig);
        recv_cr_single_sim(&src.lat, orig);
        recv_cr_single_sim(&src.dep, orig);
        recv_cr_single_sim(&src.weight, orig);
        recv_bool_single_sim(&src.is_teleseismic, orig);

        recv_i_single_sim(&src.n_rec_pair, orig);
    }

    void recv_rec_inter_sim(std::vector<SrcRec>& vrec, int orig){

        // recv receiver parameters of SrcRec class from orig
        for (auto& rec: vrec) {
            // id_src id_rec n_rec n_data
            // lon lat dep weight arr_time arr_time_ori t_adj
            // is_rec_pair n_rec_pair
            recv_i_single_sim(&rec.id_src, orig);
            recv_i_single_sim(&rec.id_rec, orig);
            recv_i_single_sim(&rec.n_rec, orig);
            recv_i_single_sim(&rec.n_data, orig);
            recv_cr_single_sim(&rec.lon, orig);
            recv_cr_single_sim(&rec.lat, orig);
            recv_cr_single_sim(&rec.dep, orig);
            recv_cr_single_sim(&rec.weight, orig);
            recv_cr_single_sim(&rec.arr_time, orig);
            recv_cr_single_sim(&rec.arr_time_ori, orig);
            recv_cr_single_sim(&rec.t_adj, orig);
            recv_bool_single_sim(&rec.is_rec_pair, orig);

            // if (is_rec_pair) id_rec_pair dep_pair lat_pair lon_pair ddt_adj_pair
            if (rec.is_rec_pair){
                // TODO : check how to recv double difference pairs
                for (int i = 0; i < 2; i++){
                    recv_i_single_sim(&rec.id_rec_pair[i],              orig);
                    recv_cr_single_sim(&rec.dep_pair[i],                orig);
                    recv_cr_single_sim(&rec.lat_pair[i],                orig);
                    recv_cr_single_sim(&rec.lon_pair[i],                orig);
                    recv_cr_single_sim(&rec.ddt_adj_pair[i],            orig);
                    recv_str_sim(rec.name_rec_pair[i],                  orig);
                    recv_cr_single_sim(&rec.station_correction_pair[i], orig);
                }

                recv_cr_single_sim(&rec.dif_arr_time,     orig);
                recv_cr_single_sim(&rec.dif_arr_time_ori, orig);

            }
        }

    }

};



#endif