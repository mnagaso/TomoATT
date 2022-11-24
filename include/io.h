#ifndef IO_H
#define IO_H

#include <fstream>
#include <vector>
#if USE_HDF5
    #include "hdf5.h"
#endif
#include "tinyxml2.h"

// to handle circular dependency
#pragma once
#include "grid.fwd.h"
//#include "source.fwd.h"
#include "io.fwd.h"

#include "utils.h"
#include "mpi_funcs.h"
#include "grid.h"
#include "config.h"

class IO_utils {

public:
    IO_utils();
    ~IO_utils();

    // initialize data output file
    void init_data_output_file();
    // finalize data/grid outpiut file (xdmf)
    void finalize_data_output_file(int);
    // change the output file name and xdmf objects
    void change_xdmf_obj(int);
    // output for updating xdmf file
    void update_xdmf_file(int);

    //
    // write functions
    //

    // xdmf write out for grid data
    void write_xdmf_file_grid();



    // write grid data to output file
    void write_grid(Grid&);
    // general function for write out data in hdf5
    void write_data_h5(Grid&, std::string&, std::string&, CUSTOMREAL*);
    // general function for write out data in hdf5
    void write_data_ascii(Grid&, std::string&, CUSTOMREAL*);
    // write true solution
    void write_true_solution(Grid&);
    // write velocity model
    void write_velocity_model_h5(Grid&);
    // write T0v
    void write_T0v(Grid&, int);
    // write u
    void write_u(Grid&);
    // write tau
    void write_tau(Grid&, int);
    // write temporal tau fields
    void write_tmp_tau(Grid&, int);
    // write result timetable T
    void write_T(Grid&, int);
    // write residual (resudual = true_solution - result)
    void write_residual(Grid&);
    // write adjoint field
    void write_adjoint_field(Grid&, int);
    // write fun_loc
    void write_fun(Grid&, int);
    // write a
    void write_a(Grid&, int);
    // write b
    void write_b(Grid&, int);
    // write c
    void write_c(Grid&, int);
    // write f
    void write_f(Grid&, int);
    // Ks
    void write_Ks(Grid&, int);
    // Ks_update
    void write_Ks_update(Grid&, int);
    // write all concerning parameters
    std::vector<CUSTOMREAL> get_grid_data(CUSTOMREAL * data);
    void write_concerning_parameters(Grid&, int);


    //
    // read functions
    //

    // read model data
    void read_model(std::string&, const char*, CUSTOMREAL*, int, int, int);
    // read Travel time from file for earthquake relocation
    void read_T(Grid&);

    void read_data_ascii(Grid&, std::string&);

private:
    // member variables
    std::string h5_output_grid_fname = "./out_data_grid.h5"; // output file name
    std::string h5_output_fname      = "./out_data.h5";      // output file name
    std::string xdmf_output_fname    = "./out_data.xmf";     // output file name
    std::string h5_group_name_grid   = "Mesh";               // group name for grid data
    std::string h5_group_name_event  = "Event";              // group name for event data
    std::string h5_group_name_data   = "Data";               // group name for event data

    std::string h5_dset_name_node_coords_x = "node_coords_x"; // dataset name for node coordinates
    std::string h5_dset_name_node_coords_y = "node_coords_y"; // dataset name for node coordinates
    std::string h5_dset_name_node_coords_z = "node_coords_z"; // dataset name for node coordinates
    std::string h5_dset_name_node_coords_p = "node_coords_p"; // dataset name for node coordinates
    std::string h5_dset_name_node_coords_t = "node_coords_t"; // dataset name for node coordinates
    std::string h5_dset_name_node_coords_r = "node_coords_r"; // dataset name for node coordinates
    std::string h5_dset_name_elem_conn     = "elem_conn";     // dataset name for element connectivity
    std::string h5_dset_name_procid        = "procid";        // dataset name for my_procs

    std::string h5_file_and_group  = h5_output_grid_fname + ":/" + h5_group_name_grid;
    std::string h5_whole_name_conn = h5_file_and_group + "/" + h5_dset_name_elem_conn;
    std::string h5_whole_name_x    = h5_file_and_group + "/" + h5_dset_name_node_coords_x;
    std::string h5_whole_name_y    = h5_file_and_group + "/" + h5_dset_name_node_coords_y;
    std::string h5_whole_name_z    = h5_file_and_group + "/" + h5_dset_name_node_coords_z;
    std::string h5_whole_name_proc = h5_file_and_group + "/" + h5_dset_name_procid;

    int nnodes_glob = 0;
    int nelms_glob  = 0;

    int array_rank_default = 2; // default array rank

    int custom_real_flag; //  custom real type

    // custom data type for real
    std::string type_name = "Float";
    std::string type_size = "4";

    // helper for xdmf/xml file io
    tinyxml2::XMLDocument* doc;
    tinyxml2::XMLNode* xdmf, *domain, *grid;
    void init_xdmf_file();
    void write_xdmf_file(int);
    void finalize_xdmf_file(int);
    void insert_data_xdmf(std::string& group, std::string& dset);
    // vectors for storing tinyxml2 objects
    std::vector<tinyxml2::XMLDocument*> doc_vec;
    std::vector<tinyxml2::XMLNode*>     xdmf_vec;
    std::vector<tinyxml2::XMLNode*>     domain_vec;
    std::vector<tinyxml2::XMLNode*>     grid_vec;
    std::vector<std::string>            fname_vec;
    std::vector<std::string>            xmfname_vec;
    void store_xdmf_obj();

#ifdef USE_HDF5
    // h5 variables
    hid_t  file_id;
    hid_t  group_id;
    hid_t  space_id;
    hid_t  dset_id;
    hid_t  plist_id, plist_id_dset;
    hid_t  file_dspace_id, mem_dspace_id;
    herr_t status;
    hid_t  plist_id_2d;
    hid_t  file_id_2d;

    // h5 low level utilities
    void h5_create_file_by_group_main(std::string&);
    void h5_open_file_by_group_main(std::string&);
    void h5_open_file_collective(std::string&);
    void h5_open_file_collective_input(std::string&);
    void h5_close_file_by_group_main();
    void h5_close_file_collective();

    void h5_create_group_by_group_main(std::string& );
    void h5_open_group_by_group_main(std::string& );
    void h5_open_group_collective(std::string& );
    void h5_close_group_by_group_main();
    void h5_close_group_collective();

    void h5_create_dataset(std::string&, int, int*, int);
    void h5_open_dataset(std::string&);
    void h5_open_dataset_no_group(std::string&);
    void h5_close_dataset();

    template <typename T>
    void h5_write_array(std::string&, int, int*, T*, int);

    template <typename T>
    void h5_read_array(std::string&, int, int*, T*, int*, bool);

    template <typename T>
    void h5_read_array_simple(std::string&,  T*);

    void read_data_h5(Grid&, CUSTOMREAL*, std::string, std::string);
#endif // HDF_HDF5

    // utilities
    std::string create_fname_ascii(std::string&);
    std::string create_fname_ascii_model(std::string&);


};

#endif // IO_H