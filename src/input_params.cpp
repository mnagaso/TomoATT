#include "input_params.h"
// #include "GaussEliminationWithPivoting.hpp"

InputParams::InputParams(std::string& input_file){

    if (world_rank == 0) {
        // parse input files
        YAML::Node config = YAML::LoadFile(input_file);

        // read domain information
        if (config["domain"]) {
            // minimum and maximum depth
            if (config["domain"]["min_max_dep"]) {
                min_dep = config["domain"]["min_max_dep"][0].as<CUSTOMREAL>();
                max_dep = config["domain"]["min_max_dep"][1].as<CUSTOMREAL>();
            }
            // minimum and maximum latitude
            if (config["domain"]["min_max_lat"]) {
                min_lat = config["domain"]["min_max_lat"][0].as<CUSTOMREAL>();
                max_lat = config["domain"]["min_max_lat"][1].as<CUSTOMREAL>();
            }
            // minimum and maximum longitude
            if (config["domain"]["min_max_lon"]) {
                min_lon = config["domain"]["min_max_lon"][0].as<CUSTOMREAL>();
                max_lon = config["domain"]["min_max_lon"][1].as<CUSTOMREAL>();
            }
            // number of grid nodes on each axis r(dep), t(lat), p(lon)
            if (config["domain"]["n_rtp"]) {
                ngrid_k = config["domain"]["n_rtp"][0].as<int>();
                ngrid_j = config["domain"]["n_rtp"][1].as<int>();
                ngrid_i = config["domain"]["n_rtp"][2].as<int>();
            }
        }

        if (config["source"]) {
            // source depth(km) lat lon
            if (config["source"]["src_dep_lat_lon"]) {
                src_dep = config["source"]["src_dep_lat_lon"][0].as<CUSTOMREAL>();
                src_lat = config["source"]["src_dep_lat_lon"][1].as<CUSTOMREAL>();
                src_lon = config["source"]["src_dep_lat_lon"][2].as<CUSTOMREAL>();
            }
            // src rec file
            if (config["source"]["src_rec_file"]) {
                src_rec_file_exist = true;
                src_rec_file = config["source"]["src_rec_file"].as<std::string>();
            }

            if (config["source"]["swap_src_rec"]) {
                int tmp_swap = config["source"]["swap_src_rec"].as<int>();
                if (tmp_swap == 1) {
                    swap_src_rec = true;
                } else {
                    swap_src_rec = false;
                }
            }
        }

        if (config["model"]) {
            // model type
            if (config["model"]["init_model_type"]) {
                init_model_type = config["model"]["init_model_type"].as<std::string>();
            }
            // model file path
            if (config["model"]["init_model_path"]) {
                init_model_path = config["model"]["init_model_path"].as<std::string>();
            }
            // model file path
            if (config["model"]["model_1d_name"]) {
                model_1d_name = config["model"]["model_1d_name"].as<std::string>();
            }
        }

        if (config["inversion"]) {
            // do inversion or not
            if (config["inversion"]["run_mode"]) {
                run_mode = config["inversion"]["run_mode"].as<int>();
            }
            // number of inversion grid
            if (config["inversion"]["n_inversion_grid"]) {
                n_inversion_grid = config["inversion"]["n_inversion_grid"].as<int>();
            }

            // number of inversion grid
            if (config["inversion"]["n_inv_dep_lat_lon"]) {
                n_inv_r = config["inversion"]["n_inv_dep_lat_lon"][0].as<int>();
                n_inv_t = config["inversion"]["n_inv_dep_lat_lon"][1].as<int>();
                n_inv_p = config["inversion"]["n_inv_dep_lat_lon"][2].as<int>();
            }

            // output_dir
            if (config["inversion"]["output_dir"]) {
                output_dir = config["inversion"]["output_dir"].as<std::string>();
            }

            // sta_correction_file
            if (config["inversion"]["sta_correction_file"]) {
                sta_correction_file_exist = true;
                sta_correction_file = config["inversion"]["sta_correction_file"].as<std::string>();
            }

            // type of input inversion grid
            if (config["inversion"]["type_dep_inv"]) {
                type_dep_inv = config["inversion"]["type_dep_inv"].as<int>();
            }
            if (config["inversion"]["type_lat_inv"]) {
                type_lat_inv = config["inversion"]["type_lat_inv"].as<int>();
            }
            if (config["inversion"]["type_lon_inv"]) {
                type_lon_inv = config["inversion"]["type_lon_inv"].as<int>();
            }

            // inversion grid positions
            if (config["inversion"]["min_max_dep_inv"]) {
                min_dep_inv = config["inversion"]["min_max_dep_inv"][0].as<CUSTOMREAL>();
                max_dep_inv = config["inversion"]["min_max_dep_inv"][1].as<CUSTOMREAL>();
            }
            // minimum and maximum latitude
            if (config["inversion"]["min_max_lat_inv"]) {
                min_lat_inv = config["inversion"]["min_max_lat_inv"][0].as<CUSTOMREAL>();
                max_lat_inv = config["inversion"]["min_max_lat_inv"][1].as<CUSTOMREAL>();
            }
            // minimum and maximum longitude
            if (config["inversion"]["min_max_lon_inv"]) {
                min_lon_inv = config["inversion"]["min_max_lon_inv"][0].as<CUSTOMREAL>();
                max_lon_inv = config["inversion"]["min_max_lon_inv"][1].as<CUSTOMREAL>();
            }

            // flexible inversion grid
            if (config["inversion"]["dep_inv"]) {
                n_inv_r_flex = config["inversion"]["dep_inv"].size();
                dep_inv = new CUSTOMREAL[n_inv_r_flex];
                for (int i = 0; i < n_inv_r_flex; i++){
                    dep_inv[i] = config["inversion"]["dep_inv"][i].as<CUSTOMREAL>();
                }
                n_inv_r_flex_read = true;
            }
            if (config["inversion"]["lat_inv"]) {
                n_inv_t_flex = config["inversion"]["lat_inv"].size();
                lat_inv = new CUSTOMREAL[n_inv_t_flex];
                for (int i = 0; i < n_inv_t_flex; i++){
                    lat_inv[i] = config["inversion"]["lat_inv"][i].as<CUSTOMREAL>();
                }
                n_inv_t_flex_read = true;
            }
            if (config["inversion"]["lon_inv"]) {
                n_inv_p_flex = config["inversion"]["lon_inv"].size();
                lon_inv = new CUSTOMREAL[n_inv_p_flex];
                for (int i = 0; i < n_inv_p_flex; i++){
                    lon_inv[i] = config["inversion"]["lon_inv"][i].as<CUSTOMREAL>();
                }
                n_inv_p_flex_read = true;
            }

            // nimber of max iteration for inversion
            if (config["inversion"]["max_iterations_inv"]) {
                max_iter_inv = config["inversion"]["max_iterations_inv"].as<int>();
            }
            // step_size
            if (config["inversion"]["step_size"]) {
                step_size_init = config["inversion"]["step_size"].as<CUSTOMREAL>();
            }
            if (config["inversion"]["step_size_sc"]) {
                step_size_init_sc = config["inversion"]["step_size_sc"].as<CUSTOMREAL>();
            }
            if (config["inversion"]["step_size_decay"]) {
                step_size_decay = config["inversion"]["step_size_decay"].as<CUSTOMREAL>();
            }
            
            // smoothing method
            if (config["inversion"]["smooth_method"]) {
                smooth_method = config["inversion"]["smooth_method"].as<int>();
                if (smooth_method > 1) {
                    std::cout << "undefined smooth_method. stop." << std::endl;
                    MPI_Finalize();
                    exit(1);
                }
            }
            // l_smooth_rtp
            if (config["inversion"]["l_smooth_rtp"]) {
                smooth_lr = config["inversion"]["l_smooth_rtp"][0].as<CUSTOMREAL>();
                smooth_lt = config["inversion"]["l_smooth_rtp"][1].as<CUSTOMREAL>();
                smooth_lp = config["inversion"]["l_smooth_rtp"][2].as<CUSTOMREAL>();
            }
            // optim_method
            if (config["inversion"]["optim_method"]) {
                optim_method = config["inversion"]["optim_method"].as<int>();
                if (optim_method > 2) {
                    std::cout << "undefined optim_method. stop." << std::endl;
                    //MPI_Finalize();
                    exit(1);
                }
            }
            // regularization weight
            if (config["inversion"]["regularization_weight"]) {
                regularization_weight = config["inversion"]["regularization_weight"].as<CUSTOMREAL>();
            }
            // max sub iteration
            if (config["inversion"]["max_sub_iterations"]) {
                max_sub_iterations = config["inversion"]["max_sub_iterations"].as<int>();
            }
        }

        if (config["inv_strategy"]) {
            // update which model parameters
            if (config["inv_strategy"]["is_inv_slowness"])
                is_inv_slowness = config["inv_strategy"]["is_inv_slowness"].as<int>();
            if (config["inv_strategy"]["is_inv_azi_ani"])
                is_inv_azi_ani = config["inv_strategy"]["is_inv_azi_ani"].as<int>();
            if (config["inv_strategy"]["is_inv_rad_ani"])
                is_inv_rad_ani = config["inv_strategy"]["is_inv_rad_ani"].as<int>();

            // taper kernel (now only for teleseismic tomography)
            if (config["inv_strategy"]["kernel_taper"]){
                kernel_taper[0] = config["inv_strategy"]["kernel_taper"][0].as<CUSTOMREAL>();
                kernel_taper[1] = config["inv_strategy"]["kernel_taper"][1].as<CUSTOMREAL>();
            }

            // station correction (now only for teleseismic data)
            if (config["inv_strategy"]["is_sta_correction"]){
                is_sta_correction = config["inv_strategy"]["is_sta_correction"].as<int>();
            }
        }


        if (config["parallel"]) {
            // number of simultaneous runs
            if(config["parallel"]["n_sims"]) {
                n_sims = config["parallel"]["n_sims"].as<int>();
            }
            // number of subdomains
            if (config["parallel"]["ndiv_rtp"]) {
                ndiv_k = config["parallel"]["ndiv_rtp"][0].as<int>();
                ndiv_j = config["parallel"]["ndiv_rtp"][1].as<int>();
                ndiv_i = config["parallel"]["ndiv_rtp"][2].as<int>();
            }
            // number of processes in each subdomain
            if (config["parallel"]["nproc_sub"]) {
                n_subprocs = config["parallel"]["nproc_sub"].as<int>();
            }
            // gpu flag
            if (config["parallel"]["use_gpu"]) {
                use_gpu = config["parallel"]["use_gpu"].as<int>();
            }
        }

        if (config["calculation"]) {
            // convergence tolerance
            if (config["calculation"]["convergence_tolerance"]) {
                conv_tol = config["calculation"]["convergence_tolerance"].as<CUSTOMREAL>();
            }
            // max number of iterations
            if (config["calculation"]["max_iterations"]) {
                max_iter = config["calculation"]["max_iterations"].as<int>();
            }
            // stencil order
            if (config["calculation"]["stencil_order"]) {
                stencil_order = config["calculation"]["stencil_order"].as<int>();
                // check if stencil_order == 999 : hybrid scheme
                if (stencil_order == 999) {
                    hybrid_stencil_order = true;
                    stencil_order = 1;
                }
            }
            // stencil type
            if (config["calculation"]["stencil_type"]) {
                stencil_type = config["calculation"]["stencil_type"].as<int>();
            }
            // sweep type
            if (config["calculation"]["sweep_type"]) {
                sweep_type = config["calculation"]["sweep_type"].as<int>();
            }
            // output file format
            if (config["calculation"]["output_file_format"]) {
                int ff_flag = config["calculation"]["output_file_format"].as<int>();
                if (ff_flag == 0){
                    #if USE_HDF5
                        output_format = OUTPUT_FORMAT_HDF5;
                    #else
                        std::cout << "output_file_format is 0, but the code is compiled without HDF5. stop." << std::endl;
                        MPI_Finalize();
                        exit(1);
                    #endif
                } else if (ff_flag == 1){
                    output_format = OUTPUT_FORMAT_ASCII;
                } else {
                    std::cout << "undefined output_file_format. stop." << std::endl;
                    MPI_Finalize();
                    exit(1);
                }
            }
        }

        if (config["output_setting"]) {
            if (config["output_setting"]["is_output_source_field"])
                is_output_source_field = config["output_setting"]["is_output_source_field"].as<int>();
            if (config["output_setting"]["is_output_model_dat"])
                is_output_model_dat = config["output_setting"]["is_output_model_dat"].as<int>();
            if (config["output_setting"]["is_verbose_output"])
                is_verbose_output = config["output_setting"]["is_verbose_output"].as<int>();
            if (config["output_setting"]["is_output_final_model"])
                is_output_final_model = config["output_setting"]["is_output_final_model"].as<int>();
        }

        if (config["debug"]) {
            if (config["debug"]["debug_mode"]) {
                int tmp_test = config["debug"]["debug_mode"].as<int>();
                if (tmp_test == 1) {
                    if_test = true;
                } else {
                    if_test = false;
                }
            }
        }

        std::cout << "min_max_dep: " << min_dep << " " << max_dep << std::endl;
        std::cout << "min_max_lat: " << min_lat << " " << max_lat << std::endl;
        std::cout << "min_max_lon: " << min_lon << " " << max_lon << std::endl;
        std::cout << "n_rtp: "    << ngrid_k << " " << ngrid_j << " " << ngrid_i << std::endl;
        std::cout << "ndiv_rtp: " << ndiv_k << " "  << ndiv_j  << " " << ndiv_i << std::endl;
        std::cout << "n_subprocs: " << n_subprocs << std::endl;
        std::cout << "n_sims: " << n_sims << std::endl;

        // set inversion grid definition if not set by user
        if (min_dep_inv <= -99999) min_dep_inv = min_dep;
        if (max_dep_inv <= -99999) max_dep_inv = max_dep;
        if (min_lat_inv <= -99999) min_lat_inv = min_lat;
        if (max_lat_inv <= -99999) max_lat_inv = max_lat;
        if (min_lon_inv <= -99999) min_lon_inv = min_lon;
        if (max_lon_inv <= -99999) max_lon_inv = max_lon;

        // allocate dummy arrays for flex inv grid
        if (!n_inv_r_flex_read)
            dep_inv = new CUSTOMREAL[n_inv_r_flex];
        if (!n_inv_t_flex_read)
            lat_inv = new CUSTOMREAL[n_inv_t_flex];
        if (!n_inv_p_flex_read)
            lon_inv = new CUSTOMREAL[n_inv_p_flex];
    }

    std::cout << "parameter file read done." << std::endl;

    synchronize_all_world();

    // broadcast all the values read
    broadcast_cr_single(min_dep, 0);
    broadcast_cr_single(max_dep, 0);
    broadcast_cr_single(min_lat, 0);
    broadcast_cr_single(max_lat, 0);
    broadcast_cr_single(min_lon, 0);
    broadcast_cr_single(max_lon, 0);

    broadcast_i_single(ngrid_i, 0);
    broadcast_i_single(ngrid_j, 0);
    broadcast_i_single(ngrid_k, 0);

    broadcast_cr_single(src_dep, 0);
    broadcast_cr_single(src_lat, 0);
    broadcast_cr_single(src_lon, 0);
    if (src_rec_file_exist == false){
        SrcRec src;
        src.id_src = 0;
        src.lat    = src_lat;
        src.lon    = src_lon;
        src.dep    = src_dep;
        src_points.push_back(src);
        src_ids_this_sim.push_back(0);
        SrcRec rec;
        rec_points.push_back({rec});
    }
    broadcast_bool_single(swap_src_rec, 0);

    broadcast_str(src_rec_file, 0);
    broadcast_str(sta_correction_file, 0);
    broadcast_str(output_dir, 0);
    broadcast_bool_single(src_rec_file_exist, 0);
    broadcast_bool_single(sta_correction_file_exist, 0);
    broadcast_str(init_model_type, 0);
    broadcast_str(init_model_path, 0);
    broadcast_i_single(run_mode, 0);
    broadcast_i_single(n_inversion_grid, 0);
    broadcast_i_single(n_inv_r, 0);
    broadcast_i_single(n_inv_t, 0);
    broadcast_i_single(n_inv_p, 0);
    broadcast_cr_single(min_dep_inv, 0);
    broadcast_cr_single(max_dep_inv, 0);
    broadcast_cr_single(min_lat_inv, 0);
    broadcast_cr_single(max_lat_inv, 0);
    broadcast_cr_single(min_lon_inv, 0);
    broadcast_cr_single(max_lon_inv, 0);

    broadcast_i_single(type_dep_inv, 0);
    broadcast_i_single(type_lat_inv, 0);
    broadcast_i_single(type_lon_inv, 0);
    broadcast_i_single(n_inv_r_flex, 0);
    broadcast_i_single(n_inv_t_flex, 0);
    broadcast_i_single(n_inv_p_flex, 0);
    if (world_rank != 0) {
        dep_inv = new CUSTOMREAL[n_inv_r_flex];
        lat_inv = new CUSTOMREAL[n_inv_t_flex];
        lon_inv = new CUSTOMREAL[n_inv_p_flex];
    }

    broadcast_cr(dep_inv,n_inv_r_flex, 0);
    broadcast_cr(lat_inv,n_inv_t_flex, 0);
    broadcast_cr(lon_inv,n_inv_p_flex, 0);

    broadcast_i_single(smooth_method, 0);
    broadcast_cr_single(smooth_lr, 0);
    broadcast_cr_single(smooth_lt, 0);
    broadcast_cr_single(smooth_lp, 0);
    broadcast_i_single(optim_method, 0);
    broadcast_i_single(max_iter_inv, 0);
    broadcast_cr_single(step_size_init, 0);
    broadcast_cr_single(step_size_init_sc, 0);
    broadcast_cr_single(step_size_decay, 0);
    broadcast_cr_single(regularization_weight, 0);
    broadcast_i_single(max_sub_iterations, 0);
    broadcast_i_single(ndiv_i, 0);
    broadcast_i_single(ndiv_j, 0);
    broadcast_i_single(ndiv_k, 0);
    broadcast_i_single(n_subprocs, 0);
    broadcast_i_single(n_sims, 0);
    broadcast_cr_single(conv_tol, 0);
    broadcast_i_single(max_iter, 0);
    broadcast_i_single(stencil_order, 0);
    broadcast_bool_single(hybrid_stencil_order, 0);
    broadcast_i_single(stencil_type, 0);
    broadcast_i_single(sweep_type, 0);
    broadcast_i_single(output_format, 0);
    broadcast_bool_single(if_test, 0);
    broadcast_i_single(use_gpu, 0);

    broadcast_bool_single(is_output_source_field, 0);
    broadcast_bool_single(is_output_model_dat, 0);
    broadcast_bool_single(is_verbose_output, 0);
    broadcast_bool_single(is_output_final_model, 0);
    broadcast_bool_single(is_inv_slowness, 0);
    broadcast_bool_single(is_inv_azi_ani, 0);
    broadcast_bool_single(is_inv_rad_ani, 0);
    broadcast_cr(kernel_taper,2,0);
    broadcast_bool_single(is_sta_correction, 0);


    // read station correction file by all processes
    if (sta_correction_file_exist) {
        for (int i_proc = 0; i_proc < world_nprocs; i_proc++){
            if (i_proc == world_rank) {
                // store all src/rec info
                parse_sta_correction_file();
            }
            synchronize_all_world();
        }
    }

    // read src rec file by all processes #TODO: check if only subdomain's main requires this
    if (src_rec_file_exist) {
        for (int i_proc = 0; i_proc < world_nprocs; i_proc++){
            if (i_proc == world_rank) {
                // store all src/rec info
                parse_src_rec_file();
            }
            synchronize_all_world();
        }
        // define src/rec file name for output
        // size_t ext = src_rec_file.find_last_of(".");
        // src_rec_file_out = src_rec_file.substr(0, ext) + "_out.dat";

        // backup original src/rec list
        src_points_back = src_points;
        rec_points_back = rec_points;

        // new version backup original src/rec list
        src_list_back_nv = src_list_nv;
        rec_list_back_nv = rec_list_nv;
        data_info_back_nv = data_info_nv;

        // check if src positions are within the domain or not (teleseismic source)
        // detected teleseismic source is separated into tele_src_points and tele_rec_points
        separate_region_and_tele_src();

        if (swap_src_rec) {
            // here only reginal events will be processed
            stdout_by_main("###### Swapping src and rec. (only regional events will be processed) ######\n");
            do_swap_src_rec();
        }

        // rearrange and analize "data_info_nv" -> data_info_smap; 
        rearrange_data_info();

        // generate src_list_prepare_nv based on data_info_smap
        generate_src_list_prepare();

        // generate syn_time_list_sr
        generate_syn_time_list();

        // concatenate resional and teleseismic src/rec points
        merge_region_and_tele_src();
    }

    

    // check contradictory settings
    check_contradictions();

    // broadcast the values to all processes
    stdout_by_main("read input file successfully.");

}


InputParams::~InputParams(){
    // free memory
    if (subdom_main) {
        for (auto& id_src: src_ids_this_sim){
            if (src_points[id_src].is_teleseismic){
                if (j_last)
                    free(src_points[id_src].arr_times_bound_N);
                if (j_first)
                    free(src_points[id_src].arr_times_bound_S);
                if (i_last)
                    free(src_points[id_src].arr_times_bound_E);
                if (i_first)
                    free(src_points[id_src].arr_times_bound_W);
                if (k_first)
                    free(src_points[id_src].arr_times_bound_Bot);

                free(src_points[id_src].is_bound_src);
            }
        }
    }

    delete [] dep_inv;
    delete [] lat_inv;
    delete [] lon_inv;

    // clear all src/rec points
    src_points.clear();
    rec_points.clear();
    src_points_back.clear();
    rec_points_back.clear();
    src_points_out.clear();
    rec_points_out.clear();
}


// return radious
CUSTOMREAL InputParams::get_src_radius() {
    if (src_rec_file_exist)
        return depth2radius(src_points[id_sim_src].dep);
    else
        return depth2radius(src_dep);
}


CUSTOMREAL InputParams::get_src_lat() {
    if (src_rec_file_exist)
        return src_points[id_sim_src].lat*DEG2RAD;
    else
        return src_lat*DEG2RAD;
}


CUSTOMREAL InputParams::get_src_lon() {
    if (src_rec_file_exist)
        return src_points[id_sim_src].lon*DEG2RAD;
    else
        return src_lon*DEG2RAD;
}


SrcRec& InputParams::get_src_point(int i_src){
    return src_points[i_src];
}


std::vector<SrcRec>& InputParams::get_rec_points(int id_src) {
    return rec_points[id_src];
}



// return radious
CUSTOMREAL InputParams::get_src_radius_nv() {
    if (src_rec_file_exist)
        return depth2radius(src_list_nv[name_sim_src].dep);
    else
        return depth2radius(src_dep);
}


CUSTOMREAL InputParams::get_src_lat_nv() {
    if (src_rec_file_exist)
        return src_list_nv[name_sim_src].lat*DEG2RAD;
    else
        return src_lat*DEG2RAD;
}


CUSTOMREAL InputParams::get_src_lon_nv() {
    if (src_rec_file_exist)
        return src_list_nv[name_sim_src].lon*DEG2RAD;
    else
        return src_lon*DEG2RAD;
}

SrcRecInfo& InputParams::get_src_point_nv(std::string name_src){
    return src_list_nv[name_src];
}

SrcRecInfo& InputParams::get_rec_point_nv(std::string name_rec){
    return rec_list_nv[name_rec];
}

std::vector<std::string> InputParams::get_rec_points_nv(std::string name_src){
    std::vector<std::string> recs;
    for(auto iter = syn_time_list_sr[name_src].begin(); iter != syn_time_list_sr[name_src].end(); iter++){
        recs.push_back(iter->first);
    }
    return recs;
}

//
// functions for processing src_rec_file
// #TODO: functions concerinng SrcRec object may be moved to another file
//

void InputParams::parse_src_rec_file(){
    // #TODO: add error handling

    std::ifstream ifs(src_rec_file);
    std::string line;
    int cc =0; // count the number of lines
    int i_src_now = 0; // count the number of srcs
    int i_rec_now = 0; // count the number of receivers
    int ndata_tmp = 0; // count the number of receivers or differential traveltime data for each source
    std::vector<SrcRec> rec_points_tmp;
    src_points.clear();
    rec_points_tmp.clear();
    rec_points.clear();

    src_list_nv.clear();
    rec_list_nv.clear();
    data_info_nv.clear();

    std::string src_name;
    CUSTOMREAL src_weight = 1.0;
    CUSTOMREAL rec_weight = 1.0;
    int src_id = -1;
    while (std::getline(ifs, line)) {
        // skip comment and empty lines
        if (line[0] == '#' || line.empty())
            continue;

        // erase the trailing space
        line.erase(line.find_last_not_of(" \n\r\t")+1);

        // parse the line with arbitrary number of spaces
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;

        while (std::getline(ss, token, ' ')) {
            if (token.size() > 0) // skip the first spaces and multiple spaces
                tokens.push_back(token);
        }

        // store values into structure
        
        
        
        if (cc == 0){
            SrcRec src;
            //src.id_src     = std::stoi(tokens[0]);
            src.id_src     = i_src_now; // MNMN: here use id_src of active source lines order of src rec file, which allow to comment out bad events.
            src.year       = std::stoi(tokens[1]);
            src.month      = std::stoi(tokens[2]);
            src.day        = std::stoi(tokens[3]);
            src.hour       = std::stoi(tokens[4]);
            src.min        = std::stoi(tokens[5]);
            src.sec        = static_cast<CUSTOMREAL>(std::stod(tokens[6]));
            src.lat        = static_cast<CUSTOMREAL>(std::stod(tokens[7])); // in degree
            src.lon        = static_cast<CUSTOMREAL>(std::stod(tokens[8])); // in degree
            src.dep        = static_cast<CUSTOMREAL>(std::stod(tokens[9])); // source in km
            src.mag        = static_cast<CUSTOMREAL>(std::stod(tokens[10]));
            src.n_data     = std::stoi(tokens[11]);
            ndata_tmp      = src.n_data;
            src.n_rec      = 0;
            src.n_rec_pair = 0;
            src.name_src = tokens[12];
            // check if tokens[13] exists, then read weight
            if (tokens.size() > 13)
                src.weight = static_cast<CUSTOMREAL>(std::stod(tokens[13]));
            else
                src.weight = 1.0; // default weight
            src_points.push_back(src);
            cc++;


            // new version of source
            SrcRecInfo src_nv;
            src_nv.id         = std::stoi(tokens[0]);
            src_nv.year       = std::stoi(tokens[1]);
            src_nv.month      = std::stoi(tokens[2]);
            src_nv.day        = std::stoi(tokens[3]);
            src_nv.hour       = std::stoi(tokens[4]);
            src_nv.min        = std::stoi(tokens[5]);
            src_nv.sec        = static_cast<CUSTOMREAL>(std::stod(tokens[6]));
            src_nv.lat        = static_cast<CUSTOMREAL>(std::stod(tokens[7])); // in degree
            src_nv.lon        = static_cast<CUSTOMREAL>(std::stod(tokens[8])); // in degree
            src_nv.dep        = static_cast<CUSTOMREAL>(std::stod(tokens[9])); // source in km
            src_nv.mag        = static_cast<CUSTOMREAL>(std::stod(tokens[10]));
            src_nv.n_data     = std::stoi(tokens[11]);
            src_nv.name       = tokens[12];

            
            if (tokens.size() > 13)
                src_weight = static_cast<CUSTOMREAL>(std::stod(tokens[13]));
            else
                src_weight = 1.0; // default weight

            if (src_list_nv.find(src_nv.name) == src_list_nv.end()){   // new source
                src_list_nv[src_nv.name] = src_nv;
            }
            
            src_id = src_nv.id;
            src_name = src_nv.name;

        } else {    

            // read single receiver or differential traveltime data
            if (tokens.size() < 11) {

                // read absolute traveltime
                SrcRec rec;
                //rec.id_src   = std::stoi(tokens[0]);
                rec.id_src   = i_src_now; // MNMN: here use id_src of active source lines order of src rec file, which allow to comment out bad events.
                //rec.id_rec   = std::stoi(tokens[1]);
                rec.id_rec   = i_rec_now; // MNMN: here use id_rec of active receiver lines order of src rec file, which allow to comment out bad stations.
                rec.name_rec = tokens[2];
                rec.lat      = static_cast<CUSTOMREAL>(std::stod(tokens[3])); // in degree
                rec.lon      = static_cast<CUSTOMREAL>(std::stod(tokens[4])); // in degree
                rec.dep      = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[5])/1000.0); // convert elevation in meter to depth in km
                rec.phase    = tokens[6];
                rec.arr_time_ori = static_cast<CUSTOMREAL>(std::stod(tokens[7])); // store read data

                // check if tokens[8] exists read weight
                if (tokens.size() > 8)
                    rec.weight = static_cast<CUSTOMREAL>(std::stod(tokens[8]));
                else
                    rec.weight = 1.0; // default weight
                rec.is_src_rec = true;

                rec_points_tmp.push_back(rec);
                cc++;
                src_points.at(src_points.size()-1).n_rec++;

                if (rec_list.find(rec.name_rec) == rec_list.end()){
                    // a new receiver
                    rec_list[rec.name_rec] = rec;
                    N_receiver += 1;
                }


                // ----- new version of receiver and data -----
                // receiver
                SrcRecInfo rec_nv;
                rec_nv.id       = std::stoi(tokens[1]);
                rec_nv.name     = tokens[2];
                rec_nv.lat      = static_cast<CUSTOMREAL>(std::stod(tokens[3])); // in degree
                rec_nv.lon      = static_cast<CUSTOMREAL>(std::stod(tokens[4])); // in degree
                rec_nv.dep      = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[5])/1000.0); // convert elevation in meter to depth in km

                if(rec_list_nv.find(rec_nv.name) == rec_list_nv.end()){     // new receiver
                    rec_list_nv[rec_nv.name] = rec_nv;
                }

                // traveltime data
                DataInfo data_nv;
                if (tokens.size() > 8)
                    rec_weight = static_cast<CUSTOMREAL>(std::stod(tokens[8]));
                else
                    rec_weight = 1.0;

                data_nv.data_weight = src_weight * rec_weight;
                data_nv.phase = tokens[6];

                data_nv.is_src_rec = true;
                data_nv.id_src      = src_id;
                data_nv.name_src    = src_name;
                data_nv.id_rec      = rec_nv.id;
                data_nv.name_rec    = rec_nv.name;
                data_nv.travel_time_obs = static_cast<CUSTOMREAL>(std::stod(tokens[7])); // store read data
                
                data_info_nv.push_back(data_nv);

            } else {

                // read differential traveltime
                SrcRec rec;
                //rec.id_src    = std::stoi(tokens[0]);
                rec.id_src   = i_src_now; // MNMN: here use id_src of active source lines order of src rec file, which allow to comment out bad events.
                rec.id_rec_pair[0] = std::stoi(tokens[1]);
                rec.name_rec_pair[0] = tokens[2];
                rec.lat_pair[0] = static_cast<CUSTOMREAL>(std::stod(tokens[3]));        // in degree
                rec.lon_pair[0] = static_cast<CUSTOMREAL>(std::stod(tokens[4]));        // in degree
                rec.dep_pair[0] = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[5])/1000.0);        // convert elevation in meter to depth in km
                rec.id_rec_pair[1] = std::stoi(tokens[6]);
                rec.name_rec_pair[1] = tokens[7];
                rec.lat_pair[1]      = static_cast<CUSTOMREAL>(std::stod(tokens[8])); // in degree
                rec.lon_pair[1]      = static_cast<CUSTOMREAL>(std::stod(tokens[9])); // in degree
                rec.dep_pair[1]      = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[10])/1000.0); // convert elevation in meter to depth in km
                rec.phase    = tokens[11];
                // rec.dif_arr_time = static_cast<CUSTOMREAL>(std::stod(tokens[12]));
                rec.dif_arr_time = 0.0;
                rec.dif_arr_time_ori = static_cast<CUSTOMREAL>(std::stod(tokens[12])); // store read data

                // check if tokens[9] exists read weight
                if (tokens.size() > 13)
                    rec.weight = static_cast<CUSTOMREAL>(std::stod(tokens[13]));
                else
                    rec.weight = 1.0; // default weight
                rec.is_rec_pair = true;
                
                if (rec_list.find(rec.name_rec_pair[0]) == rec_list.end()){
                    // a new receiver
                    SrcRec tmp_rec;
                    tmp_rec.id_rec   = std::stoi(tokens[1]);
                    tmp_rec.name_rec = tokens[2];
                    tmp_rec.lat      = static_cast<CUSTOMREAL>(std::stod(tokens[3])); // in degree
                    tmp_rec.lon      = static_cast<CUSTOMREAL>(std::stod(tokens[4])); // in degree
                    tmp_rec.dep      = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[5])/1000.0); // convert elevation in meter to depth in km
                    rec_list[rec.name_rec_pair[0]] = tmp_rec;

                    rec.station_correction_pair[0] = 0.0;
                    station_correction[rec.name_rec_pair[0]] = 0.0;
                    station_correction_kernel[rec.name_rec_pair[0]] = 0.0;
                    N_receiver += 1;
                } else {
                    // station exists in the rec_list
                    rec.station_correction_pair[0] = station_correction[rec.name_rec_pair[0]];
                    // std::cout << "station exist, " << rec.name_rec_pair[0] << ", correction: " << rec.station_correction_pair[0] << std::endl;
                }
                if (rec_list.find(rec.name_rec_pair[1]) == rec_list.end()){
                    // a new receiver
                    SrcRec tmp_rec;
                    tmp_rec.id_rec   = std::stoi(tokens[6]);
                    tmp_rec.name_rec = tokens[7];
                    tmp_rec.lat      = static_cast<CUSTOMREAL>(std::stod(tokens[8])); // in degree
                    tmp_rec.lon      = static_cast<CUSTOMREAL>(std::stod(tokens[9])); // in degree
                    tmp_rec.dep      = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[10])/1000.0); // convert elevation in meter to depth in km
                    rec_list[rec.name_rec_pair[1]] = tmp_rec;

                    station_correction[rec.name_rec_pair[1]] = 0.0;
                    station_correction_kernel[rec.name_rec_pair[1]] = 0.0;  
                    N_receiver += 1;
                } else {
                    // station exists in the rec_list
                    rec.station_correction_pair[1] = station_correction[rec.name_rec_pair[1]];
                    // std::cout << "station exist, " << rec.name_rec_pair[1] << ", correction: " << rec.station_correction_pair[1] << std::endl;
                 }

                rec_points_tmp.push_back(rec);
                cc++;
                src_points.at(src_points.size()-1).n_rec_pair++;



                // ----- new version of receiver and data -----
                // receiver
                SrcRecInfo rec_nv1;
                rec_nv1.id       = std::stoi(tokens[1]);
                rec_nv1.name     = tokens[2];
                rec_nv1.lat      = static_cast<CUSTOMREAL>(std::stod(tokens[3])); // in degree
                rec_nv1.lon      = static_cast<CUSTOMREAL>(std::stod(tokens[4])); // in degree
                rec_nv1.dep      = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[5])/1000.0); // convert elevation in meter to depth in km

                if(rec_list_nv.find(rec_nv1.name) == rec_list_nv.end()){     // new receiver
                    rec_list_nv[rec_nv1.name] = rec_nv1;
                }

                SrcRecInfo rec_nv2;
                rec_nv2.id       = std::stoi(tokens[6]);
                rec_nv2.name     = tokens[7];
                rec_nv2.lat      = static_cast<CUSTOMREAL>(std::stod(tokens[8])); // in degree
                rec_nv2.lon      = static_cast<CUSTOMREAL>(std::stod(tokens[9])); // in degree
                rec_nv2.dep      = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[10])/1000.0); // convert elevation in meter to depth in km

                if(rec_list_nv.find(rec_nv2.name) == rec_list_nv.end()){     // new receiver
                    rec_list_nv[rec_nv2.name] = rec_nv2;
                }


                // common source differential traveltime 
                DataInfo data_nv;
                if (tokens.size() > 13)
                    rec_weight = static_cast<CUSTOMREAL>(std::stod(tokens[8]));
                else
                    rec_weight = 1.0;
                data_nv.data_weight = src_weight * rec_weight;
                data_nv.phase = tokens[11];

                data_nv.is_rec_pair         = true;
                data_nv.id_src_single       = src_id;
                data_nv.name_src_single     = src_name;
                data_nv.id_rec_pair         = {rec_nv1.id, rec_nv2.id};
                data_nv.name_rec_pair       = {rec_nv1.name, rec_nv2.name};
                data_nv.cs_dif_travel_time_obs = static_cast<CUSTOMREAL>(std::stod(tokens[12])); // store read data
                
                data_info_nv.push_back(data_nv);
            }

            if (cc > ndata_tmp) {
                // go to the next source
                cc = 0;
                rec_points.push_back(rec_points_tmp);
                rec_points_tmp.clear();
                i_src_now++;
                i_rec_now = 0;
            } else {
                i_rec_now++;
            }
        }

    /*
        // print for DEBUG
        for (auto& t : tokens) {
            std::cout << t << "---";
        }
        std::cout << std::endl;
    */

    }

    // abort if number of src_points are less than n_sims
    int n_src_points = src_points.size();
    if (n_src_points < n_sims){
        std::cout << "Error: number of sources in src_rec_file is less than n_sims. Abort." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }


    // check new version of src rec data
    if (1 < 0){     // check by Chen Jing
        for(auto iter = src_list_nv.begin(); iter != src_list_nv.end(); iter++){
            std::cout   << "source id: " << iter->second.id 
                        << ", source name: " << iter->second.name
                        << std::endl;
        }

        for(auto iter = rec_list_nv.begin(); iter != rec_list_nv.end(); iter++){
            std::cout   << "receiver id: " << iter->second.id 
                        << ", receiver name: " << iter->second.name
                        << std::endl;
        } 

        for(int i = 0; i < (int)data_info_nv.size(); i++){
            if (data_info_nv[i].is_src_rec){
                std::cout   << "absolute traveltime: " << data_info_nv[i].travel_time_obs
                            << ", source name: " << data_info_nv[i].name_src
                            << ", receiver name: " << data_info_nv[i].name_rec
                            << std::endl; 
            }   
            if (data_info_nv[i].is_rec_pair){
                std::cout   << "common source differential traveltime: " << data_info_nv[i].cs_dif_travel_time_obs
                            << ", source name: " << data_info_nv[i].name_src_single
                            << ", receiver pair name: " << data_info_nv[i].name_rec_pair[0] 
                            << ", " << data_info_nv[i].name_rec_pair[1]
                            << std::endl; 
            }
        }
        std::cout << data_info_nv.size() << std::endl;
    }
}

void InputParams::parse_sta_correction_file(){
    std::ifstream ifs(sta_correction_file);
    std::string line;

    while(std::getline(ifs,line)) {
        if(line[0] == '#' || line.empty())
            continue;
        
        line.erase(line.find_last_not_of(" \n\r\t")+1);
    
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;

        while (std::getline(ss,token,' ')) {
            if (token.size() > 0) 
                tokens.push_back(token);
        }

        // store station corrections into rec_list

        std::string tmp_sta_name = tokens[0];
        CUSTOMREAL tmp_correct = static_cast<CUSTOMREAL>(std::stod(tokens[4]));

        if (rec_list.find(tmp_sta_name) == rec_list.end()){
            // new station
            SrcRec tmp_rec;
            tmp_rec.name_rec = tmp_sta_name;
            tmp_rec.lat      = static_cast<CUSTOMREAL>(std::stod(tokens[1])); // in degree
            tmp_rec.lon      = static_cast<CUSTOMREAL>(std::stod(tokens[2])); // in degree
            tmp_rec.dep      = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[3])/1000.0); // convert elevation in meter to depth in km
            rec_list[tmp_sta_name] = tmp_rec;
            station_correction[tmp_sta_name] = tmp_correct;
            station_correction_kernel[tmp_sta_name] = 0.0;
            N_receiver += 1;
        } else {
            // pre exist station
            station_correction[tmp_sta_name] = tmp_correct;
            station_correction_kernel[tmp_sta_name] = 0.0;
        }
    }

}

void InputParams::do_swap_src_rec(){

    // swap src/rec points
    // at this moment, all the sources are divided into src_points (regional) and tele_src_points (teleseismic)

    std::vector<SrcRec> new_srcs; // new src points
    std::vector<std::vector<SrcRec>> new_recs;

    // generate new source list
    for (long unsigned int i_src = 0; i_src < src_points.size(); i_src++) {
        for(long unsigned int i_rec = 0; i_rec < rec_points[i_src].size(); i_rec++) {

            if (new_srcs.size() == 0){
                new_srcs.push_back(rec_points[i_src][i_rec]);
                new_srcs.back().id_srcs_ori.push_back(rec_points[i_src][i_rec].id_src);
            } else if (new_srcs.size() != 0) {

                bool found = false;
                // check if the existing element has same name
                for (long unsigned int i_new_src = 0; i_new_src < new_srcs.size(); i_new_src++) {
                    if (new_srcs[i_new_src].name_rec.compare(rec_points[i_src][i_rec].name_rec) == 0) {
                        new_srcs[i_new_src].id_srcs_ori.push_back(rec_points[i_src][i_rec].id_src);
                        found=true;
                        break;
                    }
                    else {
                    }
                }
                if (!found) {
                    new_srcs.push_back(rec_points[i_src][i_rec]);
                    new_srcs.back().id_srcs_ori.push_back(rec_points[i_src][i_rec].id_src);
                }

            }
        }
    }

    // generate new rec list
    for (long unsigned int i_src = 0; i_src < new_srcs.size(); i_src++) {

        std::vector<SrcRec> tmp_list_recs;

        for (auto& i_src_ori : new_srcs[i_src].id_srcs_ori){

            SrcRec tmp_new_rec = src_points_back[i_src_ori];        // now src_points_back is the copy of src_points

            for (auto& tmp_rec_ori : rec_points_back[i_src_ori]){
                if (tmp_rec_ori.name_rec.compare(new_srcs[i_src].name_rec)==0) {
                    // we can use the same arrival time for a src-rec pair by principle of reciprocity (Aki & Richards, 2002).
                    tmp_new_rec.arr_time     = tmp_rec_ori.arr_time;
                    tmp_new_rec.arr_time_ori = tmp_rec_ori.arr_time_ori;
                    tmp_new_rec.id_rec_ori   = tmp_rec_ori.id_rec;
                    goto rec_found;
                }
            }

            rec_found:
                tmp_list_recs.push_back(tmp_new_rec);
        }

        new_recs.push_back(tmp_list_recs);

    }


    // backup and set new src/rec points
    // !! ONLY REGIONAL EVENTS AND RECESVERS ARE STORED IN *_back VECTORS !!
    src_points.clear();
    rec_points.clear();
    // teleseismic events are concatenate to the vectors below later
    src_points = new_srcs;
    rec_points = new_recs;



    // for new version
    std::map<std::string, SrcRecInfo> tmp_src_rec_list_nv = src_list_nv;
    src_list_nv = rec_list_nv;
    rec_list_nv = tmp_src_rec_list_nv;

    for(int i = 0; i < (int)data_info_nv.size(); i++){
        DataInfo tmp_data = data_info_nv[i];
        if (tmp_data.is_src_rec){   // absolute traveltime  ->  absolute traveltime
            tmp_data.id_src = data_info_nv[i].id_rec;
            tmp_data.name_src = data_info_nv[i].name_rec;
            tmp_data.id_rec = data_info_nv[i].id_src;
            tmp_data.name_rec = data_info_nv[i].name_src;
        } else if (tmp_data.is_rec_pair) {      // common source differential traveltime  ->  common receiver differential traveltime
            tmp_data.is_rec_pair = false;
            tmp_data.is_src_pair = true;
            tmp_data.id_src_pair = data_info_nv[i].id_rec_pair;
            tmp_data.name_src_pair = data_info_nv[i].name_rec_pair;
            tmp_data.id_rec_single = data_info_nv[i].id_src_single; 
            tmp_data.name_rec_single = data_info_nv[i].name_src_single; 
            tmp_data.cr_dif_travel_time_obs = data_info_nv[i].cs_dif_travel_time_obs;
        } else if (tmp_data.is_src_pair) {       // common receiver differential traveltime  ->  common source differential traveltime
            tmp_data.is_src_pair = false;
            tmp_data.is_rec_pair = true;
            tmp_data.id_rec_pair = data_info_nv[i].id_src_pair;
            tmp_data.name_rec_pair = data_info_nv[i].name_src_pair;
            tmp_data.id_src_single = data_info_nv[i].id_rec_single; 
            tmp_data.name_src_single = data_info_nv[i].name_rec_single; 
            tmp_data.cs_dif_travel_time_obs = data_info_nv[i].cr_dif_travel_time_obs;
        }
        data_info_nv[i] = tmp_data;
    } 

    // check new version of src rec data
    if (1 > 0){     // check by Chen Jing
        for(auto iter = src_list_nv.begin(); iter != src_list_nv.end(); iter++){
            std::cout   << "source id: " << iter->second.id 
                        << ", source name: " << iter->second.name
                        << std::endl;
        }

        for(auto iter = rec_list_nv.begin(); iter != rec_list_nv.end(); iter++){
            std::cout   << "receiver id: " << iter->second.id 
                        << ", receiver name: " << iter->second.name
                        << std::endl;
        } 

        for(int i = 0; i < (int)data_info_nv.size(); i++){
            if (data_info_nv[i].is_src_rec){
                std::cout   << "absolute traveltime: " << data_info_nv[i].travel_time_obs
                            << ", source name: " << data_info_nv[i].name_src
                            << ", receiver name: " << data_info_nv[i].name_rec
                            << std::endl; 
            } else if (data_info_nv[i].is_rec_pair){
                std::cout   << "common source differential traveltime: " << data_info_nv[i].cs_dif_travel_time_obs
                            << ", source name: " << data_info_nv[i].name_src_single
                            << ", receiver pair name: " << data_info_nv[i].name_rec_pair[0] 
                            << ", " << data_info_nv[i].name_rec_pair[1]
                            << std::endl; 
            } else if (data_info_nv[i].is_src_pair){
                std::cout   << "common receiver differential traveltime: " << data_info_nv[i].cr_dif_travel_time_obs
                            << ", source pair name: " << data_info_nv[i].name_src_pair[0]
                            << ", " << data_info_nv[i].name_src_pair[1]
                            << ", receiver name: " << data_info_nv[i].name_rec_single 
                            << std::endl; 
            } 
        }
        std::cout << data_info_nv.size() << std::endl;
    }
}

void InputParams::rearrange_data_info(){
    for(int i = 0; i < (int)data_info_nv.size(); i++){
        DataInfo data = data_info_nv[i];
        if(data.is_src_rec){    // add absolute traveltime 
            data_info_smap[data.name_src].push_back(data);
        } else if (data.is_rec_pair){   // add common source differential traveltime 
            data_info_smap[data.name_src_single].push_back(data);
        } else if (data.is_src_pair){   // add common receiver differential traveltime 
            data_info_smap[data.name_src_pair[0]].push_back(data);
            data_info_smap[data.name_src_pair[1]].push_back(data);
        }
    }
}

void InputParams::generate_src_list_prepare(){
    for(auto iter = data_info_smap.begin(); iter != data_info_smap.end(); iter++){
        for (DataInfo data : iter->second){
            if (data.is_src_pair){      // if this source has common source differential traveltime data
                // add this source and turn to the next source
                src_list_prepare_nv[iter->first] = src_list_nv[iter->first];
                break;  
            }
        }
    }
}

void InputParams::generate_syn_time_list(){
    for(int i = 0; i < (int)data_info_nv.size(); i++){
        DataInfo data = data_info_nv[i];
        if(data.is_src_rec){    // add absolute traveltime 
            syn_time_list_sr[data.name_src][data.name_rec] = 0.0;
        } else if (data.is_rec_pair){   // add common source differential traveltime 
            syn_time_list_sr[data.name_src_single][data.name_rec_pair[0]] = 0.0;
            syn_time_list_sr[data.name_src_single][data.name_rec_pair[1]] = 0.0;
        } else if (data.is_src_pair){   // add common receiver differential traveltime 
            syn_time_list_sr[data.name_src_pair[0]][data.name_rec_single] = 0.0;
            syn_time_list_sr[data.name_src_pair[1]][data.name_rec_single] = 0.0;
        }
    }
}

void InputParams::initialize_syn_time_list(){
    for(auto iter1 = syn_time_list_sr.begin(); iter1 != syn_time_list_sr.end(); iter1++){
        for(auto iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++){
            iter2->second = 0.0;
        }
    }
}
void InputParams::reduce_syn_time_list(){
    for(auto iter1 = syn_time_list_sr.begin(); iter1 != syn_time_list_sr.end(); iter1++){
        for(auto iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++){
            allreduce_cr_sim_single_inplace(iter2->second);
        }
    }
}

void InputParams::initialize_adjoint_source(){
    for(auto iter = rec_list_nv.begin(); iter != rec_list_nv.end(); iter++){
        iter->second.adjoint_source = 0;
    }
}

void InputParams::do_swap_src_rec_ver2(){

    // swap src/rec points
    // at this moment, all the sources are divided into src_points (regional) and tele_src_points (teleseismic)

    std::vector<SrcRec> new_srcs; // new src points
    std::vector<std::vector<SrcRec>> new_recs;

    // generate new source list
    for (long unsigned int i_src = 0; i_src < src_points.size(); i_src++) {
        for(long unsigned int i_rec = 0; i_rec < rec_points[i_src].size(); i_rec++) {

            if (new_srcs.size() == 0){
                if (rec_points[i_src][i_rec].is_src_rec == true){       // single source-receiver pair, absolute traveltime
                    new_srcs.push_back(rec_points[i_src][i_rec]);
                    new_srcs.back().id_srcs_ori.push_back(rec_points[i_src][i_rec].id_src);
                } else if (rec_points[i_src][i_rec].is_rec_pair == true) {      // common source  differential traveltime 
                    SrcRec src1;        // src 1
                    SrcRec src2;        // src 2
                    // build src1
                    src1.lat = rec_points[i_src][i_rec].lat_pair[0];
                    src1.lon = rec_points[i_src][i_rec].lon_pair[0];
                    src1.dep = rec_points[i_src][i_rec].dep_pair[0];
                    src1.name_rec = rec_points[i_src][i_rec].name_rec_pair[0];

                    src2.lat = rec_points[i_src][i_rec].lat_pair[1];
                    src2.lon = rec_points[i_src][i_rec].lon_pair[1];
                    src2.dep = rec_points[i_src][i_rec].dep_pair[1];
                    src2.name_rec = rec_points[i_src][i_rec].name_rec_pair[1];
                    
                    new_srcs.push_back(src1);
                    new_srcs.back().id_srcs_ori.push_back(rec_points[i_src][i_rec].id_src);     // CHS: 这个是用记录new_src中每个新源(台站) 对应的旧源(地震)，并生成新的rec
                    new_srcs.push_back(src2);
                    new_srcs.back().id_srcs_ori.push_back(rec_points[i_src][i_rec].id_src);
                }
            } else if (new_srcs.size() != 0) {
                if (rec_points[i_src][i_rec].is_src_rec == true){       // single source-receiver pair, absolute traveltime
                    bool found = false;
                    // check if the existing element has same name
                    for (long unsigned int i_new_src = 0; i_new_src < new_srcs.size(); i_new_src++) {
                        if (new_srcs[i_new_src].name_rec.compare(rec_points[i_src][i_rec].name_rec) == 0) {
                            new_srcs[i_new_src].id_srcs_ori.push_back(rec_points[i_src][i_rec].id_src);     // CHS: 附上关联地震
                            found=true;
                            break;
                        }
                        else {
                        }
                    }
                    if (!found) {
                        new_srcs.push_back(rec_points[i_src][i_rec]);
                        new_srcs.back().id_srcs_ori.push_back(rec_points[i_src][i_rec].id_src);
                    }
                } else if (rec_points[i_src][i_rec].is_rec_pair == true) {      // common source  differential traveltime 
                    for (int i = 0; i <= 1; i++){   // 
                        bool found = false;
                        // check if the existing element has same name
                        for (long unsigned int i_new_src = 0; i_new_src < new_srcs.size(); i_new_src++) {
                            if (new_srcs[i_new_src].name_rec.compare(rec_points[i_src][i_rec].name_rec_pair[i]) == 0) {
                                new_srcs[i_new_src].id_srcs_ori.push_back(rec_points[i_src][i_rec].id_src);     
                                found=true;
                                break;
                            }
                            else {
                            }
                        }
                        if (!found) {
                            SrcRec src1;
                            src1.lat = rec_points[i_src][i_rec].lat_pair[i];
                            src1.lon = rec_points[i_src][i_rec].lon_pair[i];
                            src1.dep = rec_points[i_src][i_rec].dep_pair[i];
                            src1.name_rec = rec_points[i_src][i_rec].name_rec_pair[i];

                            new_srcs.push_back(src1);
                            new_srcs.back().id_srcs_ori.push_back(rec_points[i_src][i_rec].id_src);
                        }
                    }
                }
            }
        }
    }

    // generate new rec list
    for (long unsigned int i_src = 0; i_src < new_srcs.size(); i_src++) {
        std::vector<SrcRec> tmp_list_recs;

        new_recs.push_back(tmp_list_recs);
    }

    for (long unsigned int i_src = 0; i_src < new_srcs.size(); i_src++) {

        for (auto& i_src_ori : new_srcs[i_src].id_srcs_ori){

            for (auto& tmp_rec_ori : rec_points_back[i_src_ori]){   // CHS: 循环旧源对应的数据，如果发现一致的台站(对)名，就将信息提取出来，构造新的rec
                if (tmp_rec_ori.is_src_rec == true){    // CHS: 如果这个数据是 绝对走时数据 
                    if (tmp_rec_ori.name_rec.compare(new_srcs[i_src].name_rec)==0) {
                        // we can use the same arrival time for a src-rec pair by principle of reciprocity (Aki & Richards, 2002).
                        SrcRec tmp_new_rec = src_points_back[i_src_ori];    // CHS: 找到旧源
                        tmp_new_rec.arr_time     = tmp_rec_ori.arr_time;
                        tmp_new_rec.arr_time_ori = tmp_rec_ori.arr_time_ori;
                        tmp_new_rec.id_rec_ori   = tmp_rec_ori.id_rec;      // CHS: 后面输出src rec 文件会用到，现在先不管
                        tmp_new_rec.is_src_rec = true;
                        new_recs[i_src].push_back(tmp_new_rec);
                    }
                } 
                else if (tmp_rec_ori.is_src_rec == true) {  // CHS: 如果这个数据是 common source differential traveltime 数据
                    if (tmp_rec_ori.name_rec.compare(new_srcs[i_src].name_rec_pair[0])==0) {
                        SrcRec tmp_new_rec;
                        tmp_new_rec.dif_arr_time = tmp_rec_ori.dif_arr_time;
                        tmp_new_rec.dif_arr_time_ori = tmp_rec_ori.dif_arr_time_ori;
                        tmp_new_rec.is_src_pair = true;
                        tmp_new_rec.id_rec_pair[0] = src_points_back[i_src_ori].id_src;
                        tmp_new_rec.dep_pair[0] = src_points_back[i_src_ori].dep;   // CHS: pair 的第一个是原震源，第二个是关联的另一个台站
                        tmp_new_rec.lat_pair[0] = src_points_back[i_src_ori].lat;   
                        tmp_new_rec.lon_pair[0] = src_points_back[i_src_ori].lon;   
                        tmp_new_rec.name_rec_pair[0] = src_points_back[i_src_ori].name_src;   

                        tmp_new_rec.id_rec_pair[1] = tmp_rec_ori.id_rec_pair[1];
                        tmp_new_rec.dep_pair[1] = tmp_rec_ori.dep_pair[1];   // CHS: pair 的第一个是原震源，第二个是关联的另一个台站
                        tmp_new_rec.lat_pair[1] = tmp_rec_ori.lat_pair[1];   
                        tmp_new_rec.lon_pair[1] = tmp_rec_ori.lon_pair[1];   
                        tmp_new_rec.name_rec_pair[1] = tmp_rec_ori.name_rec_pair[1];  
                        new_recs[i_src].push_back(tmp_new_rec);
                    } else if (tmp_rec_ori.name_rec.compare(new_srcs[i_src].name_rec_pair[0])==1) {
                        SrcRec tmp_new_rec;
                        tmp_new_rec.dif_arr_time = - tmp_rec_ori.dif_arr_time;
                        tmp_new_rec.dif_arr_time_ori = - tmp_rec_ori.dif_arr_time_ori;
                        tmp_new_rec.is_src_pair = true;
                        tmp_new_rec.id_rec_pair[0] = src_points_back[i_src_ori].id_src;
                        tmp_new_rec.dep_pair[0] = src_points_back[i_src_ori].dep;   // CHS: pair 的第一个是原震源，第二个是关联的另一个台站
                        tmp_new_rec.lat_pair[0] = src_points_back[i_src_ori].lat;   
                        tmp_new_rec.lon_pair[0] = src_points_back[i_src_ori].lon;   
                        tmp_new_rec.name_rec_pair[0] = src_points_back[i_src_ori].name_src;   

                        tmp_new_rec.id_rec_pair[1] = tmp_rec_ori.id_rec_pair[0];
                        tmp_new_rec.dep_pair[1] = tmp_rec_ori.dep_pair[0];   // CHS: pair 的第一个是原震源，第二个是关联的另一个台站
                        tmp_new_rec.lat_pair[1] = tmp_rec_ori.lat_pair[0];   
                        tmp_new_rec.lon_pair[1] = tmp_rec_ori.lon_pair[0];   
                        tmp_new_rec.name_rec_pair[1] = tmp_rec_ori.name_rec_pair[0];  
                        new_recs[i_src].push_back(tmp_new_rec);
                    }

                }
            }

        }

    }


    // backup and set new src/rec points
    // !! ONLY REGIONAL EVENTS AND RECESVERS ARE STORED IN *_back VECTORS !!
    src_points.clear();
    rec_points.clear();
    // teleseismic events are concatenate to the vectors below later
    src_points = new_srcs;
    rec_points = new_recs;

}



void InputParams::prepare_src_list(){
    // if (src_rec_file_exist) {

    //     int n_all_src = src_points.size();
    //     src_ids_this_sim.clear();

    //     // assign elements of src_points
    //     for (int i_src = 0; i_src < n_all_src; i_src++){
    //         if (i_src % n_sims == id_sim){
    //             src_ids_this_sim.push_back(i_src);
    //         }
    //     }

    //     // check IP.src_ids_this_sim for this rank
    //     if (myrank==0) {
    //         std::cout << id_sim << " assigned src id : ";
    //         for (auto& src_id : src_ids_this_sim) {
    //             std::cout << src_id << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    // }

    // new version
    if (src_rec_file_exist) {
        src_ids_this_sim.clear();       
        src_names_this_sim.clear();

        // assign elements of src_points
        int i_src = 0;
        for (auto iter = src_list_nv.begin(); iter != src_list_nv.end(); iter++){
            if (i_src % n_sims == id_sim){
                src_ids_this_sim.push_back(iter->second.id);
                src_names_this_sim.push_back(iter->second.name);
            }
            i_src++;
        }

        // check IP.src_ids_this_sim for this rank
        if (myrank==0) {
            std::cout << id_sim << " assigned src id(name) : ";
            for (int i = 0; i < (int)src_ids_this_sim.size(); i++) {
                std::cout << src_ids_this_sim[i] << "(" << src_names_this_sim[i] << ") ";
            }
            std::cout << std::endl;
        }


        src_ids_this_sim_prepare.clear();       
        src_names_this_sim_prepare.clear();

        // assign elements of src_points
        i_src = 0;
        for (auto iter = src_list_prepare_nv.begin(); iter != src_list_prepare_nv.end(); iter++){
            if (i_src % n_sims == id_sim){
                src_ids_this_sim_prepare.push_back(iter->second.id);
                src_names_this_sim_prepare.push_back(iter->second.name);
            }
            i_src++;
        }

        // check IP.src_ids_this_sim for this rank
        if (myrank==0) {
            std::cout << id_sim << "(prepare synthetic time) assigned src id(name) : ";
            for (int i = 0; i < (int)src_ids_this_sim_prepare.size(); i++) {
                std::cout << src_ids_this_sim_prepare[i] << "(" << src_names_this_sim_prepare[i] << ") ";
            }
            std::cout << std::endl;
        }
    }

}



void InputParams::gather_all_arrival_times_to_main(){
    for (long unsigned int i_src = 0; i_src < src_points.size(); i_src++){
        if (subdom_main && id_subdomain==0){
            // check if the target source is calculated by this simulation group
            if (std::find(src_ids_this_sim.begin(), src_ids_this_sim.end(), i_src) != src_ids_this_sim.end()) {
                // if this sim is main
                if (id_sim == 0) {
                    // do nothing
                } else {
                    // send to main simulation group
                    for (auto &rec : rec_points[i_src]) {
                        send_cr_single_sim(&(rec.arr_time), 0);
                        send_cr_single_sim(&(rec.dif_arr_time), 0);
                    }
                }
            } else {
                if (id_sim == 0) {
                    // receive
                    int id_sim_group = i_src % n_sims;

                    for (auto &rec : rec_points[i_src]) {
                        recv_cr_single_sim(&(rec.arr_time), id_sim_group);
                        recv_cr_single_sim(&(rec.dif_arr_time), id_sim_group);
                    }
                } else {
                    // do nothing
                }
            }
        }
    }
}

void InputParams::gather_all_arrival_times_to_main_nv(){

    for (auto iter = src_list_nv.begin(); iter != src_list_nv.end(); iter++){
        int id_src = iter->second.id;
        std::string name_src = iter->second.name;

        if (subdom_main && id_subdomain==0){
            // check if the target source is calculated by this simulation group
            if (std::find(src_ids_this_sim.begin(), src_ids_this_sim.end(), id_src) != src_ids_this_sim.end()) {
                if (id_sim == 0) {
                    // do nothing
                } else {
                    // send to main simulation group
                    for (auto iter2 = syn_time_list_sr[name_src].begin(); iter2 != syn_time_list_sr[name_src].end(); iter2++){
                        send_cr_single_sim(&(iter2->second), 0);
                    }
                    
                }
            } else {
                if (id_sim == 0) {
                    // receive
                    int id_sim_group = id_src % n_sims;

                    for (auto iter2 = syn_time_list_sr[name_src].begin(); iter2 != syn_time_list_sr[name_src].end(); iter2++){
                        recv_cr_single_sim(&(iter2->second), id_sim_group);
                    }
                } else {
                    // do nothing
                }
            }
        }
    }
}


void InputParams::write_station_correction_file(int i_inv){
    if(is_sta_correction && run_mode == DO_INVERSION) {  // if apply station correction
        station_correction_file_out = output_dir + "/station_correction_file_step_" + int2string_zero_fill(i_inv) +".dat";

        std::ofstream ofs;

        if (world_rank == 0 && subdom_main && id_subdomain==0){    // main processor of subdomain && the first id of subdoumains

            ofs.open(station_correction_file_out);

            ofs << "# stname " << "   lat   " << "   lon   " << "elevation   " << " station correction (s) " << std::endl;
            for(auto iter = station_correction.begin(); iter != station_correction.end(); iter++){
                std::string tmp_sta_name = iter->first;
                SrcRec rec = rec_list[tmp_sta_name];
                ofs << iter->first << " " 
                    << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << rec.lat << " "
                    << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << rec.lon << " "
                    << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << rec.dep * -1000.0 << " "
                    << std::fixed << std::setprecision(6) << std::setw(9) << std::right << std::setfill(' ') << iter->second << " "
                    << std::endl;
            }

            ofs.close();

        }
        synchronize_all_world();
    }
}

void InputParams::write_src_rec_file(int i_inv) {

    // check src and rec:
    // for (int i_proc = 0; i_proc<=world_nprocs; i_proc++){
    //     if (i_proc == world_rank){
    //         std::cout << "check src info" << std::endl;
    //         for(auto& src: src_points){
    //             std::cout << "world_rank: "<< world_rank <<", src name: " << src.name_rec << ", lat: " << src.lat << ", lon:"
    //                     << src.lon << ", dep:" << src.dep << std::endl;
    //         }
    //         std::cout << "check rec info" << std::endl;
    //         for(auto& rec: rec_points){
    //             for (auto &data: rec){
    //                 std::cout << "world_rank: "<< world_rank <<", rec name: " << data.name_src << ", lat: " << data.lat << ", lon:"
    //                     << data.lon << ", dep:" << data.dep << ", arrival time: " << data.arr_time << std::endl;
    //             }
    //         }
    //     }
    //     synchronize_all_world();
    // }


    if (src_rec_file_exist){

        std::ofstream ofs;

        // gather all arrival time info to the main process
        if (n_sims > 1)
            gather_all_arrival_times_to_main();

        // store the calculated travel time to be output
        reverse_src_rec_points();

        if (run_mode == ONLY_FORWARD)
            src_rec_file_out = output_dir + "/src_rec_file_forward.dat";
        else if (run_mode == DO_INVERSION){
            // write out source and receiver points with current inversion iteration number
            src_rec_file_out = output_dir + "/src_rec_file_step_" + int2string_zero_fill(i_inv) +".dat";
        } else if (run_mode == TELESEIS_PREPROCESS) {
            src_rec_file_out = output_dir + "/src_rec_file_teleseis_pre.dat";
        } else if (run_mode == SRC_RELOCATION) {
            src_rec_file_out = output_dir + "/src_rec_file_src_reloc.dat";
        } else {
            std::cerr << "Error: run_mode is not defined" << std::endl;
            exit(1);
        }

        //std::cout << "world_rank: " << world_rank << ", myrank: " << myrank << ", subdom_main: " << subdom_main << ", id_subdomain: " << id_subdomain << std::endl;

        for (long unsigned int i_src = 0; i_src < src_points_out.size(); i_src++){
            if (world_rank == 0 && subdom_main && id_subdomain==0){    // main processor of subdomain && the first id of subdoumains
                if (i_src == 0)
                    ofs.open(src_rec_file_out);
                else
                    ofs.open(src_rec_file_out, std::ios_base::app);

                // set output precision
                // ofs << std::fixed << std::setprecision(ASCII_OUTPUT_PRECISION);

                // format should be the same as input src_rec_file
                // source line :  id_src yearm month day hour min sec lat lon dep_km mag num_recs id_event
                ofs << std::setw(7) << std::right << std::setfill(' ') <<  i_src << " "
                    << src_points_out[i_src].year << " " << src_points_out[i_src].month << " " << src_points_out[i_src].day << " "
                    << src_points_out[i_src].hour << " " << src_points_out[i_src].min   << " "
                    << std::fixed << std::setprecision(2) << std::setw(5) << std::right << std::setfill(' ') << src_points_out[i_src].sec << " "
                    << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << src_points_out[i_src].lat << " "
                    << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << src_points_out[i_src].lon << " "
                    << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << src_points_out[i_src].dep << " "
                    << std::fixed << std::setprecision(2) << std::setw(5) << std::right << std::setfill(' ') << src_points_out[i_src].mag << " "
                    << std::setw(5) << std::right << std::setfill(' ') << src_points_out[i_src].n_data << " "
                    << src_points_out[i_src].name_src << " "
                    << std::fixed << std::setprecision(4) << std::setw(6) << std::right << std::setfill(' ') << src_points_out[i_src].weight
                    << std::endl;
                for (long unsigned int i_rec = 0; i_rec < rec_points_out[i_src].size(); i_rec++){
                    if(!rec_points_out[i_src][i_rec].is_rec_pair){
                        // receiver line : id_src id_rec name_rec lat lon elevation_m phase epicentral_distance_km arival_time
                        ofs << std::setw(7) << std::right << std::setfill(' ') << i_src << " "
                            << std::setw(4) << std::right << std::setfill(' ') << rec_points_out[i_src][i_rec].id_rec << " "
                            << rec_points_out[i_src][i_rec].name_rec << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << rec_points_out[i_src][i_rec].lat << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << rec_points_out[i_src][i_rec].lon << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << -1.0*rec_points_out[i_src][i_rec].dep*1000.0 << " "
                            << rec_points_out[i_src][i_rec].phase << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << rec_points_out[i_src][i_rec].arr_time << " "
                            << std::fixed << std::setprecision(4) << std::setw(6) << std::right << std::setfill(' ') << rec_points_out[i_src][i_rec].weight
                            << std::endl;
                    } else {
                        // receiver pair line : id_src id_rec1 name_rec1 lat1 lon1 elevation_m1 id_rec2 name_rec2 lat2 lon2 elevation_m2 phase differential_arival_time
                        ofs << std::setw(7) << std::right << std::setfill(' ') <<  i_src << " "
                            << std::setw(4) << std::right << std::setfill(' ') << rec_points_out[i_src][i_rec].id_rec_pair[0] << " "
                            << rec_points_out[i_src][i_rec].name_rec_pair[0] << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << rec_points_out[i_src][i_rec].lat_pair[0] << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << rec_points_out[i_src][i_rec].lon_pair[0] << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << -1.0*rec_points_out[i_src][i_rec].dep_pair[0]*1000.0 << " "
                            << std::setw(4) << std::right << std::setfill(' ') << rec_points_out[i_src][i_rec].id_rec_pair[1] << " "
                            << rec_points_out[i_src][i_rec].name_rec_pair[1] << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << rec_points_out[i_src][i_rec].lat_pair[1] << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << rec_points_out[i_src][i_rec].lon_pair[1] << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << -1.0*rec_points_out[i_src][i_rec].dep_pair[1]*1000.0 << " "
                            << rec_points_out[i_src][i_rec].phase << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << rec_points_out[i_src][i_rec].dif_arr_time << " "
                            << std::fixed << std::setprecision(4) << std::setw(6) << std::right << std::setfill(' ') << rec_points_out[i_src][i_rec].weight
                            << std::endl;
                    }
                }

                ofs.close();
            }
            // synchro
            synchronize_all_world();
        }

    }
}

void InputParams::write_src_rec_file_nv(int i_inv) {

    // check src and rec:
    // for (int i_proc = 0; i_proc<=world_nprocs; i_proc++){
    //     if (i_proc == world_rank){
    //         std::cout << "check src info" << std::endl;
    //         for(auto& src: src_points){
    //             std::cout << "world_rank: "<< world_rank <<", src name: " << src.name_rec << ", lat: " << src.lat << ", lon:"
    //                     << src.lon << ", dep:" << src.dep << std::endl;
    //         }
    //         std::cout << "check rec info" << std::endl;
    //         for(auto& rec: rec_points){
    //             for (auto &data: rec){
    //                 std::cout << "world_rank: "<< world_rank <<", rec name: " << data.name_src << ", lat: " << data.lat << ", lon:"
    //                     << data.lon << ", dep:" << data.dep << ", arrival time: " << data.arr_time << std::endl;
    //             }
    //         }
    //     }
    //     synchronize_all_world();
    // }


    if (src_rec_file_exist){

        std::ofstream ofs;

        // gather all arrival time info to the main process
        if (n_sims > 1)
            gather_all_arrival_times_to_main_nv();

        if (run_mode == ONLY_FORWARD)
            src_rec_file_out = output_dir + "/src_rec_file_forward.dat";
        else if (run_mode == DO_INVERSION){
            // write out source and receiver points with current inversion iteration number
            src_rec_file_out = output_dir + "/src_rec_file_step_" + int2string_zero_fill(i_inv) +".dat";
        } else if (run_mode == TELESEIS_PREPROCESS) {
            src_rec_file_out = output_dir + "/src_rec_file_teleseis_pre.dat";
        } else if (run_mode == SRC_RELOCATION) {
            src_rec_file_out = output_dir + "/src_rec_file_src_reloc.dat";
        } else {
            std::cerr << "Error: run_mode is not defined" << std::endl;
            exit(1);
        }

        int data_count = 0;
        for (auto iter = src_list_back_nv.begin(); iter != src_list_back_nv.end(); iter++){
            if (world_rank == 0 && subdom_main && id_subdomain==0){    // main processor of subdomain && the first id of subdoumains
                if (iter == src_list_back_nv.begin())
                    ofs.open(src_rec_file_out);
                else
                    ofs.open(src_rec_file_out, std::ios_base::app);

                std::string name_src = iter->first;
                SrcRecInfo src = iter->second;

                // format should be the same as input src_rec_file
                // source line :  id_src yearm month day hour min sec lat lon dep_km mag num_recs id_event
                ofs << std::setw(7) << std::right << std::setfill(' ') <<  src.id << " "
                    << src.year << " " << src.month << " " << src.day << " "
                    << src.hour << " " << src.min   << " "
                    << std::fixed << std::setprecision(2) << std::setw(5) << std::right << std::setfill(' ') << src.sec << " "
                    << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << src.lat << " "
                    << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << src.lon << " "
                    << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << src.dep << " "
                    << std::fixed << std::setprecision(2) << std::setw(5) << std::right << std::setfill(' ') << src.mag << " "
                    << std::setw(5) << std::right << std::setfill(' ') << src.n_data << " "
                    << src.name << " "
                    << std::fixed << std::setprecision(4) << std::setw(6) << std::right << std::setfill(' ') << 1.0     // the weight of source is assigned to data
                    << std::endl;
                // data line
                for(int i = 0; i < src.n_data; i++){
                    DataInfo data = data_info_back_nv[data_count];
                    if (data.is_src_rec){       // absolute traveltime data
                        std::string name_rec = data.name_rec;
                        SrcRecInfo rec = rec_list_back_nv[name_rec];
                        CUSTOMREAL travel_time;
                        if (get_is_srcrec_swap())     // do swap
                            travel_time = syn_time_list_sr[name_rec][name_src]; 
                        else // undo swap
                            travel_time = syn_time_list_sr[name_src][name_rec]; 

                        // receiver line : id_src id_rec name_rec lat lon elevation_m phase epicentral_distance_km arival_time
                        ofs << std::setw(7) << std::right << std::setfill(' ') << src.id << " "
                            << std::setw(5) << std::right << std::setfill(' ') << rec.id << " "
                            << rec.name << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << rec.lat << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << rec.lon << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << -1.0*rec.dep*1000.0 << " "
                            << data.phase << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << travel_time << " "
                            << std::fixed << std::setprecision(4) << std::setw(6) << std::right << std::setfill(' ') << data.data_weight
                            << std::endl;
                    } else if (data.is_rec_pair){   // common source differential traveltime 
                        std::string name_rec1 = data.name_rec_pair[0];
                        SrcRecInfo rec1 = rec_list_back_nv[name_rec1];
                        std::string name_rec2 = data.name_rec_pair[1];
                        SrcRecInfo rec2 = rec_list_back_nv[name_rec2];
                        CUSTOMREAL cs_dif_travel_time;
                        if (get_is_srcrec_swap())      // do swap
                            cs_dif_travel_time = syn_time_list_sr[name_rec1][name_src] - syn_time_list_sr[name_rec2][name_src]; 
                        else // undo swap
                            cs_dif_travel_time = syn_time_list_sr[name_src][name_rec1] - syn_time_list_sr[name_src][name_rec2]; 

                        // receiver pair line : id_src id_rec1 name_rec1 lat1 lon1 elevation_m1 id_rec2 name_rec2 lat2 lon2 elevation_m2 phase differential_arival_time
                        ofs << std::setw(7) << std::right << std::setfill(' ') <<  src.id << " "
                            << std::setw(5) << std::right << std::setfill(' ') <<  rec1.id << " "
                            << rec1.name << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << rec1.lat << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << rec1.lon << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << -1.0*rec1.dep*1000.0 << " "
                            << std::setw(5) << std::right << std::setfill(' ') << rec2.id << " "
                            << rec2.name << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << rec2.lat << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << rec2.lon << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << -1.0*rec2.dep*1000.0 << " "
                            << data.phase << " "
                            << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << cs_dif_travel_time << " "
                            << std::fixed << std::setprecision(4) << std::setw(6) << std::right << std::setfill(' ') << data.data_weight
                            << std::endl;
                    } else if (data.is_src_pair){   // common receiver differential traveltime 
                        // not ready
                    }

                    data_count++;
                }


                ofs.close();
            }
            // synchro
            synchronize_all_world();
        }

    }
}



void InputParams::reverse_src_rec_points(){

    // loop swapped sources
    for (long unsigned int i_src = 0; i_src < src_points.size(); i_src++){
        // swapped only the regional events && really need swap
        if (src_points[i_src].is_teleseismic == false && swap_src_rec){

            // int id_rec_orig = src_points[i_src].id_rec;
            // loop swapped receivers
            for (long unsigned int i_rec = 0; i_rec < rec_points[i_src].size(); i_rec++){
                int id_src_orig = rec_points[i_src][i_rec].id_src;
                int id_rec_orig = rec_points[i_src][i_rec].id_rec_ori; // cannot fully ecover  the original receiver id

                // store calculated arrival time in backuped receiver list
                rec_points_back[id_src_orig][id_rec_orig].arr_time     = rec_points[i_src][i_rec].arr_time;
                rec_points_back[id_src_orig][id_rec_orig].dif_arr_time = rec_points[i_src][i_rec].dif_arr_time;
                // std::cout   << "world_rank: " << world_rank << ", id_rec_orig: " << id_rec_orig << ", id_src_orig:"
                //             << id_src_orig << ", arr_time: " << rec_points_back[id_src_orig][id_rec_orig].arr_time
                //             << ", i_src:" << i_src << ", i_rec:" << i_rec
                //             << ", rec_points:" << rec_points[i_src][i_rec].arr_time <<std::endl;

                // update relocated source positions
                if (run_mode == SRC_RELOCATION) {
                    src_points_back[id_src_orig].lat = rec_points[i_src][i_rec].lat;
                    src_points_back[id_src_orig].lon = rec_points[i_src][i_rec].lon;
                    src_points_back[id_src_orig].dep = rec_points[i_src][i_rec].dep;
                }
            }
        } else {
            // teleseismic events are not swapped
            for (long unsigned int i_rec = 0; i_rec < rec_points[i_src].size(); i_rec++){
                int id_src_orig = rec_points[i_src][i_rec].id_src;
                // store calculated arrival time in backuped receiver list
                rec_points_back[id_src_orig][i_rec].arr_time = rec_points[i_src][i_rec].arr_time;
                rec_points_back[id_src_orig][i_rec].dif_arr_time = rec_points[i_src][i_rec].dif_arr_time;
            }
        }

    }

    // copy backup to backup to src_points and rec_points
    src_points_out = src_points_back;
    rec_points_out = rec_points_back;


}


void InputParams::separate_region_and_tele_src(){
    // check if the source is inside the simulation boundary

    // initialize vectors for teleseismic events
    tele_src_points.clear();
    tele_rec_points.clear();

    // temporary store src and recs in another vector
    std::vector<SrcRec> src_points_tmp = src_points;
    std::vector<std::vector<SrcRec>> rec_points_tmp = rec_points;

    // clear original src and recs
    src_points.clear();
    rec_points.clear();

    int cc=0;
    for (auto& src_point : src_points_tmp){

        if (src_point.lat < min_lat || src_point.lat > max_lat \
         || src_point.lon < min_lon || src_point.lon > max_lon \
         || src_point.dep < min_dep || src_point.dep > max_dep){

            // for teleseismic events
            src_point.is_teleseismic = true;

            // add this source to the list of tele sources
            tele_src_points.push_back(src_point);
            // add receivers to the list of tele receivers
            tele_rec_points.push_back(rec_points_tmp[cc]);

        } else {

            // for local events
            src_point.is_teleseismic = false;
            // add this source to the list of local sources
            src_points.push_back(src_point);
            // add receivers to the list of local receivers
            rec_points.push_back(rec_points_tmp[cc]);

        }
        cc++;
    }


    // for new version
    tele_src_list_nv.clear();
    src_list_nv.clear();

    for(auto iter = src_list_back_nv.begin(); iter != src_list_back_nv.end(); iter++){
        SrcRecInfo src = iter->second;
        if (src.lat < min_lat || src.lat > max_lat \
         || src.lon < min_lon || src.lon > max_lon \
         || src.dep < min_dep || src.dep > max_dep){

            // out of region (teleseismic events) 
            src.is_out_of_region = true;
            tele_src_list_nv[iter->first] = src;
        } else {
            // within region (teleseismic events) 
            src.is_out_of_region = false;
            src_list_nv[iter->first] = src;
        }
    }

}

void InputParams::merge_region_and_tele_src(){
    if(tele_src_points.size() > 0) {
        src_points.insert(src_points.end(), tele_src_points.begin(), tele_src_points.end());
        rec_points.insert(rec_points.end(), tele_rec_points.begin(), tele_rec_points.end());
    }
}

// check contradictory parameters
void InputParams::check_contradictions(){

    // if run_mode == 0 then the max_iter should be 1
    if (run_mode == ONLY_FORWARD && max_iter_inv > 1){
        std::cout << "Warning: run_mode = 0, max_iter should be 1" << std::endl;
        max_iter_inv = 1;
    }

#ifdef USE_CUDA
    if (use_gpu){

        if (sweep_type != SWEEP_TYPE_LEVEL || n_subprocs != 1){
            if(world_rank == 0) {
                std::cout << "ERROR: In GPU mode, sweep_type must be 1 and n_subprocs must be 1." << std::endl;
                std::cout << "Abort." << std::endl;
            }
            MPI_Finalize();
            exit(1);
        }
    }
#else
    if (use_gpu){
        if(world_rank == 0) {
            std::cout << "ERROR: TOMOATT is not compiled with CUDA." << std::endl;
            std::cout << "Abort." << std::endl;
        }
        MPI_Finalize();
        exit(1);
    }
#endif
}


void InputParams::allocate_memory_tele_boundaries(int np, int nt, int nr, int src_id, \
        bool i_first_in, bool i_last_in, bool j_first_in, bool j_last_in, bool k_first_in) {

    i_first = i_first_in;
    i_last  = i_last_in;
    j_first = j_first_in;
    j_last  = j_last_in;
    k_first = k_first_in;

    // allocate memory for teleseismic boundary sources
    SrcRec& src = src_points[src_id];

    // check if this src is teleseismic source
    if (src.is_teleseismic){
        // North boundary
        if (j_last)
            src.arr_times_bound_N = (CUSTOMREAL*)malloc(sizeof(CUSTOMREAL)*np*nr*N_LAYER_SRC_BOUND);
        // South boundary
        if (j_first)
            src.arr_times_bound_S = (CUSTOMREAL*)malloc(sizeof(CUSTOMREAL)*np*nr*N_LAYER_SRC_BOUND);
        // East boundary
        if (i_last)
            src.arr_times_bound_E = (CUSTOMREAL*)malloc(sizeof(CUSTOMREAL)*nt*nr*N_LAYER_SRC_BOUND);
        // West boundary
        if (i_first)
            src.arr_times_bound_W = (CUSTOMREAL*)malloc(sizeof(CUSTOMREAL)*nt*nr*N_LAYER_SRC_BOUND);
        // Bottom boundary
        if (k_first)
            src.arr_times_bound_Bot = (CUSTOMREAL*)malloc(sizeof(CUSTOMREAL)*nt*np*N_LAYER_SRC_BOUND);

        // boundary source flag
        src.is_bound_src = (bool*)malloc(sizeof(bool)*5);
    }

}

// station correction kernel
void InputParams::station_correction_update(CUSTOMREAL stepsize){
    if (!is_sta_correction)
        return;

    // station correction kernel is generated in the main process and sent the value to all other processors

    // step 1, gather all arrival time info to the main process
    if (n_sims > 1){        
        gather_all_arrival_times_to_main();
    }

    station_correction_value = new CUSTOMREAL[N_receiver];

    // do it in the main processor
    if (id_sim == 0){

        // step 2, initialize the kernel K_{\hat T_i}
        for (auto iter = station_correction_kernel.begin(); iter!=station_correction_kernel.end(); iter++){
            iter->second = 0.0;
        }

        CUSTOMREAL max_kernel = 0.0;

        // step 3, calculate the kernel
        for (long unsigned int i_src = 0; i_src < src_points.size(); i_src++){
            for (auto &rec : rec_points[i_src]) {
                // loop all data
                if(src_points[i_src].is_teleseismic ){
                    station_correction_kernel[rec.name_rec_pair[0]] += _2_CR * rec.weight * src_points[i_src].weight
                        * (rec.dif_arr_time - rec.dif_arr_time_ori + rec.station_correction_pair[0] - rec.station_correction_pair[1]);

                    station_correction_kernel[rec.name_rec_pair[1]] -= _2_CR * rec.weight * src_points[i_src].weight 
                        * (rec.dif_arr_time - rec.dif_arr_time_ori + rec.station_correction_pair[0] - rec.station_correction_pair[1]);
                    max_kernel = std::max(max_kernel,station_correction_kernel[rec.name_rec_pair[0]]);
                    max_kernel = std::max(max_kernel,station_correction_kernel[rec.name_rec_pair[1]]);
                }
            }
        }

        // step 4, update station correction 
        int tmp_count = 0;

        for (auto iter = station_correction_kernel.begin(); iter!=station_correction_kernel.end(); iter++){
            station_correction[iter->first] += iter->second / (-max_kernel) * stepsize;
            station_correction_value[tmp_count] = station_correction[iter->first];
            tmp_count++;
        }

        

    }

    // step 5, broadcast the station correction all all procesors
    broadcast_cr(station_correction_value, N_receiver, 0);
    broadcast_cr_inter_sim(station_correction_value, N_receiver, 0);
    int tmp_count = 0;
    for (auto iter = station_correction.begin(); iter!=station_correction.end(); iter++){
        iter->second = station_correction_value[tmp_count];
        tmp_count++;
    }
    

    // step 6, all processor update the station correction of receiver pair
    for (long unsigned int i_src = 0; i_src < src_points.size(); i_src++){
        for (auto &rec : rec_points[i_src]) {
            // loop all data
            if(src_points[i_src].is_teleseismic ){
                rec.station_correction_pair[0] = station_correction[rec.name_rec_pair[0]];
                rec.station_correction_pair[1] = station_correction[rec.name_rec_pair[1]];
            }
        }
    }

    // std::cout << "final station correction" << std::endl;   
    // auto iter = station_correction.begin();
    // std::cout << "world_rank: " << world_rank << ", id_sim: " << id_sim << ", id_subdomain: " << id_subdomain 
    //           << ", subdom_main" << subdom_main
    //           << ", station name: " << iter->first << ", station correction: " << iter->second << std::endl;   
    // iter = station_correction.end();
    // iter --;
    // std::cout << "world_rank: " << world_rank << ", id_sim: " << id_sim << ", id_subdomain: " << id_subdomain 
    //           << ", subdom_main" << subdom_main
    //           << ", station name: " << iter->first << ", station correction: " << iter->second << std::endl; 
                 

}
