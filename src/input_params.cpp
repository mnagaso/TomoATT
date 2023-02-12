#include "input_params.h"

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

            // output path
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

            // number of max iteration for inversion
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
            if (config["output_setting"]["is_output_in_process"])
                is_output_in_process = config["output_setting"]["is_output_in_process"].as<int>();
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

    stdout_by_rank_zero("parameter file read done.");

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
    broadcast_bool_single(is_output_in_process, 0);
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

    // check contradictory settings
    check_contradictions();

    // broadcast the values to all processes
    stdout_by_rank_zero("read input file successfully.");

}


InputParams::~InputParams(){
    // free memory
    if (subdom_main) {
        for (auto& id_src: src_ids_this_sim){
            SrcRec& src = get_src_point(id_src);
            if(src.is_teleseismic){
                if (j_last && src.arr_times_bound_N != nullptr)
                    free(src.arr_times_bound_N);
                if (j_first && src.arr_times_bound_S != nullptr)
                    free(src.arr_times_bound_S);
                if (i_last && src.arr_times_bound_E != nullptr)
                    free(src.arr_times_bound_E);
                if (i_first && src.arr_times_bound_W != nullptr)
                    free(src.arr_times_bound_W);
                if (k_first && src.arr_times_bound_Bot != nullptr)
                    free(src.arr_times_bound_Bot);

                free(src.is_bound_src);
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
        return depth2radius(get_src_point(id_sim_src).dep);
    else
        return depth2radius(src_dep);
}


CUSTOMREAL InputParams::get_src_lat() {
    if (src_rec_file_exist)
        return get_src_point(id_sim_src).lat*DEG2RAD;
    else
        return src_lat*DEG2RAD;
}


CUSTOMREAL InputParams::get_src_lon() {
    if (src_rec_file_exist)
        return get_src_point(id_sim_src).lon*DEG2RAD;
    else
        return src_lon*DEG2RAD;
}


SrcRec& InputParams::get_src_point(int i_src){

    if (subdom_main){
        for (auto& src: src_points){
            if (src.id_src == i_src)
                return src;
        }

        // if not found, return error
        std::cout << "Error: src point " << i_src << " not found!" << std::endl;
        // assigned src id
        std::cout << "Assigned src ids: ";
        for (auto& src: src_points){
            std::cout << src.id_src << " ";
        }
        std::cout << std::endl;

        exit(1);
    } else {
        // return error because non-subdom_main process should not call this function
        std::cout << "Error: non-subdom_main process should not call this function!" << std::endl;
        exit(1);
    }

}


std::vector<SrcRec>& InputParams::get_rec_points(int id_src) {
    if (subdom_main){
        for (auto& vrecs: rec_points) {
            if (vrecs[0].id_src == id_src)
                return vrecs;
        }

        // if not found, return error
        std::cout << "Error: rec points for src " << id_src << " not found!" << std::endl;
        exit(1);
    } else {
        // return error because non-subdom_main process should not call this function
        std::cout << "Error: non-subdom_main process should not call this function!" << std::endl;
        exit(1);
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

bool InputParams::get_if_src_teleseismic(int id_src) {
    bool if_src_teleseismic;

    if (subdom_main)
        if_src_teleseismic = get_src_point(id_src).is_teleseismic;

    // broadcast to all processes within simultaneous run group
    broadcast_bool_single_sub(if_src_teleseismic, 0);

    return if_src_teleseismic;
}


void InputParams::prepare_src_list(){
    //
    // only the subdom_main process of the first simultaneous run group (id_sim==0 && sim_rank==any && subdom_main) reads src/rec file
    // and stores entile src/rec list in src_points and rec_points
    // then, the subdom_main process of each simultaneous run group (id_sim==any && subdom_main==true) retains only its own src/rec objects,
    // which are actually calculated in those simultaneous run groups
    //

    // read src rec file
    if (src_rec_file_exist && id_sim==0 && subdom_main) {

        parse_src_rec_file(src_rec_file,
                           src_points,
                           rec_points,
                           rec_list,
                           station_correction,
                           station_correction_kernel);

        // total number of unique receivers
        N_receiver = rec_list.size();

        // backup original src/rec list
        src_points_back = src_points;
        rec_points_back = rec_points;

        // check if src positions are within the domain or not (teleseismic source)
        // detected teleseismic source is separated into tele_src_points and tele_rec_points
        separate_region_and_tele_src();

        if (swap_src_rec) {
            // here only reginal events will be processed
            stdout_by_main("Swapping src and rec. This may take few minutes for a large dataset (only regional events will be processed)\n");
            do_swap_src_rec(src_points, rec_points, src_points_back, rec_points_back);
        }

        // concatenate resional and teleseismic src/rec points
        merge_region_and_tele_src();
    }

    // wait
    synchronize_all_world();

    // broadcast src_points and rec_points to all the subdom_main processes of each simultaneous run group
    if (src_rec_file_exist) {
        int n_all_src=0;

        if (id_sim == 0 && subdom_main)
            n_all_src = src_points.size();

        // broadcast n_all_src to all processes
        broadcast_i_single_inter_sim(n_all_src, 0); // id_sim=0&&subdom_main to id_sim=all&&subdom_main
        broadcast_i_single_sub(n_all_src, 0); // id_sim=all&&subdom_main to id_sim=all&&all subprocesses

        // abort program  if n_all_src is zero
        if (n_all_src == 0) {
            if (subdom_main)
                stdout_by_main("ERROR: No source is found in the source file. Aborting...\n");
            exit(1);
        }

        src_ids_this_sim.clear();

        // assign elements of src_points
        for (int i_src = 0; i_src < n_all_src; i_src++){
            // id of simultaneous run group to which the i_src-th element is assigned
            int dst_id_sim = i_src % n_sims;

            if (id_sim==0){
                if (dst_id_sim == id_sim){
                    src_ids_this_sim.push_back(i_src);
                } else if (subdom_main) {
                    // send src_points[i_src] to the main process of dst_id_sim
                    send_src_inter_sim(src_points[i_src], dst_id_sim);
                    // send rec_points[i_src] to the main process of dst_id_sim
                    if (src_points[i_src].n_rec > 0) // #TODO: this is probably not supporting n_rec_pair
                        send_rec_inter_sim(rec_points[i_src], dst_id_sim);
                }
            } else {
                if (dst_id_sim == id_sim){
                    src_ids_this_sim.push_back(i_src);
                    // receive src/rec_points from the main process of dst_id_sim
                    if(subdom_main){
                        // initialize src_points[i_src]
                        src_points.push_back(SrcRec());
                        // receive src_points from the main process of dst_id_sim
                        recv_src_inter_sim(src_points.back(), 0);

                        if (src_points.back().n_rec > 0) {
                            // initialize rec_points
                            rec_points.push_back(std::vector<SrcRec>());
                            for (int i_rec = 0; i_rec < src_points.back().n_rec; i_rec++){
                                rec_points.back().push_back(SrcRec());
                            }
                            // receive rec_points[i_src] from the main process of dst_id_sim
                            recv_rec_inter_sim(rec_points.back(), 0);
                        }
                   }
                } else {
                    // do nothing
                }
            }
        }

        // check IP.src_ids_this_sim for this rank
        //if (myrank==0) {
        //    std::cout << id_sim << " assigned src id : ";
        //    for (auto& src_id : src_ids_this_sim) {
        //        std::cout << src_id << " ";
        //    }
        //    std::cout << std::endl;
        //}
    }

    // wait
    synchronize_all_world();


}


void InputParams::gather_all_arrival_times_to_main(){
    //
    // calculated arrival times at receivers are stored in rec_points of its simultaneous run group
    // this function gathers all arrival times to the subdom_main process of the first simultaneous run group (id_sim==0 && subdom_main==true)
    // from the main processes of each simultaneous run group (id_sim==any && sim_rank==0)
    //

    int n_all_src;

    if (src_rec_file_exist) {

        // broadcast n_all_src to all processes
        if (id_sim == 0 && sim_rank == 0)
            n_all_src = src_points.size();

        broadcast_i_single_inter_sim(n_all_src, 0);

        for (int i_src = 0; i_src < n_all_src; i_src++){

            //if (subdom_main && id_subdomain==0){
            if (sim_rank==0) {
                // check if the target source is calculated by this simulation group
                if (std::find(src_ids_this_sim.begin(), src_ids_this_sim.end(), i_src) != src_ids_this_sim.end()) {
                    // if this sim is main
                    if (id_sim == 0) {
                        // do nothing
                    } else {
                        // send to main simulation group
                        auto& vrecs = get_rec_points(i_src);
                        for (auto &rec : vrecs) {
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

    // check src and rec
    //for (int i_proc = 0; i_proc<=world_nprocs; i_proc++){
    //    if (i_proc == world_rank){
    //        std::cout << "check src info" << std::endl;
    //        for(auto& src: src_points){
    //            std::cout << "world_rank: "<< world_rank <<", src name: " << src.name_rec << ", lat: " << src.lat << ", lon:"
    //                    << src.lon << ", dep:" << src.dep << std::endl;
    //        }
    //        std::cout << "check rec info" << std::endl;
    //        for(auto& rec: rec_points){
    //            for (auto &data: rec){
    //                std::cout << "world_rank: "<< world_rank <<", rec name: " << data.name_src << ", lat: " << data.lat << ", lon:"
    //                    << data.lon << ", dep:" << data.dep << ", arrival time: " << data.arr_time << std::endl;
    //            }
    //        }
    //    }
    //    synchronize_all_world();
    //}


    if (src_rec_file_exist && subdom_main){

        // gather all arrival time info to the main process
        if (n_sims > 1)
            gather_all_arrival_times_to_main();

        // store the calculated travel time to be output
        if (id_sim==0 && sim_rank==0)
            reverse_src_rec_points(src_points, rec_points, \
                                   src_points_back, rec_points_back, \
                                   src_points_out, rec_points_out, \
                                   swap_src_rec, run_mode);

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

        // write out source and receiver points
        writeout_src_rec_file(src_rec_file_out, src_points_out, rec_points_out);
    }
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
            if (src_point.n_rec > 0 || src_point.n_rec_pair > 0)
                tele_rec_points.push_back(rec_points_tmp[cc]);

        } else {

            // for local events
            src_point.is_teleseismic = false;
            // add this source to the list of local sources
            src_points.push_back(src_point);
            // add receivers to the list of local receivers
            if (src_point.n_rec > 0 || src_point.n_rec_pair > 0)
                rec_points.push_back(rec_points_tmp[cc]);
        }
        cc++;
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
    SrcRec& src = get_src_point(src_id);

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

    // free memory
    delete [] station_correction_value;

}
