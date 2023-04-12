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

            // weight of different types of data
            if (config["inversion"]["abs_time_local_weight"]) {
                abs_time_local_weight = config["inversion"]["abs_time_local_weight"].as<CUSTOMREAL>();
            }
            if (config["inversion"]["cr_dif_time_local_weight"]) {
                cr_dif_time_local_weight = config["inversion"]["cr_dif_time_local_weight"].as<CUSTOMREAL>();
            }
            if (config["inversion"]["cs_dif_time_local_weight"]) {
                cs_dif_time_local_weight = config["inversion"]["cs_dif_time_local_weight"].as<CUSTOMREAL>();
            }
            if (config["inversion"]["teleseismic_weight"]) {
                teleseismic_weight = config["inversion"]["teleseismic_weight"].as<CUSTOMREAL>();
            }
            if (config["inversion"]["is_balance_data_weight"]) {
                is_balance_data_weight = config["inversion"]["is_balance_data_weight"].as<int>();
            }



            // ----- for relocation ----

            // step size of relocation
            if (config["inversion"]["step_length_src_reloc"]) {
                step_length_src_reloc = config["inversion"]["step_length_src_reloc"].as<CUSTOMREAL>();
            }
            // step size decay of relocation
            if (config["inversion"]["step_length_decay"]) {
                step_length_decay = config["inversion"]["step_length_decay"].as<CUSTOMREAL>();
            }

            // max change for earthquake location
            if (config["inversion"]["max_change_dep_lat_lon"]) {
                max_change_dep = config["inversion"]["max_change_dep_lat_lon"][0].as<CUSTOMREAL>();
                max_change_lat = config["inversion"]["max_change_dep_lat_lon"][1].as<CUSTOMREAL>();
                max_change_lon = config["inversion"]["max_change_dep_lat_lon"][2].as<CUSTOMREAL>();
            }

            // norm(grad) threshold of stopping relocation
            if (config["inversion"]["tol_gradient"]) {
                TOL_SRC_RELOC = config["inversion"]["tol_gradient"].as<CUSTOMREAL>();
            }

            // local search scheme for relocation
            if (config["inversion"]["is_ortime_local_search"]) {
                is_ortime_local_search = config["inversion"]["is_ortime_local_search"].as<int>();
            }
            // kernel =   K_t/ref_ortime_change (only for is_ortime_local_search : 1)
            if (config["inversion"]["ref_ortime_change"]) {
                ref_ortime_change = config["inversion"]["ref_ortime_change"].as<CUSTOMREAL>();
            }
            // the change of ortime do not exceed max_change (only for is_ortime_local_search : 1)
            if (config["inversion"]["max_change_ortime"]) {
                max_change_ortime = config["inversion"]["max_change_ortime"].as<CUSTOMREAL>();
            }
            // step size of ortime is :  step_length_src_reloc * step_length_ortime_rescale
            if (config["inversion"]["step_length_ortime_rescale"]) {
                step_length_ortime_rescale = config["inversion"]["step_length_ortime_rescale"].as<CUSTOMREAL>();
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
            if (config["output_setting"]["is_single_precision_output"])
                is_single_precision_output = config["output_setting"]["is_single_precision_output"].as<int>();
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

        // write parameter file to output directory
        write_params_to_file();

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

    broadcast_bool_single(src_rec_file_exist, 0);
    broadcast_bool_single(sta_correction_file_exist, 0);

    // This have to be done only after broadcast src_rec_file_exist
    if (src_rec_file_exist == false){
        SrcRecInfo src;
        src.id = 0;
        src.name = "s0";
        src.lat    = src_lat;
        src.lon    = src_lon;
        src.dep    = src_dep;
        src_map[src.name] = src;
        //src_ids_this_sim.push_back(0);
        //src_names_this_sim.push_back("s0");
        SrcRecInfo rec;
        rec.id = 0;
        rec.name = "r0";
        rec_map[rec.name] = rec;
        DataInfo data;
        data_map[src.name][rec.name].push_back(data);
    }
    broadcast_bool_single(swap_src_rec, 0);

    broadcast_str(src_rec_file, 0);
    broadcast_str(sta_correction_file, 0);
    broadcast_str(output_dir, 0);
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

    broadcast_cr_single(abs_time_local_weight, 0);
    broadcast_cr_single(cs_dif_time_local_weight, 0);
    broadcast_cr_single(cr_dif_time_local_weight, 0);
    broadcast_cr_single(teleseismic_weight, 0);
    broadcast_i_single(is_balance_data_weight, 0);

    broadcast_cr_single(max_change_dep, 0);
    broadcast_cr_single(max_change_lat, 0);
    broadcast_cr_single(max_change_lon, 0);
    broadcast_cr_single(TOL_SRC_RELOC, 0);
    broadcast_i_single(is_ortime_local_search, 0);
    broadcast_cr_single(ref_ortime_change, 0);
    broadcast_cr_single(max_change_ortime, 0);
    broadcast_cr_single(step_length_ortime_rescale, 0);

    broadcast_i_single(smooth_method, 0);
    broadcast_cr_single(smooth_lr, 0);
    broadcast_cr_single(smooth_lt, 0);
    broadcast_cr_single(smooth_lp, 0);
    broadcast_i_single(optim_method, 0);
    broadcast_i_single(max_iter_inv, 0);
    broadcast_cr_single(step_size_init, 0);
    broadcast_cr_single(step_size_init_sc, 0);
    broadcast_cr_single(step_size_decay, 0);
    broadcast_cr_single(step_length_src_reloc, 0);
    broadcast_cr_single(step_length_decay, 0);
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
    broadcast_bool_single(is_single_precision_output, 0);
    broadcast_bool_single(is_inv_slowness, 0);
    broadcast_bool_single(is_inv_azi_ani, 0);
    broadcast_bool_single(is_inv_rad_ani, 0);
    broadcast_cr(kernel_taper,2,0);
    broadcast_bool_single(is_sta_correction, 0);

    // check contradictory settings
    check_contradictions();

    // broadcast the values to all processes
    stdout_by_rank_zero("read input file successfully.");

}


InputParams::~InputParams(){
    // free memory
    if (subdom_main) {
        //for (std::string name_src : src_names_this_sim){
        //    SrcRecInfo& src = get_src_point(name_src);

        for (auto iter = src_map.begin(); iter != src_map.end(); iter++){
            SrcRecInfo& src = iter->second;
            if (src.is_out_of_region){
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

    // clear all src, rec, data
    src_map.clear();
    src_map_tele.clear();
    src_map_all.clear();
    src_map_back.clear();
    src_map_tele.clear();
    rec_map.clear();
    rec_map_tele.clear();
    rec_map_all.clear();
    rec_map_back.clear();
    data_map.clear();
    data_map_tele.clear();
    data_map_all.clear();
    data_map_back.clear();
}


void InputParams::write_params_to_file() {
    // write all the simulation parameters in another yaml file
    std::string file_name = "params_log.yaml";
    std::ofstream fout(file_name);
    fout << "version: " << 2 << std::endl;

    fout << std::endl;
    fout << "domain:" << std::endl;
    fout << "   min_max_dep: [" << min_dep << ", " << max_dep << "] # depth in km" << std::endl;
    fout << "   min_max_lat: [" << min_lat << ", " << max_lat << "] # latitude in degree" << std::endl;
    fout << "   min_max_lon: [" << min_lon << ", " << max_lon << "] # longitude in degree" << std::endl;
    fout << "   n_rtp: [" << ngrid_k << ", " << ngrid_j << ", " << ngrid_i << "] # number of nodes in depth,latitude,longitude direction" << std::endl;

    fout << std::endl;
    fout << "source:" << std::endl;
    fout << "   src_rec_file: " << src_rec_file     << " # source receiver file path" << std::endl;
    fout << "   swap_src_rec: " << int(swap_src_rec) << " # swap source and receiver (1: yes, 0: no)" << std::endl;

    fout << std::endl;
    fout << "model:" << std::endl;
    fout << "   init_model_path: " << init_model_path << " # path to initial model file " << std::endl;
    // check if model_1d_name has any characters
    if (model_1d_name.size() > 0)
        fout << "   model_1d_name: " << model_1d_name;
    else
        fout << "#   model_1d_name: " << "dummy_model_1d_name";
    fout << " # 1D model name used in teleseismic 2D solver (iasp91, ak135, user_defined is available), defined in include/1d_model.h" << std::endl;

    fout << std::endl;
    fout << "inversion:" << std::endl;
    fout << "   run_mode: "           << run_mode << " # 0 for forward simulation only, 1 for inversion" << std::endl;
    fout << "   output_dir: "         << output_dir << " # path to output director (default is ./OUTPUT_FILES/)" << std::endl;
    fout << "   optim_method: "          << optim_method << " # optimization method. 0 : grad_descent, 1 : halve-stepping, 2 : lbfgs (EXPERIMENTAL)" << std::endl;
    fout << "   max_iterations_inv: "    << max_iter_inv << " # maximum number of inversion iterations" << std::endl;
    fout << "   step_size: "             << step_size_init << " # initial step size for model update" << std::endl;
    fout << "   step_size_sc: "          << step_size_init_sc << " # ..."  << std::endl;
    fout << "   step_size_decay: "       << step_size_decay << " # ..." << std::endl;
    fout << "   smooth_method: "         << smooth_method << " # 0: multiparametrization, 1: laplacian smoothing (EXPERIMENTAL)" << std::endl;

    fout << std::endl;
    fout << "   # parameters for multiparametric inversion" << std::endl;
    fout << "   n_inversion_grid: "   << n_inversion_grid << " # number of inversion grid sets" << std::endl;
    fout << "   n_inv_dep_lat_lon: [" << n_inv_r << ", " << n_inv_t << ", " << n_inv_p << "] # number of the base inversion grid points" << std::endl;
    if (sta_correction_file_exist)
        fout << "   sta_correction_file: " << sta_correction_file;
    else
        fout << "#   sta_correction_file: " << "dummy_sta_correction_file";
    fout << " # station correction file path" << std::endl;
    fout << "   type_dep_inv: " << type_dep_inv << " # 0: uniform inversion grid, 1: flexible grid" <<std::endl;
    fout << "   type_lat_inv: " << type_lat_inv << " # 0: uniform inversion grid, 1: flexible grid" <<std::endl;
    fout << "   type_lon_inv: " << type_lon_inv << " # 0: uniform inversion grid, 1: flexible grid" <<std::endl;

    fout << std::endl;
    fout << "   # parameters for uniform inversion grid" << std::endl;
    fout << "   min_max_dep_inv: [" << min_dep_inv << ", " << max_dep_inv << "]" << " # depth in km (Radius of the earth is defined in config.h/R_earth)"  << std::endl;
    fout << "   min_max_lat_inv: [" << min_lat_inv << ", " << max_lat_inv << "]" << " # latitude in degree"  << std::endl;
    fout << "   min_max_lon_inv: [" << min_lon_inv << ", " << max_lon_inv << "]" << " # longitude in degree" << std::endl;

    fout << std::endl;
    fout << "   # parameters for flexible inversion grid" << std::endl;
    if (n_inv_r_flex_read) {
        fout << "   n_inv_r_flex: " << n_inv_r_flex << std::endl;
        fout << "   dep_inv: [";
        for (int i = 0; i < n_inv_r_flex; i++){
            fout << dep_inv[i];
            if (i != n_inv_r_flex-1)
                fout << ", ";
        }
        fout << "]" << std::endl;
    } else {
        fout << "#   n_inv_r_flex: " << "3" << std::endl;
        fout << "#   dep_inv: " << "[1, 1, 1]" << std::endl;
    }
    if (n_inv_t_flex_read) {
        fout << "   n_inv_t_flex: " << n_inv_t_flex << std::endl;
        fout << "   lat_inv: [";
        for (int i = 0; i < n_inv_t_flex; i++){
            fout << lat_inv[i];
            if (i != n_inv_t_flex-1)
                fout << ", ";
        }
        fout << "]" << std::endl;
    } else {
        fout << "#   n_inv_t_flex: " << "3" << std::endl;
        fout << "#   lat_inv: " << "[1, 1, 1]" << std::endl;
    }
    if (n_inv_p_flex_read) {
        fout << "   n_inv_p_flex: " << n_inv_p_flex << std::endl;
        fout << "   lon_inv: [";
        for (int i = 0; i < n_inv_p_flex; i++){
            fout << lon_inv[i];
            if (i != n_inv_p_flex-1)
                fout << ", ";
        }
        fout << "]" << std::endl;
    } else {
        fout << "#   n_inv_p_flex: " << "3" << std::endl;
        fout << "#   lon_inv: " << "[1, 1, 1]" << std::endl;
    }

    fout << std::endl;
    fout << "   # parameters for halve-stepping or lbfg mode" << std::endl;
    fout << "   max_sub_iterations: "    << max_sub_iterations << " # maximum number of each sub-iteration" << std::endl;
    fout << "   l_smooth_rtp: ["         << smooth_lr << ", " << smooth_lt << ", " << smooth_lp << "] # smoothing coefficients for laplacian smoothing" << std::endl;
    fout << "   regularization_weight: " << regularization_weight << " # weight value for regularization (lbfgs mode only)" << std::endl;

    fout << std::endl;
    fout << "inv_strategy: # flags for selecting the target parameters to be inversed" << std::endl;
    fout << "   is_inv_slowness: "   << int(is_inv_slowness) << " # 1: slowness value will be calculated in inversion, 0: will not be calculated" << std::endl;
    fout << "   is_inv_azi_ani: "    << int(is_inv_azi_ani)  << " # 1: azimuth anisotropy value will be calculated in inversion, 0: will not be calculated"<< std::endl;
    fout << "   is_inv_rad_ani: "    << int(is_inv_rad_ani)  << " # flag for radial anisotropy (Not implemented yet)" << std::endl;
    fout << "   kernel_taper: ["     << kernel_taper[0] << ", " << kernel_taper[1] << "]" << std::endl;
    fout << "   is_sta_correction: " << int(is_sta_correction) << std::endl;

    fout << std::endl;
    fout << "parallel: # parameters for parallel computation" << std::endl;
    fout << "   n_sims: "    << n_sims << " # number of simultanoues runs" << std::endl;
    fout << "   ndiv_rtp: [" << ndiv_k << ", " << ndiv_j << ", " << ndiv_i << "] # number of subdivision on each direction" << std::endl;
    fout << "   nproc_sub: " << n_subprocs << " # number of processors for sweep parallelization" << std::endl;
    fout << "   use_gpu: "   << int(use_gpu) << " # 1 if use gpu (EXPERIMENTAL)" << std::endl;

    fout << std::endl;
    fout << "calculation:" << std::endl;
    fout << "   convergence_tolerance: " << conv_tol << " # threshold value for checking the convergence for each forward/adjoint run"<< std::endl;
    fout << "   max_iterations: " << max_iter << " # number of maximum iteration for each forward/adjoint run" << std::endl;
    fout << "   stencil_order: " << stencil_order << " # order of stencil, 1 or 3" << std::endl;
    fout << "   stencil_type: " << stencil_type << " # 0: , 1: first-order upwind scheme (only sweep_type 0 is supported) " << std::endl;
    fout << "   sweep_type: " << sweep_type << " # 0: legacy, 1: cuthill-mckee with shm parallelization" << std::endl;
    int ff_flag=0;
    if (output_format == OUTPUT_FORMAT_HDF5) ff_flag = 0;
    else if (output_format == OUTPUT_FORMAT_ASCII) ff_flag = 1;
    else {
        std::cout << "Error: output_format is not defined!" << std::endl;
        exit(1);
    }
    fout << "   output_file_format: " << ff_flag << std::endl;

    fout << std::endl;
    fout << "output_setting:" << std::endl;
    fout << "   is_output_source_field:     " << int(is_output_source_field)         << " # output the calculated field of all sources                            1 for yes; 0 for no;  default: 1" << std::endl;
    fout << "   is_output_model_dat:        " << int(is_output_model_dat)            << " # output model_parameters_inv_0000.dat or not.                          1 for yes; 0 for no;  default: 1" << std::endl;
    fout << "   is_verbose_output:          " << int(is_verbose_output)              << " # output internal parameters, if no, only model parameters are out.     1 for yes; 0 for no;  default: 0" << std::endl;
    fout << "   is_output_final_model:      " << int(is_output_final_model)          << " # output merged final model or not.                                     1 for yes; 0 for no;  default: 1" << std::endl;
    fout << "   is_output_in_process:       " << int(is_output_in_process)           << " # output model at each inv iteration or not.                            1 for yes; 0 for no;  default: 1" << std::endl;
    fout << "   is_single_precision_output: " << int(is_single_precision_output)     << " # output results in single precision or not.                            1 for yes; 0 for no;  default: 0" << std::endl;

    //fout << std::endl;
    //fout << "debug:" << std::endl;
    //fout << "   debug_mode: " << int(if_test) << std::endl;


}


// return radious
CUSTOMREAL InputParams::get_src_radius(const std::string& name_sim_src) {
    if (src_rec_file_exist)
        return depth2radius(get_src_point(name_sim_src).dep);
    else
        return depth2radius(src_dep);
}


CUSTOMREAL InputParams::get_src_lat(const std::string& name_sim_src) {
    if (src_rec_file_exist)
        return get_src_point(name_sim_src).lat*DEG2RAD;
    else
        return src_lat*DEG2RAD;
}


CUSTOMREAL InputParams::get_src_lon(const std::string& name_sim_src) {
    if (src_rec_file_exist)
        return get_src_point(name_sim_src).lon*DEG2RAD;
    else
        return src_lon*DEG2RAD;
}


SrcRecInfo& InputParams::get_src_point(const std::string& name_src){

    if (subdom_main){
        // THIS LOOP MAKE THE SPEED 10x SLOWER
        //for (auto& src: src_map){
        //    if (src.second.name == name_src)
        //        return src.second;
        //}

        //if (src_map.find(name_src) != src_map.end()) // THIS SHOULD ALSO BE AVOIDED
            return src_map[name_src];
        //else {
            // if not found, return error
            std::cout << "Error: src name " << name_src << " not found!" << std::endl;
            // assigned src id
            std::cout << "Assigned src names to this simultanous run : ";
            for (auto& src: src_map){
                std::cout << src.second.name << " ";
            }
            std::cout << std::endl;

            exit(1);
        //}
    } else {
        // return error because non-subdom_main process should not call this function
        std::cout << "Error: non-subdom_main process should not call this function!" << std::endl;
        exit(1);
    }
}


SrcRecInfo& InputParams::get_rec_point(const std::string& name_rec) {
    if (subdom_main){

        // THIS LOOP MAKE THE SPEED 10x SLOWER
        //for (auto& rec: rec_map) {
        //    if (rec.second.name == name_rec)
        //        return rec.second;
        //}

        // check if rec_map[name_rec] exists
        //if (rec_map.find(name_rec) != rec_map.end()) // THIS SHOULD ALSO BE AVOIDED
            return rec_map[name_rec];
        //else {
        //    // if not found, return error
        //    std::cout << "Error: rec name " << name_rec << " not found!" << std::endl;
        //    exit(1);
        //}


        // if not found, return error
        std::cout << "Error: rec name " << name_rec << " not found!" << std::endl;
        exit(1);
    } else {
        // return error because non-subdom_main process should not call this function
        std::cout << "Error: non-subdom_main process should not call this function!" << std::endl;
        exit(1);
   }
}


std::vector<std::string> InputParams::get_rec_names(const std::string& name_src){
    std::vector<std::string> rec_names;

    for (auto iter = data_map[name_src].begin(); iter != data_map[name_src].end(); iter++) {
        rec_names.push_back(iter->first);
    }
    return rec_names;
}



bool InputParams::get_if_src_teleseismic(const std::string& src_name) {
    bool if_src_teleseismic;

    if (subdom_main)
        if_src_teleseismic = get_src_point(src_name).is_out_of_region;
    else
        if_src_teleseismic = false; // need to be broadcasted from subdom_main later

    // broadcast to all processes within simultaneous run group
    broadcast_bool_single_sub(if_src_teleseismic, 0);

    return if_src_teleseismic;
}


void InputParams::prepare_src_map(){
    //
    // only the subdom_main process of the first simultaneous run group (id_sim==0 && sim_rank==any && subdom_main) reads src/rec file
    // and stores entile src/rec list in src_points and rec_points
    // then, the subdom_main process of each simultaneous run group (id_sim==any && subdom_main==true) retains only its own src/rec objects,
    // which are actually calculated in those simultaneous run groups
    //

    // read src rec file
    if (src_rec_file_exist && id_sim==0 && subdom_main) {

        parse_src_rec_file(src_rec_file, \
                           src_map_all, \
                           rec_map_all, \
                           data_map_all, \
                           src_id2name_all, \
                           rec_id2name_back);

        // read station correction file by all processes
        if (sta_correction_file_exist) {
            // store all src/rec info
            parse_sta_correction_file(sta_correction_file,
                                      rec_map_all);
        }

        // copy backups (KEEPED AS THE STATE BEFORE SWAPPING SRC AND REC)
        src_map_back     = src_map_all;
        rec_map_back     = rec_map_all;
        data_map_back    = data_map_all;
        src_id2name_back = src_id2name_all;

        // check if src positions are within the domain or not (teleseismic source)
        // detected teleseismic source is separated into tele_src_points and tele_rec_points
        std::cout << "separate regional and teleseismic src/rec points" << std::endl;
        separate_region_and_tele_src_rec_data(src_map_back, rec_map_back, data_map_back,
                                              src_map_all,  rec_map_all,  data_map_all,
                                              src_map_tele, rec_map_tele, data_map_tele,
                                              data_type,
                                              N_abs_local_data,
                                              N_cr_dif_local_data,
                                              N_cs_dif_local_data,
                                              N_teleseismic_data,
                                              N_data,
                                              min_lat, max_lat, min_lon, max_lon, min_dep, max_dep);

        if (swap_src_rec) {
            // here only reginal events will be processed
            stdout_by_main("Swapping src and rec. This may take few minutes for a large dataset (only regional events will be processed)\n");
            do_swap_src_rec(src_map_all, rec_map_all, data_map_all, src_id2name_all);
        }

        // concatenate resional and teleseismic src/rec points
        //
        // src_map  = src_map_all  + src_map_tele
        // rec_map  = rec_map_all  + rec_map_tele
        // data_map = data_map_all + data_map_tele
        merge_region_and_tele_src(src_map_all,  rec_map_all,  data_map_all,
                                  src_map_tele, rec_map_tele, data_map_tele);

        // abort if number of src_points are less than n_sims
        int n_src_points = src_map_all.size();
        if (n_src_points < n_sims){
            std::cout << "Error: number of sources in src_rec_file is less than n_sims. Abort.1" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

    } // end of if (src_rec_file_exist && id_sim==0 && subdom_main)

    // wait
    synchronize_all_world();


    // to all the subdom_main processes of each simultaneous run group
    // #TODO: not updated yet for new srcrec !!!
    if (src_rec_file_exist) {
        // # TODO: check if this can be placed

        if (world_rank==0)
            std::cout << "\nsource assign to processors\n" <<std::endl;

        // divide and distribute the data below to each simultaneous run group:
        //  src_map,
        //  rec_map,
        //  data_map,
        //  src_id2name,
        distribute_src_rec_data(src_map_all,
                                rec_map_all,
                                data_map_all,
                                src_id2name_all,
                                src_map,
                                rec_map,
                                data_map,
                                src_id2name);

        // now src_id2name_all  includes all src names of after swapping src and rec
        //     src_id2name      includes only src names of this simultaneous run group
        //     src_id2name_back includes only src names of this simultaneous run group before swapping src and rec

        // create source list for common source double difference traveltime
        generate_src_map_with_common_source(data_map, src_map_comm_src, src_id2name_comm_src);

        // prepare source list for teleseismic source
        prepare_src_map_for_2d_solver(src_map_tele, src_id2name_2d, src_map_2d);

        synchronize_all_world();

        if (world_rank==0)
            std::cout << "end parse src_rec file" << std::endl;
    } // end of if src_rec_file_exists
}


// generate a list of events which involve common source double difference traveltime
void InputParams::generate_src_map_with_common_source(std::map<std::string, std::map<std::string, std::vector<DataInfo>>>& data_map_tmp,
                                                      std::map<std::string, SrcRecInfo>& src_map_comm_src_tmp,
                                                      std::vector<std::string>& src_id2name_comm_src_tmp){
    // for earthquake having common receiver differential traveltime, the synthetic traveltime should be computed first at each iteration
    for(auto iter = data_map_tmp.begin(); iter != data_map_tmp.end(); iter++){
        for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++){
            for (auto& data: iter2->second){
                if (data.is_src_pair){
                    // add this source and turn to the next source
                    src_map_comm_src_tmp[iter->first] = src_map[iter->first];
                    // add this source to the list of sources that will be looped in each iteration
                    // if the source is not in the list, the synthetic traveltime will not be computed
                    if (std::find(src_id2name_comm_src_tmp.begin(), src_id2name_comm_src_tmp.end(), iter->first) == src_id2name_comm_src_tmp.end())
                        src_id2name_comm_src_tmp.push_back(iter->first);

                    break;
                }
            }
        }
    }

    // at this line, only subdom_main process of each simultaneous run group has src_map_comm_src_tmp
    // and src_id2name_comm_src_tmp
    // the following lines are to broadcast src_map_comm_src_tmp and src_id2name_comm_src_tmp to all processes
    // of each simultaneous run group
    //if (id_sim != 0 || !subdom_main)
    //    src_id2name_comm_src_tmp.clear();
    // broadcast the size of src_id2name_comm_src_tmp
    int n_src_id2name_comm_src_tmp = src_id2name_comm_src_tmp.size();
    broadcast_i_single_sub(n_src_id2name_comm_src_tmp, 0); // inter-sim

    for (int i = 0; i < n_src_id2name_comm_src_tmp; i++){
        // broadcast the source name
        std::string src_name;
        if (subdom_main){
            src_name = src_id2name_comm_src_tmp[i];
        }

        broadcast_str_sub(src_name, 0);

        if (!subdom_main)
            src_id2name_comm_src_tmp.push_back(src_name);

    }

    // check if this sim group has common source double difference traveltime
    if (src_map_comm_src_tmp.size() > 0){
        src_pair_exists = true;
    }

    // flag if any src_pair exists
    allreduce_bool_inplace_inter_sim(&src_pair_exists, 1); // inter-sim
    allreduce_bool_inplace(&src_pair_exists, 1); // intra-sim
    allreduce_bool_inplace_sub(&src_pair_exists, 1); // intra-subdom

}

void InputParams::initialize_adjoint_source(){
    for(auto iter = rec_map.begin(); iter != rec_map.end(); iter++){
        iter->second.adjoint_source = 0;
    }
}

void InputParams::set_adjoint_source(std::string name_rec, CUSTOMREAL adjoint_source){
    if (rec_map.find(name_rec) != rec_map.end()){
        rec_map[name_rec].adjoint_source = adjoint_source;
    } else {
        std::cout << "error !!!, undefined receiver name when adding adjoint source: " << name_rec << std::endl;
    }
}


// gather all arrival times to main simultaneous run group
// then store them in data_map_all
void InputParams::gather_all_arrival_times_to_main(){

    for (int id_src = 0; id_src < nsrc_total; id_src++){

        // id of simulation group for this source
        int id_sim_group = select_id_sim_for_src(id_src, n_sims);

        // broadcast source name
        std::string name_src;
        if (subdom_main && id_subdomain==0 && id_sim==0){
            name_src = src_id2name_all[id_src];
        }
        broadcast_str_inter_and_intra_sim(name_src, 0);

        if (subdom_main && id_subdomain==0){

            if (id_sim_group==0) { // the simultaneous run group == 0 is the main simultaneous run group
                if (id_sim == 0) {
                    // copy arrival time to data_info_back
                    for (auto iter = data_map[name_src].begin(); iter != data_map[name_src].end(); iter++){
                        for(int i_data = 0; i_data < (int)iter->second.size(); i_data++){
                            // store travel time in all datainfo element of each src-rec pair
                            data_map_all[name_src][iter->first].at(i_data).travel_time = iter->second.at(i_data).travel_time;
                        }
                    }
                } else {
                    // do nothing
                }
            } else { // this source is calculated in the other simultaneous run group than the main

                // the main group receives the arrival time from the other simultaneous run group
                if (id_sim == 0) {
                    // number of receivers
                    int nrec;
                    recv_i_single_sim(&nrec, id_sim_group);

                    // loop over receivers
                    for (int id_rec = 0; id_rec < nrec; id_rec++){
                        // receive name of receiver
                        std::string name_rec;
                        recv_str_sim(name_rec, id_sim_group); ///////

                        // receive number of data
                        int ndata;
                        recv_i_single_sim(&ndata, id_sim_group);

                        // exit if the number of data is not the same with data_map_all
                        if (ndata != (int)data_map_all[name_src][name_rec].size()){
                            std::cout << "ERROR: the number of data calculated sub simultaneous group is not the same with data_map_all" << std::endl;
                            std::cout << "name_src = " << name_src << ", name_rec = " << name_rec << std::endl;
                            std::cout << "ndata = " << ndata << ", data_map_all[name_src][name_rec].size() = " << data_map_all[name_src][name_rec].size() << std::endl;
                            exit(1);   // for rec_3 src_0, ndata = 1   data_map_all[][].size() = 3
                        }

                        // then receive travel time
                        for (auto& data: data_map_all[name_src][name_rec])
                            recv_cr_single_sim(&(data.travel_time), id_sim_group);
                    }

                // the non-main simultaneous run group sends the arrival time to the main group
                } else if (id_sim == id_sim_group) {
                    // send number of receivers
                    int nrec = data_map[name_src].size();
                    send_i_single_sim(&nrec, 0);

                    for (auto iter = data_map[name_src].begin(); iter != data_map[name_src].end(); iter++){
                        // send name of receiver
                        send_str_sim(iter->first, 0);

                        // send number of data
                        int ndata = iter->second.size();
                        send_i_single_sim(&ndata, 0);

                        // then send travel time
                        for (auto& data: iter->second)
                            send_cr_single_sim(&(data.travel_time), 0);
                    }
                } else {
                    // do nothing
                }
            }

        } // end if (subdom_main && id_subdomain==0)

    } // end for  id_src

}


// gther tau_opt to main simultaneous run group
void InputParams::gather_rec_info_to_main(){

    // broadcast total number of recenver to all procs
    int nrec_total= rec_map_all.size();
    broadcast_i_single_inter_and_intra_sim(nrec_total, 0);

    std::vector<std::string> name_rec_all;

    if (subdom_main){
        if (id_sim==0){
            // assigne tau_opt to rec_map_all from its own rec_map
            for (auto iter = rec_map.begin(); iter != rec_map.end(); iter++){
                rec_map_all[iter->first].tau_opt = iter->second.tau_opt;
            }

            // make a list of receiver names
            for (auto iter = rec_map_all.begin(); iter != rec_map_all.end(); iter++){
                name_rec_all.push_back(iter->first);
            }
        }

        for (int irec =  0; irec < nrec_total; irec++){

            // broadcast name of receiver
            std::string name_rec;
            if (id_sim==0){
                name_rec = name_rec_all[irec];
            }
            broadcast_str_inter_sim(name_rec, 0);

            int        rec_counter = 0;
            CUSTOMREAL tau_tmp=0.0;
            // copy value if rec_map[name_rec] exists
            if (rec_map.find(name_rec) != rec_map.end()){
                tau_tmp = rec_map[name_rec].tau_opt;
                rec_counter = 1;
            }

            // reduce counter and tau_tmp
            allreduce_rec_map_var(rec_counter);
            allreduce_rec_map_var(tau_tmp);

            // assign tau_opt to rec_map_all
            if (rec_counter > 0){
                rec_map_all[name_rec].tau_opt = tau_tmp / (CUSTOMREAL)rec_counter;
            }

       }

    } // end of subdom_main


}

// gather traveltimes and calculate differences of synthetic data
void InputParams::gather_traveltimes_and_calc_syn_diff(){

    if (!src_pair_exists) return; // nothing to share

    // gather all synthetic traveltimes to main simultaneous run group
    gather_all_arrival_times_to_main();

    synchronize_all_world();

    if(subdom_main) {

        int mpi_tag_send=9999;
        int mpi_tag_end=9998;

        // main process calculates differences of synthetic data and send them to other processes
        if (id_sim==0){
            int n_total_src_pair = 0;

            // calculate differences of synthetic data
            for (auto iter = data_map_all.begin(); iter != data_map_all.end(); iter++){
                for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++){
                    for (auto& data: iter2->second){
                        if (data.is_src_pair){
                            data.cr_dif_travel_time = data_map_all[data.name_src_pair[0]][data.name_rec].at(0).travel_time \
                                                    - data_map_all[data.name_src_pair[1]][data.name_rec].at(0).travel_time;

                            n_total_src_pair++;

                       }
                    }
                }
            }

            // send differences of synthetic data to other processes
            for (int id_src = 0; id_src < nsrc_total; id_src++){ // MNMN: looping over all of data.id_src_pair[0]

                // id of simulation group for this source
                int id_sim_group = select_id_sim_for_src(id_src, n_sims);

                std::string name_src1 = src_id2name_all[id_src];

                // iterate over receivers
                for (auto iter = data_map_all[name_src1].begin(); iter != data_map_all[name_src1].end(); iter++){
                    std::string name_rec = iter->first;

                    // iterate over data
                    for (int i_data = 0; i_data < (int)iter->second.size(); i_data++){
                        auto& data = iter->second.at(i_data);

                        if (data.is_src_pair){
                            std::string name_src2 = data.name_src_pair[1];

                            if (id_sim_group == 0) {
                                // this source is calculated in the main simultaneous run group
                                set_cr_dif_to_src_pair(data_map[name_src1][name_rec], name_src2, data.cr_dif_travel_time);
                            } else {
                                // send signal with dummy int
                                int dummy = 0;
                                MPI_Send(&dummy, 1, MPI_INT, id_sim_group, mpi_tag_send, inter_sim_comm);

                                // send name_src1
                                send_str_sim(name_src1, id_sim_group);
                                // send name_src2
                                send_str_sim(name_src2, id_sim_group);
                                // send name_rec
                                send_str_sim(name_rec, id_sim_group);
                                // send index of data
                                send_i_single_sim(&i_data, id_sim_group);
                                // send travel time difference
                                send_cr_single_sim(&(data.cr_dif_travel_time), id_sim_group);
                            }
                        }
                    }
                }
            } // end for id_src

            // send end signal (start from 1 because 0 is main process)
            for (int id_sim = 1; id_sim < n_sims; id_sim++){
                // send dummy integer
                int dummy = 0;
                MPI_Send(&dummy, 1, MPI_INT, id_sim, mpi_tag_end, inter_sim_comm);
            }

        // un-main process receives differences of synthetic data from main process
        } else if (id_sim!=0) {
            while (true) {

                // wait with mpi probe
                MPI_Status status;
                MPI_Probe(0, MPI_ANY_TAG, inter_sim_comm, &status);

                // receive signal with dummy int
                int dummy = 0;
                MPI_Recv(&dummy, 1, MPI_INT, 0, MPI_ANY_TAG, inter_sim_comm, &status);

                // if this signal is for sending data
                if (status.MPI_TAG == mpi_tag_send) {

                    std::string name_src1, name_src2, name_rec;

                    // receive src_name1
                    recv_str_sim(name_src1, 0);
                    // receive src_name2
                    recv_str_sim(name_src2, 0);
                    // receive rec_name
                    recv_str_sim(name_rec, 0);
                    // receive index of data
                    int i_data = 0;
                    recv_i_single_sim(&i_data, 0);
                    // receive travel time difference
                    CUSTOMREAL tmp_ttd = 0;
                    recv_cr_single_sim(&(tmp_ttd), 0);
                    set_cr_dif_to_src_pair(data_map[name_src1][name_rec], name_src2, tmp_ttd);


                // if this signal is for terminating the wait loop
                } else if (status.MPI_TAG == mpi_tag_end) {
                    break;
                }

            }
        }
    }

   synchronize_all_world();

}



void InputParams::write_station_correction_file(int i_inv){
    if(is_sta_correction && run_mode == DO_INVERSION) {  // if apply station correction
        station_correction_file_out = output_dir + "/station_correction_file_step_" + int2string_zero_fill(i_inv) +".dat";

        std::ofstream ofs;

        if (world_rank == 0 && subdom_main && id_subdomain==0){    // main processor of subdomain && the first id of subdoumains

            ofs.open(station_correction_file_out);

            ofs << "# stname " << "   lat   " << "   lon   " << "elevation   " << " station correction (s) " << std::endl;
            for(auto iter = rec_map_back.begin(); iter != rec_map_back.end(); iter++){
                SrcRecInfo  rec      = iter->second;
                std::string name_rec = rec.name;

                ofs << rec.name << " "
                    << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << rec.lat << " "
                    << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << rec.lon << " "
                    << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << rec.dep * -1000.0 << " "
                    << std::fixed << std::setprecision(6) << std::setw(9) << std::right << std::setfill(' ') << rec.sta_correct << " "
                    << std::endl;
            }

            ofs.close();

        }
        synchronize_all_world();
    }
}

void InputParams::write_src_rec_file(int i_inv) {

    if (src_rec_file_exist){

        std::ofstream ofs;

        // gather all arrival time info to the main process (need to call even n_sim=1)
        gather_all_arrival_times_to_main();

        // gather tau_opt info
        if (run_mode == SRC_RELOCATION)
            gather_rec_info_to_main();

        // write only by the main processor of subdomain && the first id of subdoumains
        if (world_rank == 0 && subdom_main && id_subdomain==0){

            if (run_mode == ONLY_FORWARD)
                src_rec_file_out = output_dir + "/src_rec_file_forward.dat";
            else if (run_mode == DO_INVERSION){
                // write out source and receiver points with current inversion iteration number
                src_rec_file_out = output_dir + "/src_rec_file_step_" + int2string_zero_fill(i_inv) +".dat";
            } else if (run_mode == TELESEIS_PREPROCESS) {
                src_rec_file_out = output_dir + "/src_rec_file_teleseis_pre.dat";
            } else if (run_mode == SRC_RELOCATION) {
                src_rec_file_out = output_dir + "/src_rec_file_src_reloc_syn.dat";
            } else {
                std::cerr << "Error: run_mode is not defined" << std::endl;
                exit(1);
            }

            // open file
            ofs.open(src_rec_file_out);

            for (int i_src = 0; i_src < (int)src_id2name_back.size(); i_src++){

                std::string name_src = src_id2name_back[i_src];
                SrcRecInfo  src      = src_map_back[name_src];

                // format should be the same as input src_rec_file
                // source line :  id_src year month day hour min sec lat lon dep_km mag num_recs id_event
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
                for (auto& name_rec: rec_id2name_back[i_src]){
                    std::vector<DataInfo> v_data;

                    // CALCULATED DATA IS STORED IN data_map_all

                    v_data = data_map_back[name_src][name_rec];

                    for (const auto& data_ori : v_data){

                        // data type flag
                        bool src_rec_data  = data_ori.is_src_rec;
                        bool src_pair_data = data_ori.is_src_pair;
                        bool rec_pair_data = data_ori.is_rec_pair;

                        DataInfo& data = const_cast<DataInfo&>(data_ori); // dummy copy

                        // absolute traveltime data
                        if (src_rec_data){
                            if (get_is_srcrec_swap()) // reverse swap src and rec
                                data = get_data_src_rec(data_map_all[name_rec][name_src]);
                            else // do not swap
                                data = get_data_src_rec(data_map_all[name_src][name_rec]);

                            SrcRecInfo  rec         = rec_map_back[name_rec];
                            CUSTOMREAL  travel_time = data.travel_time;

                            // receiver line : id_src id_rec name_rec lat lon elevation_m phase epicentral_distance_km arival_time
                            ofs << std::setw(7) << std::right << std::setfill(' ') << i_src << " "
                                << std::setw(5) << std::right << std::setfill(' ') << rec.id << " "
                                << rec.name << " "
                                << std::fixed << std::setprecision(4) << std::setw(9) << rec.lat << " "
                                << std::fixed << std::setprecision(4) << std::setw(9) << rec.lon << " "
                                << std::fixed << std::setprecision(4) << std::setw(9) << -1.0*rec.dep*1000.0 << " "
                                << data.phase << " "
                                << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << travel_time << " "
                                << std::fixed << std::setprecision(4) << std::setw(6) << std::right << std::setfill(' ') << data.data_weight
                                << std::endl;

                        // common source differential traveltime
                        } else if (rec_pair_data){
                            if (get_is_srcrec_swap()) // reverse swap src and rec
                                data = get_data_src_pair(data_map_all[name_rec][name_src]);
                            else // do not swap
                                data = get_data_rec_pair(data_map_all[name_src][name_rec]);

                            std::string name_rec1, name_rec2;
                            CUSTOMREAL  cs_dif_travel_time;

                            if (get_is_srcrec_swap()) { // do reverse swap
                                name_rec1 = data.name_src_pair[0];
                                name_rec2 = data.name_src_pair[1];
                                cs_dif_travel_time = data.cr_dif_travel_time;
                            } else { // do not swap
                                cs_dif_travel_time = data.cs_dif_travel_time;
                            }

                            SrcRecInfo& rec1 = rec_map_back[name_rec1];
                            SrcRecInfo& rec2 = rec_map_back[name_rec2];

                            // receiver pair line : id_src id_rec1 name_rec1 lat1 lon1 elevation_m1 id_rec2 name_rec2 lat2 lon2 elevation_m2 phase differential_arival_time
                            ofs << std::setw(7) << std::right << std::setfill(' ') <<  i_src << " "
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

                        // common receiver differential traveltime
                        } else if (src_pair_data){
                            // TODO: implement this later
                        }

                    } // end of for (const auto& data : v_data)
                } // end of for (auto iter = data_map_back[name_src].begin(); iter != data_map_back[name_src].end(); iter++)

            } // end of for (int i_src = 0; i_src < (int)src_name_list.size(); i_src++)

            // close file
            ofs.close();

            // only for source relocation, output relocated observational data for tomography
            if (run_mode == SRC_RELOCATION) {
                src_rec_file_out = output_dir + "/src_rec_file_src_reloc_obs.dat";

                // open file
                ofs.open(src_rec_file_out);

                for (int i_src = 0; i_src < (int)src_id2name_back.size(); i_src++){

                    std::string name_src = src_id2name_back[i_src];
                    SrcRecInfo  src      = src_map_back[name_src];

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
                    for (auto& name_rec: rec_id2name_back[i_src]){

                        std::vector<DataInfo> v_data;

                        v_data = data_map_back[name_src][name_rec];

                        for (const auto& data_ori : v_data){

                            // data type flag
                            bool src_rec_data  = data_ori.is_src_rec;
                            bool src_pair_data = data_ori.is_src_pair;
                            bool rec_pair_data = data_ori.is_rec_pair;

                            DataInfo& data = const_cast<DataInfo&>(data_ori); // dummy copy

                            // absolute traveltime data
                            if (src_rec_data){

                                if (get_is_srcrec_swap()) // reverse swap src and rec
                                    data = get_data_src_rec(data_map_all[name_rec][name_src]);
                                else {// retern error
                                    std::cerr << "Error: src_rec_data should not be used in src relocation mode!" << std::endl;
                                    exit(1);
                                }

                                SrcRecInfo rec = rec_map_back[name_rec];
                                CUSTOMREAL travel_time_obs = data.travel_time_obs - rec_map_all[name_src].tau_opt;

                                //std::cout << "src_rec_data: " << name_src << " " << name_rec << " " << data.travel_time_obs << " " << rec_map_all[name_src].tau_opt << " " << travel_time_obs << std::endl;

                                // receiver line : id_src id_rec name_rec lat lon elevation_m phase epicentral_distance_km arival_time
                                ofs << std::setw(7) << std::right << std::setfill(' ') << src.id << " "
                                    << std::setw(5) << std::right << std::setfill(' ') << rec.id << " "
                                    << rec.name << " "
                                    << std::fixed << std::setprecision(4) << std::setw(9) << rec.lat << " "
                                    << std::fixed << std::setprecision(4) << std::setw(9) << rec.lon << " "
                                    << std::fixed << std::setprecision(4) << std::setw(9) << -1.0*rec.dep*1000.0 << " "
                                    << data.phase << " "
                                    << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ') << travel_time_obs << " "
                                    << std::fixed << std::setprecision(4) << std::setw(6) << std::right << std::setfill(' ') << data.data_weight
                                    << std::endl;

                            // common source differential traveltime
                            } else if (rec_pair_data){
                                if (get_is_srcrec_swap()) // reverse swap src and rec
                                    data = get_data_src_pair(data_map_all[name_rec][name_src]);
                                else // do not swap
                                    data = get_data_rec_pair(data_map_all[name_src][name_rec]);

                                std::string name_rec1, name_rec2;
                                CUSTOMREAL  cs_dif_travel_time;

                                if (get_is_srcrec_swap()) { // do reverse swap
                                    name_rec1 = data.name_src_pair[0];
                                    name_rec2 = data.name_src_pair[1];

                                    cs_dif_travel_time = data.cr_dif_travel_time;
                                } else { // do not swap
                                    name_rec1 = data.name_rec_pair[0];
                                    name_rec2 = data.name_rec_pair[1];
                                    cs_dif_travel_time = data.cs_dif_travel_time;
                                }

                                SrcRecInfo& rec1 = rec_map_back[name_rec1];
                                SrcRecInfo& rec2 = rec_map_back[name_rec2];

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
                            } else if (src_pair_data){   // common receiver differential traveltime
                                // not ready
                            }

                        }
                    } // end of data loop

                } // end of src loop

                // close file
                ofs.close();

            } // end of run_mode == SRC_RELOCATION
        } // end if (world_rank == 0)


    } // end of src_rec_file_exist

    // wait till the main process finishes to write the output file
    synchronize_all_world();

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


void InputParams::allocate_memory_tele_boundaries(int np, int nt, int nr, std::string name_src, \
        bool i_first_in, bool i_last_in, bool j_first_in, bool j_last_in, bool k_first_in) {

    i_first = i_first_in;
    i_last  = i_last_in;
    j_first = j_first_in;
    j_last  = j_last_in;
    k_first = k_first_in;

    // allocate memory for teleseismic boundary sources
    SrcRecInfo& src = get_src_point(name_src);

    // check if this src is teleseismic source
    if (src.is_out_of_region){
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


// station correction kernel (need revise)
void InputParams::station_correction_update(CUSTOMREAL stepsize){
    if (!is_sta_correction)
        return;

    // station correction kernel is generated in the main process and sent the value to all other processors

    // step 1, gather all arrival time info to the main process
    gather_all_arrival_times_to_main();

    // do it in the main processor
    if (id_sim == 0){

        // step 2 initialize the kernel K_{\hat T_i}
        for (auto iter = rec_map.begin(); iter != rec_map.end(); iter++){
            iter->second.sta_correct_kernel = 0.0;
        }
        CUSTOMREAL max_kernel = 0.0;

        // step 3, calculate the kernel
        for (auto it_src = data_map_all.begin(); it_src != data_map_all.end(); it_src++){
            for (auto  it_rec = it_src->second.begin(); it_rec != it_src->second.end(); it_rec++){
                for (const auto& data : it_rec->second){

                    // absolute traveltime
                    if (data.is_src_rec){
                        std::cout << "teleseismic data, absolute traveltime is not supported now" << std::endl;

                    // common receiver differential traveltime
                    } else if (data.is_src_pair) {
                        std::cout << "teleseismic data, common receiver differential traveltime is not supported now" << std::endl;

                    // common source differential traveltime
                    } else if (data.is_rec_pair) {
                        std::string name_src   = data.name_src;
                        std::string name_rec1  = data.name_rec_pair[0];
                        std::string name_rec2  = data.name_rec_pair[1];

                        CUSTOMREAL syn_dif_time = get_data_rec_pair(data_map_all[name_src][name_rec1]).travel_time \
                                                - get_data_rec_pair(data_map_all[name_src][name_rec2]).travel_time;
                        CUSTOMREAL obs_dif_time = data.cs_dif_travel_time_obs;
                        rec_map[name_rec1].sta_correct_kernel += _2_CR *(syn_dif_time - obs_dif_time \
                                    + rec_map[name_rec1].sta_correct - rec_map[name_rec2].sta_correct)*data.weight;
                        rec_map[name_rec2].sta_correct_kernel -= _2_CR *(syn_dif_time - obs_dif_time \
                                    + rec_map[name_rec1].sta_correct - rec_map[name_rec2].sta_correct)*data.weight;
                        max_kernel = std::max(max_kernel,rec_map[name_rec1].sta_correct_kernel);
                        max_kernel = std::max(max_kernel,rec_map[name_rec2].sta_correct_kernel);
                    }

                }
            }
        }

        // step 4, update station correction
        for (auto iter = rec_map.begin(); iter!=rec_map.end(); iter++){
            iter->second.sta_correct += iter->second.sta_correct_kernel / (-max_kernel) * stepsize;
        }
    }

    // step 5, broadcast the station correction all all procesors
    for (auto iter = rec_map.begin(); iter!=rec_map.end(); iter++){
        broadcast_cr_single_inter_and_intra_sim(iter->second.sta_correct,0);
    }

}

void InputParams::modift_swapped_source_location() {
    for(auto iter = rec_map.begin(); iter != rec_map.end(); iter++){
        src_map_back[iter->first].lat   =   iter->second.lat;
        src_map_back[iter->first].lon   =   iter->second.lon;
        src_map_back[iter->first].dep   =   iter->second.dep;
        src_map_back[iter->first].sec   =   iter->second.sec + iter->second.tau_opt;
    }
}


template <typename T>
void InputParams::allreduce_rec_map_var(T& var){

    T tmp_var = (T)var;

    if (subdom_main){
        // allreduce sum the variable var of rec_map[name_rec] to all
        // some process has rec_map[name_rec], some process does not have it

        // step 1, gather all the variable var

        // allreduce the variable var to the main process
        // if T is CUSTOMREAL
        if (std::is_same<T, CUSTOMREAL>::value){
            CUSTOMREAL v = tmp_var; // for compiler warning
            allreduce_cr_sim_single_inplace(v);
            tmp_var = v; // for compiler warning
        // if T is int
        } else if (std::is_same<T, int>::value){
            int v = tmp_var; // for compiler warning
            allreduce_i_sim_single_inplace(v);
            tmp_var = v; // for compiler warning
        } else {
            //error
            std::cout << "error in allreduce_rec_map_var" << std::endl;
            exit(1);
        }
    }

    // assign the value to the variable var
    var = tmp_var;

//    // broadcast the variable var to all subdomains within the same simultaneous run group
//    if (std::is_same<T, CUSTOMREAL>::value){
//        CUSTOMREAL v = var; // for compiler warning
//        broadcast_cr_single(v,0);
//        var = v; // for compiler warning
//    // if T is int
//    } else if (std::is_same<T, int>::value){
//        int v = var; // for compiler warning
//        broadcast_i_single(v,0);
//        var = v; // for compiler warning
//    } else {
//        //error
//        std::cout << "error in allreduce_rec_map_var" << std::endl;
//        exit(1);
//    }

}


//
// communication for unevenly distributed receiver map
//
void InputParams::allreduce_rec_map_tau_opt(){
    if(subdom_main){
        // send total number of rec_map_all.size() to all processors
        int n_rec_all;
        std::vector<std::string> name_rec_all;
        if (id_sim == 0){
            n_rec_all = rec_map_all.size();
            for (auto iter = rec_map_all.begin(); iter != rec_map_all.end(); iter++){
                name_rec_all.push_back(iter->first);
            }
        }

        // broadcast n_rec_all to all processors
        broadcast_i_single_inter_sim(n_rec_all,0);

        for (int i_rec = 0; i_rec < n_rec_all; i_rec++){
            // broadcast name_rec_all[i_rec] to all processors
            std::string name_rec;
            if (id_sim == 0)
                name_rec = name_rec_all[i_rec];

            broadcast_str_inter_sim(name_rec,0);

            // check if the tau_opt of rec_map_all[name_rec] is needed
            bool is_stop = false;
            if (rec_map.find(name_rec) != rec_map.end()){
                is_stop = rec_map[name_rec].is_stop;
            }

            // allreduce
            allreduce_bool_single_inplace_sim(is_stop);

            // stop allreduce of tau_opt if is_stop is true (no further addition of tau_opt is needed)
            if (is_stop)
                continue;

            // allreduce the tau_opt of rec_map_all[name_rec] to all processors
            if (rec_map.find(name_rec) != rec_map.end()){
                allreduce_rec_map_var(rec_map[name_rec].tau_opt);
            } else {
                CUSTOMREAL dummy = 0;
                allreduce_rec_map_var(dummy);
            }
        }
    }
}


void InputParams::allreduce_rec_map_sum_weight(){
    if(subdom_main){
        // send total number of rec_map_all.size() to all processors
        int n_rec_all;
        std::vector<std::string> name_rec_all;
        if (id_sim == 0){
            n_rec_all = rec_map_all.size();
            for (auto iter = rec_map_all.begin(); iter != rec_map_all.end(); iter++){
                name_rec_all.push_back(iter->first);
            }
        }

        // broadcast n_rec_all to all processors
        broadcast_i_single_inter_sim(n_rec_all,0);

        for (int i_rec = 0; i_rec < n_rec_all; i_rec++){
            // broadcast name_rec_all[i_rec] to all processors
            std::string name_rec;
            if (id_sim == 0)
                name_rec = name_rec_all[i_rec];

            broadcast_str_inter_sim(name_rec,0);

            // allreduce the sum_weight of rec_map_all[name_rec] to all processors
            if (rec_map.find(name_rec) != rec_map.end()){
                allreduce_rec_map_var(rec_map[name_rec].sum_weight);
            } else {
                CUSTOMREAL dummy = 0;
                allreduce_rec_map_var(dummy);
            }
        }
    }
}


void InputParams::allreduce_rec_map_vobj_src_reloc(){
    if(subdom_main){
        // send total number of rec_map_all.size() to all processors
        int n_rec_all;
        std::vector<std::string> name_rec_all;
        if (id_sim == 0){
            n_rec_all = rec_map_all.size();
            for (auto iter = rec_map_all.begin(); iter != rec_map_all.end(); iter++){
                name_rec_all.push_back(iter->first);
            }
        }

        // broadcast n_rec_all to all processors
        broadcast_i_single_inter_sim(n_rec_all,0);

        for (int i_rec = 0; i_rec < n_rec_all; i_rec++){
            // broadcast name_rec_all[i_rec] to all processors
            std::string name_rec;
            if (id_sim == 0)
                name_rec = name_rec_all[i_rec];

            broadcast_str_inter_sim(name_rec,0);

            // allreduce the vobj_src_reloc of rec_map_all[name_rec] to all processors
            if (rec_map.find(name_rec) != rec_map.end()){
                allreduce_rec_map_var(rec_map[name_rec].vobj_src_reloc);
            } else {
                CUSTOMREAL dummy = 0;
                allreduce_rec_map_var(dummy);
            }
        }
    }
}


void InputParams::allreduce_rec_map_grad_tau(){
    if(subdom_main){
        // send total number of rec_map_all.size() to all processors
        int n_rec_all;
        std::vector<std::string> name_rec_all;
        if (id_sim == 0){
            n_rec_all = rec_map_all.size();
            for (auto iter = rec_map_all.begin(); iter != rec_map_all.end(); iter++){
                name_rec_all.push_back(iter->first);
            }
        }

        // broadcast n_rec_all to all processors
        broadcast_i_single_inter_sim(n_rec_all,0);

        for (int i_rec = 0; i_rec < n_rec_all; i_rec++){
            // broadcast name_rec_all[i_rec] to all processors
            std::string name_rec;
            if (id_sim == 0)
                name_rec = name_rec_all[i_rec];

            broadcast_str_inter_sim(name_rec,0);

            // allreduce the grad_tau of rec_map_all[name_rec] to all processors
            if (rec_map.find(name_rec) != rec_map.end()){
                allreduce_rec_map_var(rec_map[name_rec].grad_tau);
            } else {
                CUSTOMREAL dummy = 0;
                allreduce_rec_map_var(dummy);
            }
        }
    }
}


void InputParams::allreduce_rec_map_grad_chi_ijk(){
    if(subdom_main){
        // send total number of rec_map_all.size() to all processors
        int n_rec_all;
        std::vector<std::string> name_rec_all;
        if (id_sim == 0){
            n_rec_all = rec_map_all.size();
            for (auto iter = rec_map_all.begin(); iter != rec_map_all.end(); iter++){
                name_rec_all.push_back(iter->first);
            }
        }

        // broadcast n_rec_all to all processors
        broadcast_i_single_inter_sim(n_rec_all,0);

        for (int i_rec = 0; i_rec < n_rec_all; i_rec++){
            // broadcast name_rec_all[i_rec] to all processors
            std::string name_rec;
            if (id_sim == 0)
                name_rec = name_rec_all[i_rec];

            broadcast_str_inter_sim(name_rec,0);

            // allreduce the grad_chi_ijk of rec_map_all[name_rec] to all processors
            if (rec_map.find(name_rec) != rec_map.end()){
                allreduce_rec_map_var(rec_map[name_rec].grad_chi_i);
                allreduce_rec_map_var(rec_map[name_rec].grad_chi_j);
                allreduce_rec_map_var(rec_map[name_rec].grad_chi_k);
            } else {
                CUSTOMREAL dummy = 0;
                allreduce_rec_map_var(dummy);
                dummy = 0;
                allreduce_rec_map_var(dummy);
                dummy = 0;
                allreduce_rec_map_var(dummy);
            }
        }
    }
}