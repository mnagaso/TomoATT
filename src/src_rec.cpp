#include "src_rec.h"


//
// functions for processing src_rec_file
//

void parse_src_rec_file(std::string& src_rec_file, \
                        std::map<std::string, SrcRecInfo>& src_list, \
                        std::map<std::string, SrcRecInfo>& rec_list, \
                        std::vector<DataInfo>& data_info, \
                        std::vector<std::string>& src_name_list){

    // start timer
    std::string timer_name = "parse_src_rec_file";
    Timer timer(timer_name);

    std::ifstream ifs; // dummy for all processes except world_rank 0
    std::stringstream ss_whole; // for parsing the whole file
    std::stringstream ss; // for parsing each line

    // only world_rank 0 reads the file
    if (sim_rank == 0){
        ifs.open(src_rec_file);
        // abort if file does not exist
        if (!ifs.is_open()){
            std::cerr << "Error: src_rec_file " << src_rec_file << " does not exist!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        ss_whole << ifs.rdbuf();
        ifs.close();
    }

    std::string line;
    int cc = 0; // count the number of lines
    int i_src_now = 0; // count the number of srcs
    int i_rec_now = 0; // count the number of receivers
    int ndata_tmp = 0; // count the number of receivers or differential traveltime data for each source
    src_list.clear();
    rec_list.clear();
    data_info.clear();

    std::string src_name;
    CUSTOMREAL src_weight = 1.0;
    CUSTOMREAL rec_weight = 1.0;
    int src_id = -1;

    while (true) {

        bool end_of_file = false;
        bool skip_this_line = false;

        line.clear(); // clear the line before use

        // read a line
        if (sim_rank == 0){
            if (!std::getline(ss_whole, line))
                end_of_file = true;
        }

        // broadcast end_of_file
        broadcast_bool_single(end_of_file, 0);

        if (end_of_file)
            break;

        // skip comment and empty lines
        if (sim_rank == 0){
            if (line[0] == '#' || line.empty())
                skip_this_line = true;
        }

        broadcast_bool_single(skip_this_line, 0);

        if (skip_this_line)
            continue;

        // parse the line
        int ntokens = 0;
        std::string token;
        std::vector<std::string> tokens;

        if (sim_rank==0){
            // erase the trailing space
            line.erase(line.find_last_not_of(" \n\r\t")+1);

            // parse the line with arbitrary number of spaces
            ss.clear(); // clear the stringstream before use
            ss << line;

            while (std::getline(ss, token, ' ')) {
                if (token.size() > 0) // skip the first spaces and multiple spaces
                    tokens.push_back(token);
            }

            // length of tokens
            ntokens = tokens.size();
        }

        // broadcast ntokens
        broadcast_i_single(ntokens, 0);
        // broadcast tokens
        for (int i=0; i<ntokens; i++){
            if (sim_rank == 0)
                token = tokens[i];
            broadcast_str(token, 0);
            if (sim_rank != 0)
                tokens.push_back(token);
        }

        try { // check failure of parsing line by line

            // store values into structure
            if (cc == 0){
                SrcRecInfo src;

                src.id         = std::stoi(tokens[0]);
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
                src.name       = tokens[12];
                cc++;

                // check if tokens[13] exists, then read weight
                if (tokens.size() > 13)
                    src_weight = static_cast<CUSTOMREAL>(std::stod(tokens[13]));
                else
                    src_weight = 1.0; // default weight

                // new source detected by its name
                if (src_list.find(src.name) == src_list.end())
                    src_list[src.name] = src;

                src_id = src.id;
                src_name = src.name;

                ndata_tmp = src.n_data;
                src_name_list.push_back(src_name);

                // source with no receiver is allowed for solver_only model
                if (ndata_tmp==0) {
                    // go to the next source
                    cc = 0;
                    // timer
                    if (i_src_now % 100 == 0 && world_rank == 0) {
                        std::cout << "reading source " << i_src_now << " finished in " << timer.get_t() << " seconds. dt = " << timer.get_t_delta() << " seconds. \n";
                    }
                }

            } else {

                // read single receiver or differential traveltime data
                if (tokens.size() < 11) {

                    SrcRecInfo rec;

                    rec.id       = std::stoi(tokens[1]);
                    rec.name     = tokens[2];
                    rec.lat      = static_cast<CUSTOMREAL>(std::stod(tokens[3])); // in degree
                    rec.lon      = static_cast<CUSTOMREAL>(std::stod(tokens[4])); // in degree
                    rec.dep      = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[5])/1000.0); // convert elevation in meter to depth in km

                    // new receiver detected by its name
                    if(rec_list.find(rec.name) == rec_list.end())
                        rec_list[rec.name] = rec;

                    // traveltime data
                    DataInfo data;
                    if (tokens.size() > 8)
                        rec_weight = static_cast<CUSTOMREAL>(std::stod(tokens[8]));
                    else
                        rec_weight = 1.0; // default weight

                    data.data_weight = src_weight * rec_weight;
                    data.weight      = data.weight * abs_time_local_weight;
                    data.phase       = tokens[6];

                    data.is_src_rec      = true;
                    data.id_src          = src_id;
                    data.name_src        = src_name;
                    data.id_rec          = rec.id;
                    data.name_rec        = rec.name;
                    data.travel_time_obs = static_cast<CUSTOMREAL>(std::stod(tokens[7])); // store read data

                    data_info.push_back(data);
                    cc++;

                } else {

                    // read differential traveltime
                    SrcRecInfo rec;
                    rec.id   = std::stoi(tokens[1]);
                    rec.name = tokens[2];
                    rec.lat = static_cast<CUSTOMREAL>(std::stod(tokens[3])); // in degree
                    rec.lon = static_cast<CUSTOMREAL>(std::stod(tokens[4])); // in degree
                    rec.dep = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[5])/1000.0); // convert elevation in meter to depth in km

                    // new receiver detected by its name
                    if(rec_list.find(rec.name) == rec_list.end())
                        rec_list[rec.name] = rec;

                    SrcRecInfo rec2;
                    rec2.id   = std::stoi(tokens[6]);
                    rec2.name = tokens[7];
                    rec2.lat = static_cast<CUSTOMREAL>(std::stod(tokens[8])); // in degree
                    rec2.lon = static_cast<CUSTOMREAL>(std::stod(tokens[9])); // in degree
                    rec2.dep = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[10])/1000.0); // convert elevation in meter to depth in km

                    // new receiver detected by its name
                    if(rec_list.find(rec2.name) == rec_list.end())
                        rec_list[rec2.name] = rec2;

                    // common source differential traveltime
                    DataInfo data;
                    if (tokens.size() > 13)
                        rec_weight = static_cast<CUSTOMREAL>(std::stod(tokens[13]));
                    else
                        rec_weight = 1.0; // default weight

                    data.data_weight = src_weight * rec_weight;
                    data.weight      = data.weight * cs_dif_time_local_weight;
                    data.phase       = tokens[11];

                    data.is_rec_pair      = true;
                    data.id_src_single   = src_id;
                    data.name_src_single = src_name;
                    data.id_rec_pair   = {rec.id, rec2.id};
                    data.name_rec_pair = {rec.name, rec2.name};
                    data.cs_dif_travel_time_obs = static_cast<CUSTOMREAL>(std::stod(tokens[12])); // store read data

                    data_info.push_back(data);
                    cc++;

                }

                if (cc > ndata_tmp) {
                    // go to the next source
                    cc = 0;
                    i_src_now++;
                    i_rec_now = 0;

                    // timer
                    if (i_src_now % 1000 == 0 && world_rank == 0) {
                        std::cout << "reading source " << i_src_now << " finished in " << timer.get_t() << " seconds. dt = " << timer.get_t_delta() << " seconds. \n";
                    }
                } else {
                    i_rec_now++;
                }
            }

        } catch (std::invalid_argument& e) {
                std::cout << "Error: invalid argument in src_rec_file. Abort." << std::endl;
                std::cout << "problematic line: \n\n" << line << std::endl;

                MPI_Abort(MPI_COMM_WORLD, 1);
        }

    /*
        // print for DEBUG
        for (auto& t : tokens) {
            std::cout << t << "---";
        }
        std::cout << std::endl;
    */

    } // end of while loop

    // abort if number of src_points are less than n_sims
    int n_src_points = src_list.size();
    if (n_src_points < n_sims){
        std::cout << "Error: number of sources in src_rec_file is less than n_sims. Abort." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // indicate elapsed time
    std::cout << "Total elapsed time for reading src_rec_file: " << timer.get_t() << " seconds.\n";

    // check new version of src rec data
    if (if_verbose){
        for(auto iter = src_list.begin(); iter != src_list.end(); iter++){
            std::cout   << "source id: "     << iter->second.id
                        << ", source name: " << iter->second.name
                        << std::endl;
        }

        for(auto iter = rec_list.begin(); iter != rec_list.end(); iter++){
            std::cout   << "receiver id: "     << iter->second.id
                        << ", receiver name: " << iter->second.name
                        << std::endl;
        }

        for(int i = 0; i < (int)data_info.size(); i++){
            if (data_info[i].is_src_rec){
                std::cout   << "absolute traveltime: " << data_info[i].travel_time_obs
                            << ", source name: "       << data_info[i].name_src
                            << ", receiver name: "     << data_info[i].name_rec
                            << std::endl;
            }
            if (data_info[i].is_rec_pair){
                std::cout   << "common source differential traveltime: " << data_info[i].cs_dif_travel_time_obs
                            << ", source name: "                         << data_info[i].name_src_single
                            << ", receiver pair name: "                  << data_info[i].name_rec_pair[0]
                            << ", "                                      << data_info[i].name_rec_pair[1]
                            << std::endl;
            }
        }
        std::cout << data_info.size() << std::endl;
    }
}


void parse_sta_correction_file(std::string& sta_correction_file, \
                               std::map<std::string, SrcRecInfo>& rec_list){

    // read station correction file
    std::ifstream ifs;
    std::stringstream ss, ss_whole;

    if (sim_rank == 0){
        ifs.open(sta_correction_file);
        if (!ifs.is_open()){
            std::cout << "Error: cannot open sta_correction_file. Abort." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        ss_whole << ifs.rdbuf();
        ifs.close();
    }

    std::string line;

    while(true){
        bool end_of_file = false;
        bool skip_this_line = false;

        line.clear();

        // read a line
        if (sim_rank == 0){
            if (!std::getline(ss_whole,line))
                end_of_file = true;
        }

        broadcast_bool_single(end_of_file, 0);

        if (end_of_file)
            break;

        // skip this line if it is a comment line
        if (sim_rank == 0){
            if (line[0] == '#' || line.empty())
                skip_this_line = true;
        }

        broadcast_bool_single(skip_this_line, 0);

        if (skip_this_line)
            continue;

        // parse the line
        int ntokens = 0;
        std::string token;
        std::vector<std::string> tokens;

        if (sim_rank == 0){
            // erase the last space
            line.erase(line.find_last_not_of(" \n\r\t")+1);

            // parse the line with arbitrary number of spaces
            ss.clear();
            ss << line;

            while (std::getline(ss,token,' ')) {
                if (token.size() > 0)
                    tokens.push_back(token);
            }
            // number of tokens
            ntokens = tokens.size();
        }

        // broadcast ntokens
        broadcast_i_single(ntokens, 0);
        // broadcast tokens
        for (int i=0; i<ntokens; i++){
            if (sim_rank == 0)
                token = tokens[i];
            broadcast_str(token, 0);
            if (sim_rank != 0)
                tokens.push_back(token);
        }

        try { // check failure of parsion line by line

            // store station corrections into rec_list
            std::string tmp_sta_name = tokens[0];

            if (rec_list.find(tmp_sta_name) == rec_list.end()){
                // new station
                SrcRecInfo tmp_rec;
                tmp_rec.name = tmp_sta_name;
                tmp_rec.lat  = static_cast<CUSTOMREAL>(std::stod(tokens[1])); // in degree
                tmp_rec.lon  = static_cast<CUSTOMREAL>(std::stod(tokens[2])); // in degree
                tmp_rec.dep  = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[3])/1000.0); // convert elevation in meter to depth in km
                tmp_rec.sta_correct = static_cast<CUSTOMREAL>(std::stod(tokens[4]));
                tmp_rec.sta_correct_kernel = 0.0;
                rec_list[tmp_sta_name] = tmp_rec;
            } else {
                // pre exist station
                rec_list[tmp_sta_name].sta_correct = static_cast<CUSTOMREAL>(std::stod(tokens[4]));
                rec_list[tmp_sta_name].sta_correct_kernel = 0.0;
            }
        } catch (std::invalid_argument& e) {
            std::cout << "Error: invalid argument in sta_correction_file. Abort." << std::endl;
            std::cout << "problematic line: \n\n" << line << std::endl;

            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
}


void do_swap_src_rec(std::map<std::string, SrcRecInfo> &src_list, \
                     std::map<std::string, SrcRecInfo> &rec_list, \
                     std::vector<DataInfo>             &data_info) {

    // swap src/rec points
    // at this moment, all the sources are divided into src_points (regional) and tele_src_points (teleseismic)

    // Start timer
    std::string timer_name = "swap_src_rec";
    Timer timer(timer_name);

    std::map<std::string, SrcRecInfo> tmp_src_rec_list = src_list;
    src_list = rec_list;
    rec_list = tmp_src_rec_list;

    for(int i = 0; i < (int)data_info.size(); i++){
        DataInfo tmp_data = data_info[i];
        if (tmp_data.is_src_rec){ // absolute traveltime  ->  absolute traveltime
            tmp_data.id_src   = data_info[i].id_rec;
            tmp_data.name_src = data_info[i].name_rec;
            tmp_data.id_rec   = data_info[i].id_src;
            tmp_data.name_rec = data_info[i].name_src;
        } else if (tmp_data.is_rec_pair) { // common source differential traveltime  ->  common receiver differential traveltime
            tmp_data.is_rec_pair            = false;
            tmp_data.is_src_pair            = true;
            tmp_data.id_src_pair            = data_info[i].id_rec_pair;
            tmp_data.name_src_pair          = data_info[i].name_rec_pair;
            tmp_data.id_rec_single          = data_info[i].id_src_single;
            tmp_data.name_rec_single        = data_info[i].name_src_single;
            tmp_data.cr_dif_travel_time_obs = data_info[i].cs_dif_travel_time_obs;
        } else if (tmp_data.is_src_pair) { // common receiver differential traveltime  ->  common source differential traveltime
            tmp_data.is_src_pair            = false;
            tmp_data.is_rec_pair            = true;
            tmp_data.id_rec_pair            = data_info[i].id_src_pair;
            tmp_data.name_rec_pair          = data_info[i].name_src_pair;
            tmp_data.id_src_single          = data_info[i].id_rec_single;
            tmp_data.name_src_single        = data_info[i].name_rec_single;
            tmp_data.cs_dif_travel_time_obs = data_info[i].cr_dif_travel_time_obs;
        }
        data_info[i] = tmp_data;
    }

    // swap total_data_weight
    CUSTOMREAL tmp = total_cr_dif_local_data_weight;
    total_cr_dif_local_data_weight = total_cs_dif_local_data_weight;
    total_cs_dif_local_data_weight = tmp;

    // check new version of src rec data
    if (if_verbose){     // check by Chen Jing
        for(auto iter = src_list.begin(); iter != src_list.end(); iter++){
            std::cout   << "source id: " << iter->second.id
                        << ", source name: " << iter->second.name
                        << std::endl;
        }

        for(auto iter = rec_list.begin(); iter != rec_list.end(); iter++){
            std::cout   << "receiver id: " << iter->second.id
                        << ", receiver name: " << iter->second.name
                        << std::endl;
        }

        for(int i = 0; i < (int)data_info.size(); i++){
            if (data_info[i].is_src_rec){
                std::cout   << "absolute traveltime: " << data_info[i].travel_time_obs
                            << ", source name: " << data_info[i].name_src
                            << ", receiver name: " << data_info[i].name_rec
                            << std::endl;
            } else if (data_info[i].is_rec_pair){
                std::cout   << "common source differential traveltime: " << data_info[i].cs_dif_travel_time_obs
                            << ", source name: " << data_info[i].name_src_single
                            << ", receiver pair name: " << data_info[i].name_rec_pair[0]
                            << ", " << data_info[i].name_rec_pair[1]
                            << std::endl;
            } else if (data_info[i].is_src_pair){
                std::cout   << "common receiver differential traveltime: " << data_info[i].cr_dif_travel_time_obs
                            << ", source pair name: " << data_info[i].name_src_pair[0]
                            << ", " << data_info[i].name_src_pair[1]
                            << ", receiver name: " << data_info[i].name_rec_single
                            << std::endl;
            }
        }
        std::cout << data_info.size() << std::endl;
    }
    // indicate elapsed time
    std::cout << "Total elapsed time for swapping src rec: " << timer.get_t() << " seconds.\n";
}



