#include "src_rec.h"


//
// functions for processing src_rec_file
//

void parse_src_rec_file(std::string& src_rec_file, \
                        std::map<std::string, SrcRecInfo>& src_map, \
                        std::map<std::string, SrcRecInfo>& rec_map, \
                        std::map<std::string, std::map<std::string, std::vector<DataInfo>>>& data_map, \
                        std::vector<std::string>& src_name_list, \
                        std::vector<std::vector<std::vector<std::string>>>& srcrec_name_list){

    // start timer
    std::string timer_name = "parse_src_rec_file";
    Timer timer(timer_name);

    std::ifstream ifs;          // dummy for all processes except world_rank 0
    std::stringstream ss_whole; // for parsing the whole file
    std::stringstream ss;       // for parsing each line

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
    int cc = 0;        // count the number of lines
    int i_src_now = 0; // count the number of srcs
    int ndata_tmp = 0; // count the number of receivers or differential traveltime data for each source
    src_map.clear();
    rec_map.clear();
    data_map.clear();

    std::string src_name;
    CUSTOMREAL src_weight = 1.0;
    CUSTOMREAL rec_weight = 1.0;
    int src_id = -1;

    // temporary receiver name list for each source
    std::vector<std::vector<std::string>> rec_name_list;

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
            SrcRecInfo src;
            if (cc == 0){
                
                src.id     = std::stoi(tokens[0]);
                src.year   = std::stoi(tokens[1]);
                src.month  = std::stoi(tokens[2]);
                src.day    = std::stoi(tokens[3]);
                src.hour   = std::stoi(tokens[4]);
                src.min    = std::stoi(tokens[5]);
                src.sec    = static_cast<CUSTOMREAL>(std::stod(tokens[6]));
                src.lat    = static_cast<CUSTOMREAL>(std::stod(tokens[7])); // in degree
                src.lon    = static_cast<CUSTOMREAL>(std::stod(tokens[8])); // in degree
                src.dep    = static_cast<CUSTOMREAL>(std::stod(tokens[9])); // source in km
                src.mag    = static_cast<CUSTOMREAL>(std::stod(tokens[10]));
                src.n_data = std::stoi(tokens[11]);
                src.name   = tokens[12];
                cc++;

                // check if tokens[13] exists, then read weight
                if (tokens.size() > 13)
                    src_weight = static_cast<CUSTOMREAL>(std::stod(tokens[13]));
                else
                    src_weight = 1.0; // default weight

                // new source detected by its name
                // TODO: add error check for duplicated source name (but different event info)
                if (src_map.find(src.name) == src_map.end())
                    src_map[src.name] = src;

                src_id   = src.id;
                src_name = src.name;

                ndata_tmp = src.n_data;
                src_name_list.push_back(src_name); // store order of sources in the file

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
                    // store receiver name of onle receiver line in src rec file
                    std::vector<std::string> rec_name_list_one_line;

                    SrcRecInfo rec;

                    rec.id   = std::stoi(tokens[1]);
                    rec.name = tokens[2];
                    rec.lat  = static_cast<CUSTOMREAL>(std::stod(tokens[3])); // in degree
                    rec.lon  = static_cast<CUSTOMREAL>(std::stod(tokens[4])); // in degree
                    rec.dep  = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[5])/1000.0); // convert elevation in meter to depth in km

                    // new receiver detected by its name
                    if(rec_map.find(rec.name) == rec_map.end())
                        rec_map[rec.name] = rec;

                    // store temporary receiver name list for each source
                    rec_name_list_one_line.push_back(rec.name);

                    // traveltime data
                    DataInfo data;
                    if (tokens.size() > 8)
                        rec_weight = static_cast<CUSTOMREAL>(std::stod(tokens[8]));
                    else
                        rec_weight = 1.0; // default weight

                    data.data_weight = src_weight * rec_weight;
                    data.weight      = data.data_weight * abs_time_local_weight;
                    data.phase       = tokens[6];

                    data.is_src_rec      = true;
                    data.id_src          = src_id;
                    data.name_src        = src_name;
                    data.id_rec          = rec.id;
                    data.name_rec        = rec.name;
                    data.travel_time_obs = static_cast<CUSTOMREAL>(std::stod(tokens[7])); // store read data

                    data_map[data.name_src][data.name_rec].push_back(data);

                    // store receiver name of onle receiver line in src rec file
                    rec_name_list.push_back(rec_name_list_one_line);

                    cc++;

                } else {
                    // read common source differential traveltime (cs_dif) or common receiver differential traveltime (cr_dif) 

                    std::vector<std::string> rec_name_list_one_line;

                    // read differential traveltime
                    SrcRecInfo rec;
                    rec.id   = std::stoi(tokens[1]);
                    rec.name = tokens[2];
                    rec.lat  = static_cast<CUSTOMREAL>(std::stod(tokens[3])); // in degree
                    rec.lon  = static_cast<CUSTOMREAL>(std::stod(tokens[4])); // in degree
                    rec.dep  = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[5])/1000.0); // convert elevation in meter to depth in km

                    // new receiver detected by its name
                    if(rec_map.find(rec.name) == rec_map.end())
                        rec_map[rec.name] = rec;

                    // store temporary receiver name list for each source
                    rec_name_list_one_line.push_back(rec.name);

                    
                    // differential traveltime data
                    DataInfo data;
                    if (tokens.size() > 13)
                        rec_weight = static_cast<CUSTOMREAL>(std::stod(tokens[13]));
                    else
                        rec_weight = 1.0; // default weight

                    data.data_weight = src_weight * rec_weight;
                    data.phase       = tokens[11];
                  
                    //data.id_src_single          = src_id;
                    //data.name_src_single        = src_name;
                    // use common variables with src-rec data
                    data.id_src          = src_id;
                    data.name_src        = src_name;

                    // store the id and name of the first receiver (just used for key of the data_map)
                    data.id_rec          = rec.id;
                    data.name_rec        = rec.name;
                    
                    // determine this data is cr_dif or cs_dif
                    bool is_cr_dif = tokens[11].find("cr")!=std::string::npos;

                    if (is_cr_dif) {
                        // cr_dif data
                        SrcRecInfo src2;
                        src2.id   = std::stoi(tokens[6]);
                        src2.name = tokens[7];
                        src2.lat  = static_cast<CUSTOMREAL>(std::stod(tokens[8])); // in degree
                        src2.lon  = static_cast<CUSTOMREAL>(std::stod(tokens[9])); // in degree
                        src2.dep  = static_cast<CUSTOMREAL>(std::stod(tokens[10])); // convert elevation in meter to depth in km

                        // new source detected by its name
                        if (src_map.find(src2.name) == src_map.end())
                            src_map[src2.name] = src2;

                        // store temporary receiver(source) name list for each source
                        rec_name_list_one_line.push_back(src2.name);

                        // common receiver differential traveltime data
                        data.is_src_pair            = true;
                        data.id_src_pair            = {src.id, src2.id};
                        data.name_src_pair          = {src.name, src2.name};
                        data.cr_dif_travel_time_obs = static_cast<CUSTOMREAL>(std::stod(tokens[12])); // store read data
                        data_map[data.name_src_pair[0]][data.name_rec].push_back(data); // USE ONE-DATAMAP-FOR-ONE-SRCREC-LINE

                        data.weight = data.data_weight * cr_dif_time_local_weight;
                    } else {
                        // cs_dif data     
                        SrcRecInfo rec2;                   
                        rec2.id   = std::stoi(tokens[6]);
                        rec2.name = tokens[7];
                        rec2.lat  = static_cast<CUSTOMREAL>(std::stod(tokens[8])); // in degree
                        rec2.lon  = static_cast<CUSTOMREAL>(std::stod(tokens[9])); // in degree
                        rec2.dep  = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[10])/1000.0); // convert elevation in meter to depth in km

                        // new receiver detected by its name
                        if(rec_map.find(rec2.name) == rec_map.end())
                            rec_map[rec2.name] = rec2;

                        // store temporary receiver name list for each source
                        rec_name_list_one_line.push_back(rec2.name);

                        // common source differential traveltime data
                        data.is_rec_pair            = true;
                        data.id_rec_pair            = {rec.id, rec2.id};
                        data.name_rec_pair          = {rec.name, rec2.name};
                        data.cs_dif_travel_time_obs = static_cast<CUSTOMREAL>(std::stod(tokens[12])); // store read data
                        data_map[data.name_src][data.name_rec_pair[0]].push_back(data); // USE ONE-DATAMAP-FOR-ONE-SRCREC-LINE
                        //data_map[data.name_src][data.name_rec_pair[1]].push_back(data); // TODO: check if name_rec_pair[1] should be stored as well

                        data.weight = data.data_weight * cs_dif_time_local_weight;
                    }
                    
                    // store receiver name of onle receiver line in src rec file
                    rec_name_list.push_back(rec_name_list_one_line);

                    cc++;
                }

                if (cc > ndata_tmp) {
                    // go to the next source
                    cc = 0;
                    i_src_now++;

                    // store the receiver name list for the source
                    srcrec_name_list.push_back(rec_name_list);
                    // clear the temporary receiver name list
                    rec_name_list.clear();

                    // timer
                    if (i_src_now % 1000 == 0 && world_rank == 0) {
                        std::cout << "reading source " << i_src_now << " finished in " << timer.get_t() << " seconds. dt = " << timer.get_t_delta() << " seconds. \n";
                    }
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

    // indicate the number of sources and receivers, data_info
    if (world_rank == 0){
        std::cout << "\nReading src_rec_file finished." << std::endl;
        std::cout << "number of sources: "   << src_map.size()  << std::endl;
        std::cout << "number of receivers: " << rec_map.size()  << std::endl;
        std::cout << "number of data: "      << data_map.size() << "\n" << std::endl;
    }

    // indicate elapsed time
    std::cout << "Total elapsed time for reading src_rec_file: " << timer.get_t() << " seconds.\n";

    // check new version of src rec data
    if (if_verbose){
        for(auto iter = src_map.begin(); iter != src_map.end(); iter++){
            std::cout   << "source id: "     << iter->second.id
                        << ", source name: " << iter->second.name
                        << std::endl;
        }

        for(auto iter = rec_map.begin(); iter != rec_map.end(); iter++){
            std::cout   << "receiver id: "     << iter->second.id
                        << ", receiver name: " << iter->second.name
                        << std::endl;
        }

        for (auto iter = data_map.begin(); iter != data_map.end(); iter++){
            for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++){
                for (const auto& data : iter2->second) {
                    if (data.is_src_rec) {
                        std::cout   << "source name: "     << data.name_src
                                    << ", receiver name: " << data.name_rec
                                    << ", traveltime: "    << data.travel_time_obs
                                    << std::endl;
                    } else if (data.is_rec_pair) {
                        std::cout   << "source name: "          << data.name_src
                                    << ", receiver pair name: " << data.name_rec_pair[0]
                                    << ", "                     << data.name_rec_pair[1]
                                    << ", traveltime: "         << data.cs_dif_travel_time_obs
                                    << std::endl;
                    }
                }
           }
        }

    }
}


void parse_sta_correction_file(std::string& sta_correction_file, \
                               std::map<std::string, SrcRecInfo>& rec_map){

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

            // store station corrections into rec_map
            std::string tmp_sta_name = tokens[0];

            if (rec_map.find(tmp_sta_name) == rec_map.end()){
                // new station
                SrcRecInfo tmp_rec;
                tmp_rec.name = tmp_sta_name;
                tmp_rec.lat  = static_cast<CUSTOMREAL>(std::stod(tokens[1])); // in degree
                tmp_rec.lon  = static_cast<CUSTOMREAL>(std::stod(tokens[2])); // in degree
                tmp_rec.dep  = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[3])/1000.0); // convert elevation in meter to depth in km
                tmp_rec.sta_correct = static_cast<CUSTOMREAL>(std::stod(tokens[4]));
                tmp_rec.sta_correct_kernel = 0.0;
                rec_map[tmp_sta_name] = tmp_rec;
            } else {
                // pre exist station
                rec_map[tmp_sta_name].sta_correct = static_cast<CUSTOMREAL>(std::stod(tokens[4]));
                rec_map[tmp_sta_name].sta_correct_kernel = 0.0;
            }
        } catch (std::invalid_argument& e) {
            std::cout << "Error: invalid argument in sta_correction_file. Abort." << std::endl;
            std::cout << "problematic line: \n\n" << line << std::endl;

            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
}


void separate_region_and_tele_src_rec_data(std::map<std::string, SrcRecInfo>                                  &src_map_back,
                                           std::map<std::string, SrcRecInfo>                                  &rec_map_back,
                                           std::map<std::string, std::map<std::string,std::vector<DataInfo>>> &data_map_back,
                                           std::map<std::string, SrcRecInfo>                                  &src_map,
                                           std::map<std::string, SrcRecInfo>                                  &rec_map,
                                           std::map<std::string, std::map<std::string,std::vector<DataInfo>>> &data_map,
                                           std::map<std::string, SrcRecInfo>                                  &src_map_tele,
                                           std::map<std::string, SrcRecInfo>                                  &rec_map_tele,
                                           std::map<std::string, std::map<std::string,std::vector<DataInfo>>> &data_map_tele,
                                           std::map<std::string, int> &data_type,
                                           int                        &N_abs_local_data,
                                           int                        &N_cr_dif_local_data,
                                           int                        &N_cs_dif_local_data,
                                           int                        &N_teleseismic_data,
                                           int                        &N_data,
                                           const CUSTOMREAL min_lat, const CUSTOMREAL max_lat,
                                           const CUSTOMREAL min_lon, const CUSTOMREAL max_lon,
                                           const CUSTOMREAL min_dep, const CUSTOMREAL max_dep){
    // check if the source is inside the simulation boundary

    // initialize vectors for teleseismic events
    rec_map_tele.clear();
    src_map_tele.clear();
    data_map_tele.clear();

    // clear original src, rec, data
    src_map.clear();
    rec_map.clear();
    data_map.clear();

    // divide source list
    //
    // src_map_back ____> src_map
    //              \___> src_map_tele (teleseismic events)
    //
    for(auto iter = src_map_back.begin(); iter != src_map_back.end(); iter++){
        SrcRecInfo src = iter->second;
        if (src.lat < min_lat || src.lat > max_lat \
         || src.lon < min_lon || src.lon > max_lon \
         || src.dep < min_dep || src.dep > max_dep){

            // out of region (teleseismic events)
            src.is_out_of_region       = true;
            src_map_tele[iter->first] = src;
        } else {
            // within region (local events)
            src.is_out_of_region  = false;
            src_map[iter->first] = src;
        }
    }

    // divide receiver and data list
    //
    // rec_map_back  ____> rec_map
    //               \___> rec_map_tele
    // data_map_back____> data_map
    //              \___> data_map_tele
    //
    // iterate over event
    for (auto it_src = data_map_back.begin(); it_src != data_map_back.end(); it_src++){
        std::string name_src = it_src->first;
        // iterate over stations
        for (auto it_rec = it_src->second.begin(); it_rec != it_src->second.end(); it_rec++){
            // loop over all data belongs to each source-receiver combination
            for (auto& data: it_rec->second){

                // absolute traveltime
                if(data.is_src_rec){
                    std::string name_src = data.name_src;
                    std::string name_rec = data.name_rec;

                    // if name_src is out of region
                    if(src_map_tele.find(name_src) != src_map_tele.end()){
                        total_teleseismic_data_weight += data.data_weight;
                        data.weight                    = data.data_weight * teleseismic_weight;
                        data_map_tele[name_src][name_rec].push_back(data);
                        rec_map_tele[name_rec]            = rec_map_back[name_rec];
                        data_type["tele"]                 = 1;

                    // if name_src is in the region
                    } else {
                        total_abs_local_data_weight += data.data_weight;
                        data_map[name_src][name_rec].push_back(data);
                        rec_map[name_rec]            = rec_map_back[name_rec];
                        data_type["abs"]             = 1;
                    }

                // common receiver differential traveltime
                } else if (data.is_src_pair){
                    std::string name_src1 = data.name_src_pair[0];
                    std::string name_src2 = data.name_src_pair[1];
                    std::string name_rec  = data.name_rec;

                    // if both sources are out of region
                    if(src_map_tele.find(name_src1) != src_map_tele.end() \
                    && src_map_tele.find(name_src2) != src_map_tele.end()){
                        total_teleseismic_data_weight     += data.data_weight;
                        data.weight                        = data.data_weight * teleseismic_weight;
                        data_map_tele[name_src1][name_rec].push_back(data);
                        rec_map_tele[name_rec]             = rec_map_back[name_rec];
                        data_type["tele"]                  = 1;

                    // if both sources is in the region
                    } else if (src_map.find(name_src1) != src_map.end() \
                            && src_map.find(name_src2) != src_map.end() ) {
                        total_cr_dif_local_data_weight += data.data_weight;
                        data_map[name_src1][name_rec].push_back(data);
                        rec_map[name_rec]               = rec_map_back[name_rec];
                        data_type["cr_dif"]             = 1;

                    } else {
                        std::cout << "ERROR data: common receiver differential time, but one teleseismic source, one local source";
                        exit(1);
                    }

                // common source differential traveltime
                } else if (data.is_rec_pair){
                    std::string name_src  = data.name_src;
                    std::string name_rec1 = data.name_rec_pair[0];
                    std::string name_rec2 = data.name_rec_pair[1];

                    // if name_src is out of region
                    if(src_map_tele.find(name_src) != src_map_tele.end() ){
                        total_teleseismic_data_weight     += data.data_weight;
                        data.weight                        = data.data_weight * teleseismic_weight;
                        data_map_tele[name_src][name_rec1].push_back(data);
                        rec_map_tele[name_rec1]            = rec_map_back[name_rec1];
                        rec_map_tele[name_rec2]            = rec_map_back[name_rec2];
                        data_type["tele"]                  = 1;

                    // if name_src is in the region
                    } else {
                        total_cs_dif_local_data_weight += data.data_weight;
                        data_map[name_src][name_rec1].push_back(data);
                        rec_map[name_rec1]              = rec_map_back[name_rec1];
                        rec_map[name_rec2]              = rec_map_back[name_rec2];
                        data_type["cs_dif"]             = 1;
                    }
                }

            } // end loop std::vector<datainfo>
        } // end it_rec
    } // end it_src

    //
    // balance the data weight
    //
    if (is_balance_data_weight){
        for(auto it_src = data_map.begin(); it_src != data_map.end(); it_src++){
            for(auto it_rec = it_src->second.begin(); it_rec != it_src->second.end(); it_rec++){
                for(auto& data: it_rec->second){
                    // absolute traveltime
                   if(data.is_src_rec){
                       data.weight = data.weight / total_abs_local_data_weight;

                   // common receiver differential traveltime
                   } else if (data.is_src_pair){
                       data.weight = data.weight / total_cr_dif_local_data_weight;

                   // common source differential traveltime
                   } else if (data.is_rec_pair){
                       data.weight = data.weight / total_cs_dif_local_data_weight;
                   }
                }
            }
        }

        // teleseismic data
        for(auto it_src = data_map_tele.begin(); it_src != data_map_tele.end(); it_src++){
            for(auto it_rec = it_src->second.begin(); it_rec != it_src->second.end(); it_rec++){
                for (auto& data: it_rec->second){
                    data.weight = data.weight / total_teleseismic_data_weight;
                }
            }
        }
    }

    //
    // count the number of data
    //
    // local data
    for(auto it_src = data_map.begin(); it_src != data_map.end(); it_src++){
        for(auto it_rec = it_src->second.begin(); it_rec != it_src->second.end(); it_rec++){
            for (auto& data: it_rec->second){
                // absolute traveltime
                if(data.is_src_rec){
                    N_abs_local_data += 1;

                // common receiver differential traveltime
                } else if (data.is_src_pair){
                    N_cr_dif_local_data += 1;

                // common source differential traveltime
                } else if (data.is_rec_pair){
                    N_cs_dif_local_data += 1;
                }
            }
        }
    }

    // teleseismic data
    for(auto it_src = data_map_tele.begin(); it_src != data_map_tele.end(); it_src++){
        for(auto it_rec = it_src->second.begin(); it_rec != it_src->second.end(); it_rec++){
            N_teleseismic_data += it_rec->second.size(); // add the number of data
        }
    }

    // N_data is the total number of data
    N_data = N_abs_local_data + N_cr_dif_local_data + N_cs_dif_local_data + N_teleseismic_data;

    // std::cout << "N_abs_local_data: " << N_abs_local_data << ", N_cr_dif_local_data" << N_cr_dif_local_data
    //           << ", N_cs_dif_local_data: " << N_cs_dif_local_data << ", N_teleseismic_data: " << N_teleseismic_data
    //           << std::endl
    //           << std::endl;

}


void do_swap_src_rec(std::map<std::string, SrcRecInfo> &src_map, \
                     std::map<std::string, SrcRecInfo> &rec_map, \
                     std::map<std::string, std::map<std::string, std::vector<DataInfo>>> &data_map,
                     std::vector<std::string> &src_name_list) {

    // swap src/rec points
    // at this moment, all the sources are divided into src_points (regional) and tele_src_points (teleseismic)

    // Start timer
    std::string timer_name = "swap_src_rec";
    Timer timer(timer_name);

    std::map<std::string, SrcRecInfo> tmp_src_rec_map = src_map;
    src_map = rec_map;
    rec_map = tmp_src_rec_map;

    std::map<std::string, std::map<std::string, std::vector<DataInfo>>> tmp_data_map;// = data_map;

    // for each element of src_map, count the number of rec_map with the same value of

    for (auto it_src = data_map.begin(); it_src != data_map.end(); it_src++){ // loop over src
        for (auto it_rec = it_src->second.begin(); it_rec != it_src->second.end(); it_rec++){ // loop over rec
            for (const auto& data: it_rec->second){ // loop over datainfo

                DataInfo tmp_data = data;

                // absolute traveltime  ->  absolute traveltime
                if (tmp_data.is_src_rec){
                    tmp_data.id_src   = data.id_rec;
                    tmp_data.name_src = data.name_rec;
                    tmp_data.id_rec   = data.id_src;
                    tmp_data.name_rec = data.name_src;

                    tmp_data_map[tmp_data.name_src][tmp_data.name_rec].push_back(tmp_data);

                // common source differential traveltime  ->  common receiver differential traveltime
                } else if (tmp_data.is_rec_pair) {
                    tmp_data.is_rec_pair            = false;
                    tmp_data.is_src_pair            = true;
                    tmp_data.id_src_pair            = data.id_rec_pair;
                    tmp_data.name_src_pair          = data.name_rec_pair;
                    tmp_data.id_rec                 = data.id_src;
                    tmp_data.name_rec               = data.name_src;
                    tmp_data.cr_dif_travel_time_obs = data.cs_dif_travel_time_obs;

                    tmp_data.name_rec = tmp_data.name_rec;
                    tmp_data.id_rec   = tmp_data.id_rec;

                    tmp_data.id_rec_pair            = {-1, -1};
                    tmp_data.name_rec_pair          = {"unknown", "unknown"};
                    //tmp_data.id_src_single          = -1;

                    tmp_data.name_src = tmp_data.name_src_pair[0];
                    tmp_data_map[tmp_data.name_src_pair[0]][tmp_data.name_rec].push_back(tmp_data);
                    tmp_data.name_src = tmp_data.name_src_pair[1];
                    tmp_data_map[tmp_data.name_src_pair[1]][tmp_data.name_rec].push_back(tmp_data);

                // common receiver differential traveltime  ->  common source differential traveltime
                } else if (tmp_data.is_src_pair) {
                    tmp_data.is_src_pair            = false;
                    tmp_data.is_rec_pair            = true;
                    tmp_data.id_rec_pair            = data.id_src_pair;
                    tmp_data.name_rec_pair          = data.name_src_pair;
                    tmp_data.id_src                 = data.id_rec;
                    tmp_data.name_src               = data.name_rec;
                    tmp_data.cs_dif_travel_time_obs = data.cr_dif_travel_time_obs;

                    tmp_data.id_src_pair            = {-1, -1};
                    tmp_data.name_src_pair          = {"unknown", "unknown"};
                    tmp_data.id_rec                 = -1;

                    tmp_data.name_rec = tmp_data.name_rec_pair[0];
                    tmp_data_map[tmp_data.name_src][tmp_data.name_rec_pair[0]].push_back(tmp_data);
                    tmp_data.name_rec = tmp_data.name_rec_pair[1];
                    tmp_data_map[tmp_data.name_src][tmp_data.name_rec_pair[1]].push_back(tmp_data);
                }

           }
        }
    }

    // replace data_map with swapped data map
    // TODO tmp_data_map has strange elements here
    data_map = tmp_data_map;

    // set n_data (number of receivers for each source)
    for (auto it_src = src_map.begin(); it_src != src_map.end(); it_src++){
        it_src->second.n_data = data_map[it_src->second.name].size();
    }

    // swap total_data_weight
    CUSTOMREAL tmp = total_cr_dif_local_data_weight;
    total_cr_dif_local_data_weight = total_cs_dif_local_data_weight;
    total_cs_dif_local_data_weight = tmp;

    // create new src_id2name_all from swapped src_map
    src_name_list.clear();
    for (auto iter = src_map.begin(); iter != src_map.end(); iter++){
        src_name_list.push_back(iter->first);
    }

    // check new version of src rec data
    if (if_verbose){     // check by Chen Jing
        std::cout << "do swap sources and receivers" << std::endl;

        for(auto iter = src_map.begin(); iter != src_map.end(); iter++){
            std::cout   << "source id: " << iter->second.id
                        << ", source name: " << iter->second.name
                        << std::endl;
        }

        for(auto iter = rec_map.begin(); iter != rec_map.end(); iter++){
            std::cout   << "receiver id: " << iter->second.id
                        << ", receiver name: " << iter->second.name
                        << std::endl;
        }

        for(auto it_src = data_map.begin(); it_src != data_map.end(); it_src++){
            for(auto it_rec = it_src->second.begin(); it_rec != it_src->second.end(); it_rec++){
                for(const auto& data: it_rec->second){

                    if (data.is_src_rec){
                        std::cout   << "key: it_src, it_rec: " << it_src->first << ", " << it_rec->first
                                    << "absolute traveltime: " << data.travel_time_obs
                                    << ", source name: "       << data.name_src
                                    << ", receiver name: "     << data.name_rec
                                    << std::endl;
                    } else if (data.is_rec_pair){
                        std::cout   << "key: it_src, it_rec: "  << it_src->first << ", " << it_rec->first
                                    << "common source differential traveltime: " << data.cs_dif_travel_time_obs
                                    << ", source name: "                         << data.name_src
                                    << ", receiver pair name: "                  << data.name_rec_pair[0]
                                    << ", "                                      << data.name_rec_pair[1]
                                    << std::endl;
                    } else if (data.is_src_pair){
                        std::cout   << "key: it_src, it_rec: "    << it_src->first << ", " << it_rec->first
                                    << "common receiver differential traveltime: " << data.cr_dif_travel_time_obs
                                    << ", source pair name: "                      << data.name_src_pair[0]
                                    << ", "                                        << data.name_src_pair[1]
                                    << ", receiver name: "                         << data.name_rec
                                    << std::endl;
                    }
                }
            }
        }
        std::cout << data_map.size() << std::endl;
    }

    // indicate elapsed time
    std::cout << "Total elapsed time for swapping src rec: " << timer.get_t() << " seconds.\n";
}




// merge the teleseismic data lsit into the local data list
void merge_region_and_tele_src(std::map<std::string, SrcRecInfo> &src_map,
                               std::map<std::string, SrcRecInfo> &rec_map,
                               std::map<std::string, std::map<std::string, std::vector<DataInfo>>> &data_map,
                               std::map<std::string, SrcRecInfo> &src_map_tele,
                               std::map<std::string, SrcRecInfo> &rec_map_tele,
                               std::map<std::string, std::map<std::string, std::vector<DataInfo>>> &data_map_tele){
    if(src_map_tele.size() > 0) {
        for (auto iter = src_map_tele.begin(); iter != src_map_tele.end(); iter++)
            src_map[iter->first] = iter->second;

        for (auto iter = rec_map_tele.begin(); iter != rec_map_tele.end(); iter++)
            rec_map[iter->first] = iter->second;

        data_map.insert(data_map_tele.begin(), data_map_tele.end());
    }

    // give new event id (accourding to the name)
    int id_count = 0;
    for (auto iter = src_map.begin(); iter != src_map.end(); iter ++){
        iter->second.id = id_count;
        id_count += 1;
    }
}


// distribute the source/receiver list and data list to all the processors
void distribute_src_rec_data(std::map<std::string, SrcRecInfo>&                                   src_map, \
                             std::map<std::string, SrcRecInfo>&                                   rec_map, \
                             std::map<std::string, std::map<std::string, std::vector<DataInfo>>>& data_map, \
                             std::vector<std::string>&                                            src_name_list, \
                             std::map<std::string, SrcRecInfo>&                                   src_map_this_sim, \
                             std::map<std::string, SrcRecInfo>&                                   rec_map_this_sim, \
                             std::map<std::string, std::map<std::string, std::vector<DataInfo>>>& data_map_this_sim, \
                             std::vector<std::string>&                                            src_name_list_this_sim){

    // number of total sources
    int n_src = 0;

    // number of total sources
    if (id_sim == 0 && subdom_main)
        n_src = src_map.size();

    // broadcast the number of sources to all the processors
    broadcast_i_single_inter_and_intra_sim(n_src, 0);

    // store the total number of sources
    nsrc_total = n_src;

    // assign sources to each simulutaneous run group
    for (int i_src = 0; i_src < n_src; i_src++) {

        // id of simulutaneous run group to which the i_src-th source belongs
        int dst_id_sim = select_id_sim_for_src(i_src, n_sims);

        // broadcast the source name
        std::string src_name;
        if (id_sim == 0 && subdom_main){
            src_name = src_name_list[i_src];
        }
        broadcast_str_inter_and_intra_sim(src_name, 0);

        if (id_sim==0){ // sender

            if (dst_id_sim == id_sim){ // this source belongs to this simulutaneous run group

                if (subdom_main){
                    // src
                    src_map_this_sim[src_name] = src_map[src_name];
                    // data
                    for (auto iter = data_map[src_name].begin(); iter != data_map[src_name].end(); iter++){
                        // rec by data
                        rec_map_this_sim[iter->first] = rec_map[iter->first];
                        for (auto& data : iter->second){
                            data_map_this_sim[src_name][iter->first].push_back(data);

                            // store the second receiver for rec_pair
                            if (data.is_rec_pair)
                                rec_map_this_sim[data.name_rec_pair[1]] = rec_map[data.name_rec_pair[1]];
                        }
                    }
                } // end if (subdom_main)

                // add the source name to the source name list
                src_name_list_this_sim.push_back(src_name);

            } else { // this source belongs to non-main simulutaneous run group
                if (subdom_main){
                    // send src
                    send_src_info_inter_sim(src_map[src_name], dst_id_sim);

                    // send number of receivers belonging to src
                    int n_rec = data_map[src_name].size();
                    send_i_single_sim(&n_rec, dst_id_sim);

                    if (n_rec > 0){

                        // send data_map[name_i_src] to the main process of dst_id_sim
                        for (auto iter = data_map[src_name].begin(); iter != data_map[src_name].end(); iter++){

                            // send rec_map[name_i_rec] to the main process of dst_id_sim
                            send_rec_info_inter_sim(rec_map[iter->first], dst_id_sim);

                            // send data_map[name_i_src].size() to the main process of dst_id_sim
                            int n_data = iter->second.size();
                            send_i_single_sim(&n_data, dst_id_sim);
                            // send data
                            for (auto& data : iter->second) {
                                send_data_info_inter_sim(data, dst_id_sim);

                                // send the second receiver for rec_pair
                                if (data.is_rec_pair)
                                    send_rec_info_inter_sim(rec_map[data.name_rec_pair[1]], dst_id_sim);
                            }

                        }
                    } // if (n_data > 0)
                } // end if (subdom_main)

            }
        } else { // receive
            if (dst_id_sim == id_sim){

                // receive src/rec_points from the main process of dst_id_sim
                if(subdom_main){
                    // prepare SrcRecInfo object for receiving the contents
                    SrcRecInfo tmp_SrcInfo;

                    // receive src
                    recv_src_info_inter_sim(tmp_SrcInfo, 0);

                    // add the received src_map to the src_map
                    src_map_this_sim[tmp_SrcInfo.name] = tmp_SrcInfo;

                    // recv number of receivers belonging to src
                    int n_rec = 0;
                    recv_i_single_sim(&n_rec, 0);

                    // receive data_info from the main process of dst_id_sim
                    for (int i_srcrec = 0; i_srcrec < n_rec; i_srcrec++){

                        // recv data_map[name_i_src].size() from the main process of dst_id_sim
                        SrcRecInfo tmp_RecInfo;
                        recv_rec_info_inter_sim(tmp_RecInfo, 0);

                        rec_map_this_sim[tmp_RecInfo.name] = tmp_RecInfo;

                        int n_data = 0;
                        recv_i_single_sim(&n_data, 0);

                        for (int i_data = 0; i_data < n_data; i_data++){
                            // prepare DataInfo object for receiving the contents
                            DataInfo tmp_DataInfo;
                            // receive data_info from the main process of dst_id_sim
                            recv_data_info_inter_sim(tmp_DataInfo, 0);
                            // add the received data_info to the data_info
                            data_map_this_sim[tmp_DataInfo.name_src][tmp_DataInfo.name_rec].push_back(tmp_DataInfo);

                            // store the second receiver for rec_pair
                            if (tmp_DataInfo.is_rec_pair)
                                recv_rec_info_inter_sim(rec_map_this_sim[tmp_DataInfo.name_rec_pair[1]], 0);
                        }

                    } // end of for i_srcrec
                } // end of if(subdom_main)

                // add the source name to the source name list
                src_name_list_this_sim.push_back(src_name);

            } else {
                // do nothing
            }

        } // end of if (id_sim==0)


    } // end of for i_src

    // check IP.src_ids_this_sim for this rank
    if (myrank==0 && if_verbose) {
        std::cout << id_sim << "assigned src id(name) : ";
        for (auto iter = src_map_this_sim.begin(); iter != src_map_this_sim.end(); iter++){
            std::cout << iter->second.id << "(" << iter->second.name << ") ";
        }
        std::cout << std::endl;
    }

}


void prepare_src_map_for_2d_solver(std::map<std::string, SrcRecInfo>& src_map_tele, \
                                   std::vector<std::string>& src_id2name_2d, \
                                   std::map<std::string, SrcRecInfo>& src_map_2d) {
    // src_map_tele: map of teleseismic src objects, (only the main process has the information.)
    // src_id2name_2d: list of src name assigned to this simultaneous run group
    // src_map_2d: src map assigned to this simultaneous run group

    std::map<std::string, SrcRecInfo> tmp_src_map_unique;
    std::vector<std::string> tmp_src_name_list_unique;

    // at first, make a depth-unique source list in the main process from src_map_tele
    if (id_sim==0 && subdom_main) {

        for (auto iter = src_map_tele.begin(); iter != src_map_tele.end(); iter++){
            std::string tmp_name = iter->second.name;

            // check if there is no element in tmp_src_map_unique with the same iter->second.depth
            bool if_unique = true;
            for (auto iter2 = tmp_src_map_unique.begin(); iter2 != tmp_src_map_unique.end(); iter2++){
                if (iter2->second.dep == iter->second.dep){
                    if_unique = false;
                    break;
                }
            }

            if (if_unique) {
                tmp_src_map_unique[tmp_name] = iter->second;
                tmp_src_name_list_unique.push_back(tmp_name);
            }
        }
    }

    // broadcast the number of unique sources to all processes
    int n_src_unique = 0;
    if (id_sim==0 && subdom_main) n_src_unique = tmp_src_map_unique.size();
    MPI_Bcast(&n_src_unique, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (subdom_main){
        // iterate over all the unique sources
        for (int i_src_unique = 0; i_src_unique < n_src_unique; i_src_unique++){
            int dst_id_sim = select_id_sim_for_src(i_src_unique, n_sims);

            if (id_sim==0){
                if (dst_id_sim==id_sim){
                    // store
                    src_map_2d[tmp_src_name_list_unique[i_src_unique]] = tmp_src_map_unique[tmp_src_name_list_unique[i_src_unique]];
                    src_id2name_2d.push_back(tmp_src_name_list_unique[i_src_unique]);
                } else {
                    // send to dst_id_sim
                    send_src_info_inter_sim(tmp_src_map_unique[tmp_src_name_list_unique[i_src_unique]], dst_id_sim);
                }
            } else {
                if (dst_id_sim==id_sim){
                    // receive from 0
                    SrcRecInfo tmp_src;
                    recv_src_info_inter_sim(tmp_src, 0);
                    src_map_2d[tmp_src.name] = tmp_src;
                    src_id2name_2d.push_back(tmp_src.name);
                } else {
                    // do nothing
                }
            }
        }
   } // end of if (subdom_main)

    // print the number of sources ssigned to this simultaneous run group
    for (int i_sim=0; i_sim < n_sims; i_sim++){
        if (id_sim==i_sim && subdom_main && id_subdomain==0){
            std::cout << "id_sim = " << id_sim << " : " << src_map_2d.size() << " 2d-sources assigned." << std::endl;
        }
        synchronize_all_world();
    }

}


//
// Belows are the function for send/receiving SrcRecInfo, DataInfo object
// so when you add the member in those class, it will be necessary to modify
// the functions below for sharing the new member.
//

void send_src_info_inter_sim(SrcRecInfo &src, int dest){

    send_i_single_sim(&src.id, dest);
    send_i_single_sim(&src.year , dest);
    send_i_single_sim(&src.month, dest);
    send_i_single_sim(&src.day  , dest);
    send_i_single_sim(&src.hour , dest);
    send_i_single_sim(&src.min  , dest);
    send_cr_single_sim(&src.sec, dest);
    send_cr_single_sim(&src.lat, dest);
    send_cr_single_sim(&src.lon, dest);
    send_cr_single_sim(&src.dep, dest);
    send_cr_single_sim(&src.mag, dest);
    send_i_single_sim(&src.n_data, dest);
    send_str_sim(src.name, dest);
    send_bool_single_sim(&src.is_out_of_region, dest);

}

void recv_src_info_inter_sim(SrcRecInfo &src, int orig){

    recv_i_single_sim(&src.id, orig);
    recv_i_single_sim(&src.year , orig);
    recv_i_single_sim(&src.month, orig);
    recv_i_single_sim(&src.day  , orig);
    recv_i_single_sim(&src.hour , orig);
    recv_i_single_sim(&src.min  , orig);
    recv_cr_single_sim(&src.sec, orig);
    recv_cr_single_sim(&src.lat, orig);
    recv_cr_single_sim(&src.lon, orig);
    recv_cr_single_sim(&src.dep, orig);
    recv_cr_single_sim(&src.mag, orig);
    recv_i_single_sim(&src.n_data, orig);
    recv_str_sim(src.name, orig);
    recv_bool_single_sim(&src.is_out_of_region, orig);
}


void send_rec_info_inter_sim(SrcRecInfo &rec, int dest){

    send_i_single_sim(&rec.id, dest);
    send_str_sim(rec.name, dest);
    send_cr_single_sim(&rec.lon, dest);
    send_cr_single_sim(&rec.lat, dest);
    send_cr_single_sim(&rec.dep, dest);
}


void recv_rec_info_inter_sim(SrcRecInfo &rec, int orig){

    recv_i_single_sim(&rec.id, orig);
    recv_str_sim(rec.name, orig);
    recv_cr_single_sim(&rec.lon, orig);
    recv_cr_single_sim(&rec.lat, orig);
    recv_cr_single_sim(&rec.dep, orig);
}

void send_data_info_inter_sim(DataInfo &data, int dest){

    send_cr_single_sim(&data.data_weight, dest);
    send_cr_single_sim(&data.weight, dest);
    send_str_sim(data.phase, dest);

    send_i_single_sim(&data.id_src, dest);
    send_str_sim(data.name_src, dest);
    send_i_single_sim(&data.id_rec, dest);
    send_str_sim(data.name_rec, dest);

    send_bool_single_sim(&data.is_src_rec,  dest);
    send_bool_single_sim(&data.is_rec_pair, dest);
    send_bool_single_sim(&data.is_src_pair, dest);

    if (data.is_src_rec){
       send_cr_single_sim(&data.travel_time_obs, dest);
    }

    if (data.is_rec_pair){
        //send_i_single_sim(&data.id_src_single, dest);
        //send_str_sim(data.name_src_single, dest);
        for (int ipair=0; ipair<2; ipair++){
            send_i_single_sim(&data.id_rec_pair[ipair], dest);
            send_str_sim(data.name_rec_pair[ipair], dest);
        }
        send_cr_single_sim(&data.cs_dif_travel_time_obs, dest);
    }

    if (data.is_src_pair){
        send_i_single_sim(&data.id_rec, dest);
        send_str_sim(data.name_rec, dest);
        for (int ipair=0; ipair<2; ipair++){
            send_i_single_sim(&data.id_src_pair[ipair], dest);
            send_str_sim(data.name_src_pair[ipair], dest);
        }
        send_cr_single_sim(&data.cr_dif_travel_time_obs, dest);
    }
}


void recv_data_info_inter_sim(DataInfo &data, int orig){

    recv_cr_single_sim(&data.data_weight, orig);
    recv_cr_single_sim(&data.weight, orig);
    recv_str_sim(data.phase, orig);

    recv_i_single_sim(&data.id_src, orig);
    recv_str_sim(data.name_src, orig);
    recv_i_single_sim(&data.id_rec, orig);
    recv_str_sim(data.name_rec, orig);

    recv_bool_single_sim(&data.is_src_rec, orig);
    recv_bool_single_sim(&data.is_rec_pair, orig);
    recv_bool_single_sim(&data.is_src_pair, orig);

    if (data.is_src_rec){
       recv_cr_single_sim(&data.travel_time_obs, orig);
    }

    if (data.is_rec_pair){
        //recv_i_single_sim(&data.id_src_single, orig);
        //recv_str_sim(data.name_src_single, orig);
        for (int ipair=0; ipair<2; ipair++){
            recv_i_single_sim(&data.id_rec_pair[ipair], orig);
            recv_str_sim(data.name_rec_pair[ipair], orig);
        }
        recv_cr_single_sim(&data.cs_dif_travel_time_obs, orig);
    }

    if (data.is_src_pair){
        recv_i_single_sim(&data.id_rec, orig);
        recv_str_sim(data.name_rec, orig);
        for (int ipair=0; ipair<2; ipair++){
            recv_i_single_sim(&data.id_src_pair[ipair], orig);
            recv_str_sim(data.name_src_pair[ipair], orig);
        }
        recv_cr_single_sim(&data.cr_dif_travel_time_obs, orig);
    }

}
