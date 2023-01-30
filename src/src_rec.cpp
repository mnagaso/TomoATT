#include "src_rec.h"


//
// functions for processing src_rec_file
//

void parse_src_rec_file(std::string& src_rec_file, \
                        std::vector<SrcRec>& src_points, \
                        std::vector<std::vector<SrcRec>>& rec_points){

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
    int cc =0; // count the number of lines
    int i_src_now = 0; // count the number of srcs
    int i_rec_now = 0; // count the number of receivers
    int ndata_tmp = 0; // count the number of receivers or differential traveltime data for each source
    std::vector<SrcRec> rec_points_tmp;
    src_points.clear();
    rec_points_tmp.clear();
    rec_points.clear();

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
                SrcRec src;
                //src.id_src     = std::stoi(tokens[0]);
                src.id_src     = i_src_now; // MNMN: here use id_src of active source lines order of src rec file, which allow to comment out bad events.
                src.id_src_ori = src.id_src; // use for swapping srcs
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
                    rec_points_tmp.push_back(rec);
                    cc++;
                    src_points.at(src_points.size()-1).n_rec++;

                } else {

                    // read differential traveltime
                    SrcRec rec;
                    //rec.id_src    = std::stoi(tokens[0]);
                    rec.id_src           = i_src_now; // MNMN: here use id_src of active source lines order of src rec file, which allow to comment out bad events.
                    rec.id_rec_pair[0]   = std::stoi(tokens[1]);
                    rec.name_rec_pair[0] = tokens[2];
                    rec.lat_pair[0]      = static_cast<CUSTOMREAL>(std::stod(tokens[3])); // in degree
                    rec.lon_pair[0]      = static_cast<CUSTOMREAL>(std::stod(tokens[4])); // in degree
                    rec.dep_pair[0]      = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[5])/1000.0); // convert elevation in meter to depth in km
                    rec.id_rec_pair[1]   = std::stoi(tokens[6]);
                    rec.name_rec_pair[1] = tokens[7];
                    rec.lat_pair[1]      = static_cast<CUSTOMREAL>(std::stod(tokens[8])); // in degree
                    rec.lon_pair[1]      = static_cast<CUSTOMREAL>(std::stod(tokens[9])); // in degree
                    rec.dep_pair[1]      = static_cast<CUSTOMREAL>(-1.0*std::stod(tokens[10])/1000.0); // convert elevation in meter to depth in km
                    rec.phase            = tokens[11];
                    // rec.dif_arr_time = static_cast<CUSTOMREAL>(std::stod(tokens[12]));
                    rec.dif_arr_time     = 0.0;
                    rec.dif_arr_time_ori = static_cast<CUSTOMREAL>(std::stod(tokens[12])); // store read data

                    // check if tokens[9] exists read weight
                    if (tokens.size() > 13)
                        rec.weight = static_cast<CUSTOMREAL>(std::stod(tokens[13]));
                    else
                        rec.weight = 1.0; // default weight
                    rec.is_rec_pair = true;
                    rec_points_tmp.push_back(rec);
                    cc++;
                    src_points.at(src_points.size()-1).n_rec_pair++;

                }

                if (cc > ndata_tmp) {
                    // go to the next source
                    cc = 0;
                    rec_points.push_back(rec_points_tmp);
                    rec_points_tmp.clear();
                    i_src_now++;
                    i_rec_now = 0;

                    // timer
                    if (i_src_now % 100 == 0 && world_rank == 0) {
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
    int n_src_points = src_points.size();
    if (n_src_points < n_sims){
        std::cout << "Error: number of sources in src_rec_file is less than n_sims. Abort." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // indicate elapsed time
    std::cout << "Total elapsed time for reading src_rec_file: " << timer.get_t() << " seconds.\n";
}


void do_swap_src_rec(std::vector<SrcRec>& src_points, std::vector<std::vector<SrcRec>>& rec_points, \
                     std::vector<SrcRec>& src_points_back, std::vector<std::vector<SrcRec>>& rec_points_back){

    // swap src/rec points
    // at this moment, all the sources are divided into src_points (regional) and tele_src_points (teleseismic)

    // Start timer
    std::string timer_name = "swap_src_rec";
    Timer timer(timer_name);

    std::vector<SrcRec> new_srcs; // new src points
    std::vector<std::vector<SrcRec>> new_recs;

    // generate new source list
    for (long unsigned int i_src = 0; i_src < src_points.size(); i_src++) {
        for(long unsigned int i_rec = 0; i_rec < rec_points[i_src].size(); i_rec++) {

            if (new_srcs.size() == 0){
                // add receiver to the new source list
                new_srcs.push_back(rec_points[i_src][i_rec]);
                // record the original source id of this receiver
                new_srcs.back().id_srcs_ori.push_back(rec_points[i_src][i_rec].id_src);
            } else if (new_srcs.size() != 0) {
                // check if the new source list already has the same receiver
                bool found = false;
                for (long unsigned int i_new_src = 0; i_new_src < new_srcs.size(); i_new_src++) {
                    if (new_srcs[i_new_src].name_rec.compare(rec_points[i_src][i_rec].name_rec) == 0) {
                        // add the original source id of this receiver
                        new_srcs[i_new_src].id_srcs_ori.push_back(rec_points[i_src][i_rec].id_src);
                        found=true;
                        break;
                    }
                    else {
                    }
                }
                // if not found, add the receiver to the new source list
                if (!found) {
                    // add receiver to the new source list
                    new_srcs.push_back(rec_points[i_src][i_rec]);
                    // record the original source id of this receiver
                    new_srcs.back().id_srcs_ori.push_back(rec_points[i_src][i_rec].id_src);
                }

            }
        }
    }

    // generate new rec list
    for (long unsigned int i_src = 0; i_src < new_srcs.size(); i_src++) {

        // set new id_src
        new_srcs[i_src].id_src = i_src;

        std::vector<SrcRec> tmp_list_recs;

        for (auto& i_src_ori : new_srcs[i_src].id_srcs_ori){
            // copy the original source object for temporal use
            SrcRec tmp_new_rec = src_points_back[i_src_ori];

            // loop over all the receivers of the original source
            for (auto& tmp_rec_ori : rec_points_back[i_src_ori]){
                // check if the receiver is the same as the new source
                if (tmp_rec_ori.name_rec.compare(new_srcs[i_src].name_rec)==0) {
                    // we can use the same arrival time for a src-rec pair by principle of reciprocity (Aki & Richards, 2002).
                    tmp_new_rec.arr_time     = tmp_rec_ori.arr_time;
                    tmp_new_rec.arr_time_ori = tmp_rec_ori.arr_time_ori;
                    tmp_new_rec.id_rec_ori   = tmp_rec_ori.id_rec;
                    tmp_new_rec.id_src       = i_src; // over write the id_src of the original source for the new source
                    goto rec_found;
                }
            }

            rec_found:
                tmp_list_recs.push_back(tmp_new_rec);
        }

        new_recs.push_back(tmp_list_recs);

        // set n_rec
        new_srcs[i_src].n_rec = new_recs[i_src].size();

    }


    // backup and set new src/rec points
    // !! ONLY REGIONAL EVENTS AND RECEIVERS ARE STORED IN *_back VECTORS !!
    src_points.clear();
    rec_points.clear();
    // teleseismic events are concatenate to the vectors below later
    src_points = new_srcs;
    rec_points = new_recs;

    // indicate elapsed time
    std::cout << "Total elapsed time for swapping src rec: " << timer.get_t() << " seconds.\n";
}


void reverse_src_rec_points(std::vector<SrcRec>& src_points,      std::vector<std::vector<SrcRec>>& rec_points, \
                            std::vector<SrcRec>& src_points_back, std::vector<std::vector<SrcRec>>& rec_points_back, \
                            std::vector<SrcRec>& src_points_out,  std::vector<std::vector<SrcRec>>& rec_points_out, \
                            const bool& swap_src_rec, const int& run_mode){

    // loop swapped sources
    for (long unsigned int i_src = 0; i_src < src_points.size(); i_src++){
        // swapped only the regional events && really need swap
        if (src_points[i_src].is_teleseismic == false && swap_src_rec){

            // int id_rec_orig = src_points[i_src].id_rec;
            // loop swapped receivers
            for (long unsigned int i_rec = 0; i_rec < rec_points[i_src].size(); i_rec++){
                int id_src_orig = rec_points[i_src][i_rec].id_src_ori;
                int id_rec_orig = rec_points[i_src][i_rec].id_rec_ori; // cannot fully recover the original receiver id

                // store calculated arrival time in backuped receiver list
                rec_points_back[id_src_orig][id_rec_orig].arr_time     = rec_points[i_src][i_rec].arr_time; ///////////
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


void writeout_src_rec_file(std::string&                      src_rec_file_out, \
                           std::vector<SrcRec>&              src_points_out, \
                           std::vector<std::vector<SrcRec>>& rec_points_out){

    std::ofstream ofs;

    for (long unsigned int i_src = 0; i_src < src_points_out.size(); i_src++){
        if (world_rank == 0 && id_subdomain==0){    // main processor of subdomain && the first id of subdoumains
            if (i_src == 0)
                ofs.open(src_rec_file_out);
            else
                ofs.open(src_rec_file_out, std::ios_base::app);

            // set output precision
            // ofs << std::fixed << std::setprecision(ASCII_OUTPUT_PRECISION);

            // format should be the same as input src_rec_file
            // source line :  id_src yearm month day hour min sec lat lon dep_km mag num_recs id_event weight
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
                    // receiver line : id_src id_rec name_rec lat lon elevation_m phase epicentral_distance_km arival_time weight
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
    }



}