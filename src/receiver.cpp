#include "receiver.h"


Receiver::Receiver() {
}


Receiver::~Receiver() {
}


void Receiver::interpolate_and_store_arrival_times_at_rec_position(InputParams& IP, Grid& grid, const std::string& name_sim_src) {
    if(subdom_main){
        // get receiver positions from input parameters
        std::vector<std::string> name_receivers = IP.get_rec_names(name_sim_src);

        // calculate the travel time of the receiver by interpolation
        for(std::string name_rec: name_receivers) {
            IP.data_map[name_sim_src][name_rec].travel_time = interpolate_travel_time(grid, IP, name_sim_src, name_rec);
        }
    }
}


std::vector<CUSTOMREAL> Receiver::calculate_adjoint_source(InputParams& IP, const std::string& name_sim_src) {

    CUSTOMREAL obj           = 0.0;
    CUSTOMREAL obj_abs       = 0.0;
    CUSTOMREAL obj_cs_dif    = 0.0;
    CUSTOMREAL obj_cr_dif    = 0.0;
    CUSTOMREAL obj_tele      = 0.0;
    CUSTOMREAL misfit        = 0.0;
    CUSTOMREAL misfit_abs    = 0.0;
    CUSTOMREAL misfit_cs_dif = 0.0;
    CUSTOMREAL misfit_cr_dif = 0.0;
    CUSTOMREAL misfit_tele   = 0.0;

    // obj, obj_abs, obj_cs_dif, obj_cr_dif, obj_tele, misfit, misfit_abs, misfit_cs_dif, misfit_cr_dif, misfit_tele
    std::vector<CUSTOMREAL> allsum = std::vector<CUSTOMREAL>(20);

    // rec.adjoint_source = 0
    IP.initialize_adjoint_source();

    int count = 0;

    if(subdom_main){

        // loop all data related to this source MNMN: use reference(auto&) to avoid copy
        for (auto it_src = IP.data_map[name_sim_src].begin(); it_src != IP.data_map[name_sim_src].end(); ++it_src) {
            // reference to the target data
            DataInfo& data = it_src->second;

            //
            // absolute traveltime
            //
            if (data.is_src_rec){

                // skip if this data is not belongs to this simulation group
                // (this happens only the first group, because the first group retains all data)
                if (data.name_src != name_sim_src) continue;

                std::string name_src      = data.name_src;
                std::string name_rec      = data.name_rec;
                CUSTOMREAL syn_time       = data.travel_time;
                CUSTOMREAL obs_time       = data.travel_time_obs;
                CUSTOMREAL adjoint_source = IP.get_rec_point(name_rec).adjoint_source + (syn_time - obs_time)*data.weight;
                IP.set_adjoint_source(name_rec, adjoint_source); // set adjoint source to rec_map[name_rec]

                // contribute misfit
                obj     += 1.0 * my_square(syn_time - obs_time)*data.weight;
                misfit  += 1.0 * my_square(syn_time - obs_time);

                if (IP.get_src_point(name_src).is_out_of_region || IP.get_rec_point(name_rec).is_out_of_region){
                    obj_tele        +=  1.0 * my_square(syn_time - obs_time)*data.weight;
                    misfit_tele     +=  1.0 * my_square(syn_time - obs_time);
                } else{
                    obj_abs         +=  1.0 * my_square(syn_time - obs_time)*data.weight;
                    misfit_abs      +=  1.0 * my_square(syn_time - obs_time);
                }

            //
            // common receiver differential traveltime
            //
            } else if (data.is_src_pair) {

                std::string name_src1 = data.name_src_pair[0];
                std::string name_src2 = data.name_src_pair[1];
                std::string name_rec  = data.name_rec_single;

                // skip if this data is not belongs to this simulation group
                // (this happens only the first group, because the first group retains all data)
                if (name_sim_src != name_src1 && name_sim_src != name_src2) continue;

                if (name_sim_src == name_src1){

                    CUSTOMREAL syn_dif_time   = IP.data_map[name_src1][name_rec].travel_time - IP.data_map[name_src2][name_rec].travel_time;
                    CUSTOMREAL obs_dif_time   = data.cr_dif_travel_time_obs;
                    CUSTOMREAL adjoint_source = IP.get_rec_point(name_rec).adjoint_source + (syn_dif_time - obs_dif_time)*data.weight;

                    IP.set_adjoint_source(name_rec, adjoint_source);

                    // contribute misfit
                    obj     += 0.5 * my_square(syn_dif_time - obs_dif_time)*data.weight;
                    misfit  += 0.5 * my_square(syn_dif_time - obs_dif_time);

                    if (IP.get_src_point(name_src1).is_out_of_region || \
                        IP.get_src_point(name_src2).is_out_of_region || \
                        IP.get_rec_point(name_rec).is_out_of_region){
                        obj_tele        += 0.5 * my_square(syn_dif_time - obs_dif_time)*data.weight;
                        misfit_tele     += 0.5 * my_square(syn_dif_time - obs_dif_time);   // because a pair sf sources are counted twice, thus * 0.5
                    } else{
                        obj_cr_dif      += 0.5 * my_square(syn_dif_time - obs_dif_time)*data.weight;
                        misfit_cr_dif   += 0.5 * my_square(syn_dif_time - obs_dif_time);
                    }

                } else if (name_sim_src == name_src2) {

                    CUSTOMREAL syn_dif_time   = IP.data_map[name_src2][name_rec].travel_time - IP.data_map[name_src1][name_rec].travel_time;
                    CUSTOMREAL obs_dif_time   = - data.cr_dif_travel_time_obs;
                    CUSTOMREAL adjoint_source = IP.get_rec_point(name_rec).adjoint_source + (syn_dif_time - obs_dif_time)*data.weight;

                    IP.set_adjoint_source(name_rec, adjoint_source);

                    // contribute misfit
                    obj     += 0.5 * my_square(syn_dif_time - obs_dif_time)*data.weight;
                    misfit  += 0.5 * my_square(syn_dif_time - obs_dif_time);

                    if (IP.get_src_point(name_src1).is_out_of_region || \
                        IP.get_src_point(name_src2).is_out_of_region || \
                        IP.get_rec_point(name_rec).is_out_of_region){
                        obj_tele        +=  0.5 * my_square(syn_dif_time - obs_dif_time)*data.weight;
                        misfit_tele     +=  0.5 * my_square(syn_dif_time - obs_dif_time);
                    } else{
                        obj_cr_dif      +=  0.5 * my_square(syn_dif_time - obs_dif_time)*data.weight;
                        misfit_cr_dif   +=  0.5 * my_square(syn_dif_time - obs_dif_time);
                    }

                } else {
                    std::cout << "error match of data in function: calculate_adjoint_source() " << std::endl;
                }

            //
            // common source differential traveltime
            //
            } else if (data.is_rec_pair) {

                std::string name_src  = data.name_src_single;
                std::string name_rec1 = data.name_rec_pair[0];
                std::string name_rec2 = data.name_rec_pair[1];

                CUSTOMREAL syn_dif_time = IP.data_map[name_src][name_rec1].travel_time - IP.data_map[name_src][name_rec2].travel_time;
                CUSTOMREAL obs_dif_time = data.cs_dif_travel_time_obs;

                CUSTOMREAL adjoint_source = IP.get_rec_point(name_rec1).adjoint_source + (syn_dif_time - obs_dif_time)*data.weight;
                IP.set_adjoint_source(name_rec1, adjoint_source);

                adjoint_source = IP.get_rec_point(name_rec2).adjoint_source - (syn_dif_time - obs_dif_time)*data.weight;
                IP.set_adjoint_source(name_rec2, adjoint_source);

                // contribute misfit
                obj     += 1.0 * my_square(syn_dif_time - obs_dif_time)*data.weight;
                misfit  += 1.0 * my_square(syn_dif_time - obs_dif_time);

                if (IP.get_src_point(name_src).is_out_of_region || \
                    IP.get_rec_point(name_rec1).is_out_of_region || \
                    IP.get_rec_point(name_rec2).is_out_of_region){
                    obj_tele        += 1.0 * my_square(syn_dif_time - obs_dif_time)*data.weight;
                    misfit_tele     += 1.0 * my_square(syn_dif_time - obs_dif_time);
                } else{
                    obj_cs_dif      += 1.0 * my_square(syn_dif_time - obs_dif_time)*data.weight;
                    misfit_cs_dif   += 1.0 * my_square(syn_dif_time - obs_dif_time);
                }
            }
            count += 1;

        }
    }


    // share the calculated objective function value to other processes
    broadcast_cr_single_sub(obj,0);
    broadcast_cr_single_sub(obj_abs,0);
    broadcast_cr_single_sub(obj_cs_dif,0);
    broadcast_cr_single_sub(obj_cr_dif,0);
    broadcast_cr_single_sub(obj_tele,0);
    broadcast_cr_single_sub(misfit,0);
    broadcast_cr_single_sub(misfit_abs,0);
    broadcast_cr_single_sub(misfit_cs_dif,0);
    broadcast_cr_single_sub(misfit_cr_dif,0);
    broadcast_cr_single_sub(misfit_tele,0);

    allsum[0] = obj;
    allsum[1] = obj_abs;
    allsum[2] = obj_cs_dif;
    allsum[3] = obj_cr_dif;
    allsum[4] = obj_tele;
    allsum[5] = misfit;
    allsum[6] = misfit_abs;
    allsum[7] = misfit_cs_dif;
    allsum[8] = misfit_cr_dif;
    allsum[9] = misfit_tele;

    return allsum;

}


CUSTOMREAL Receiver::interpolate_travel_time(Grid& grid, InputParams& IP, std::string name_src, std::string name_rec) {
    // calculate the travel time of the receiver by 3d linear interpolation

    // get the reference for a receiver
    const SrcRecInfo& rec = IP.get_rec_point(name_rec);

    // copy some parameters
    CUSTOMREAL delta_lon = grid.get_delta_lon();
    CUSTOMREAL delta_lat = grid.get_delta_lat();
    CUSTOMREAL delta_r   = grid.get_delta_r();

    // store receiver position in radian
    CUSTOMREAL rec_lon = rec.lon*DEG2RAD;
    CUSTOMREAL rec_lat = rec.lat*DEG2RAD;
    CUSTOMREAL rec_r = depth2radius(rec.dep); // r in km

    // check if the receiver is in this subdomain
    bool is_in_subdomain = false;
    if (grid.get_lon_min_loc() <= rec_lon && rec_lon < grid.get_lon_max_loc() && \
        grid.get_lat_min_loc() <= rec_lat && rec_lat < grid.get_lat_max_loc() && \
        grid.get_r_min_loc()   <= rec_r   && rec_r   < grid.get_r_max_loc()   ) {
        is_in_subdomain = true;
    }

    // check the rank where the source is located
    int rec_rank = -1;
    bool* rec_flags = new bool[nprocs];
    allgather_bool_single(&is_in_subdomain, rec_flags);
    for (int i = 0; i < nprocs; i++) {
        if (rec_flags[i]) {
            rec_rank = i;
            break; // this break means that the first subdomain is used if the receiver is in multiple subdomains (ghost layer)
        }
    }
    delete[] rec_flags;
     // check if the receiver is in the global domain
    if (rec_rank == -1) {
        std::cout << "Error: the receiver is not in the global domain" << std::endl;
        // print rec
        std::cout << "name_src: " << name_src << " name_rec: " << name_rec << " depth: " << rec.dep << " lat: " << rec.lat << " lon: " << rec.lon << std::endl;
        // print boundary
        //std::cout << "lon min max rec: " << grid.get_lon_min_loc() << " " << grid.get_lon_max_loc() << " " << rec_lon << std::endl;
        //std::cout << "lat min max rec: " << grid.get_lat_min_loc() << " " << grid.get_lat_max_loc() << " " << rec_lat << std::endl;
        //std::cout << "r min max rec: " << grid.get_r_min_loc() << " " << grid.get_r_max_loc() << " " << rec_r << std::endl;

        std::cout << "lon+bound min max rec: " << (grid.get_lon_min_loc() - delta_lon)*RAD2DEG     << " " << (grid.get_lon_max_loc() + delta_lon)*RAD2DEG     << " " << rec_lon*RAD2DEG    << std::endl;
        std::cout << "lat+bound min max rec: " << (grid.get_lat_min_loc() - delta_lat)*RAD2DEG     << " " << (grid.get_lat_max_loc() + delta_lat)*RAD2DEG     << " " << rec_lat*RAD2DEG    << std::endl;
        std::cout << "r+bound min max rec: "   << radius2depth(grid.get_r_min_loc()   - delta_r  ) << " " << radius2depth(grid.get_r_max_loc()   + delta_r  ) << " " << radius2depth(rec_r)<< std::endl;
        exit(1);
    }

    CUSTOMREAL vinterp = 0.0;

    if (is_in_subdomain) {
        // calculate the interpolated travel time and broadcast it

        // descretize source position (LOCAL) ID)
        int i_rec = std::floor((rec_lon - grid.get_lon_min_loc()) / delta_lon);
        int j_rec = std::floor((rec_lat - grid.get_lat_min_loc()) / delta_lat);
        int k_rec = std::floor((rec_r   - grid.get_r_min_loc())   / delta_r);

        // discretized receiver position
        CUSTOMREAL dis_rec_lon = grid.p_loc_1d[i_rec];
        CUSTOMREAL dis_rec_lat = grid.t_loc_1d[j_rec];
        CUSTOMREAL dis_rec_r   = grid.r_loc_1d[k_rec];

        // relative position errors
        CUSTOMREAL e_lon = std::min({_1_CR,(rec_lon - dis_rec_lon)/delta_lon});
        CUSTOMREAL e_lat = std::min({_1_CR,(rec_lat - dis_rec_lat)/delta_lat});
        CUSTOMREAL e_r   = std::min({_1_CR,(rec_r   - dis_rec_r)  /delta_r});

//        if(if_verbose){
//            std::cout << "(rec_lon - dis_rec_lon)/dlon: " << (rec_lon - dis_rec_lon)/delta_lon << std::endl;
//            std::cout << "(rec_lat - dis_rec_lat)/dlat: " << (rec_lat - dis_rec_lat)/delta_lat << std::endl;
//            std::cout << "(rec_r   - dis_rec_r  )/dr: "   << (rec_r   - dis_rec_r  )/delta_r << std::endl;
//            std::cout << "rec_lon, dis_rec_lon, delta_lon: " << rec_lon << " " << dis_rec_lon << " " << delta_lon << std::endl;
//            std::cout << "rec_lat, dis_rec_lat, delta_lat: " << rec_lat << " " << dis_rec_lat << " " << delta_lat << std::endl;
//            std::cout << "rec_r, dis_rec_r, delta_r  : " << rec_r << ", " << dis_rec_r << ", " <<  delta_r  << std::endl;
//            std::cout << "loc_K. k_rec: " << loc_K << "," << k_rec << std::endl;
//            std::cout << "r_loc_1d[loc_K-1]: " << grid.r_loc_1d[loc_K-1] << std::endl;
//        }

        int i_rec_p1 = i_rec + 1;
        int j_rec_p1 = j_rec + 1;
        int k_rec_p1 = k_rec + 1;

        // exclude the points if they are out of the domain
        if (i_rec_p1 > loc_I-1 \
         || j_rec_p1 > loc_J-1 \
         || k_rec_p1 > loc_K-1) {
            // exit(1) as the source is out of the domain
            std::cout << "Error: the receiver is out of the domain" << std::endl;
            std::cout << "name_src: " << name_src << " name_rec: " << name_rec << " depth: " << rec.dep << " lat: " << rec.lat << " lon: " << rec.lon << std::endl;
            std::cout << "lon min max rec: " << grid.get_lon_min_loc()*RAD2DEG << " " << grid.get_lon_max_loc()*RAD2DEG << " " << rec_lon*RAD2DEG << std::endl;
            std::cout << "lat min max rec: " << grid.get_lat_min_loc()*RAD2DEG << " " << grid.get_lat_max_loc()*RAD2DEG << " " << rec_lat*RAD2DEG << std::endl;
            std::cout << "r min max rec: " << radius2depth(grid.get_r_min_loc()) << " " << radius2depth(grid.get_r_max_loc()) << " " << radius2depth(rec_r) << std::endl;
            std::cout << "i_rec: " << i_rec << " j_rec: " << j_rec << " k_rec: " << k_rec << std::endl;
            std::cout << "i_rec_p1: " << i_rec_p1 << " j_rec_p1: " << j_rec_p1 << " k_rec_p1: " << k_rec_p1 << std::endl;
            std::cout << "loc_I-1: " << loc_I-1 << " loc_J-1: " << loc_J-1 << " loc_K-1: " << loc_K-1 << std::endl;
            //finalize_mpi();
            exit(1);
         }


        vinterp = (_1_CR - e_lon) * (_1_CR - e_lat) * (_1_CR - e_r) * grid.T_loc[I2V(i_rec,   j_rec,   k_rec)]   \
                +          e_lon  * (_1_CR - e_lat) * (_1_CR - e_r) * grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec)]   \
                + (_1_CR - e_lon) *          e_lat  * (_1_CR - e_r) * grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec)]   \
                + (_1_CR - e_lon) * (_1_CR - e_lat) *          e_r  * grid.T_loc[I2V(i_rec,   j_rec,   k_rec_p1)] \
                +          e_lon  *          e_lat  * (_1_CR - e_r) * grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec)]   \
                +          e_lon  * (_1_CR - e_lat) *          e_r  * grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec_p1)] \
                + (_1_CR - e_lon) *          e_lat  *          e_r  * grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec_p1)] \
                +          e_lon  *          e_lat  *          e_r  * grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec_p1)];

        //std::cout << "DEBUG near and vinterp : " << grid.T_loc[I2V(i_rec,j_rec,k_rec)] << ", " << vinterp << std::endl;
        // std::cout << "times: " << grid.T_loc[I2V(i_rec,   j_rec,   k_rec)] << ", "
        //           << grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec)] << ", "
        //           << grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec)] << ", "
        //           << grid.T_loc[I2V(i_rec,   j_rec,   k_rec_p1)] << ", "
        //           << grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec)] << ", "
        //           << grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec_p1)] << ", "
        //           << grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec_p1)] << ", "
        //           << grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec_p1)] << ", "
        //           << std::endl;
        // std::cout << "rec: " << rec.name << ", lat: " << rec_lat
        //           << ", lon: " << rec_lon << ", dep: " << rec.dep
        //           << ", time: " << vinterp
        //           << std::endl;

        // broadcast interpolated travel time
        broadcast_cr_single(vinterp, rec_rank);

    } else {
        // receive the calculated traveltime
        broadcast_cr_single(vinterp, rec_rank);
    }

    // return the calculated travel time
    return vinterp;
}



void Receiver::init_vars_src_reloc(InputParams& IP){
    if (subdom_main) {

        // calculate gradient of travel time at each receiver (swapped source)
        for (auto iter = IP.rec_map.begin(); iter != IP.rec_map.end(); iter++) {
            iter->second.tau_opt                    = _0_CR;
            iter->second.grad_chi_i                 = _0_CR;
            iter->second.grad_chi_j                 = _0_CR;
            iter->second.grad_chi_k                 = _0_CR;
            iter->second.sum_weight                 = _0_CR;
            iter->second.vobj_src_reloc_old         = iter->second.vobj_src_reloc;
            iter->second.vobj_src_reloc             = _0_CR;
            iter->second.vobj_grad_norm_src_reloc   = _0_CR;
        }

    }
}

void Receiver::calculate_T_gradient(InputParams& IP, Grid& grid, const std::string& name_sim_src){

    if(subdom_main){
        // calculate gradient of travel time at each receiver (swapped source)
        for (auto iter = IP.data_map[name_sim_src].begin(); iter != IP.data_map[name_sim_src].end(); iter++){

            std::string name_rec = iter->first;
            calculate_T_gradient_one_rec(grid, IP.rec_map[name_rec]);

        // for(auto iter = IP.rec_map.begin(); iter != IP.rec_map.end(); iter++){
        //     std::cout << "DTi (lon) is: " << iter->second.DTi << std::endl;
        //     std::cout << "DTj (lat) is: " << iter->second.DTj << std::endl;
        //     std::cout << "DTk (dep) is: " << iter->second.DTk << std::endl;
        // }
        }
    }

}


void Receiver::calculate_T_gradient_one_rec(Grid& grid, SrcRecInfo& rec){

    // calculate the travel time of the receiver by 3d linear interpolation

    // copy some parameters
    CUSTOMREAL delta_lon = grid.get_delta_lon();
    CUSTOMREAL delta_lat = grid.get_delta_lat();
    CUSTOMREAL delta_r   = grid.get_delta_r();

    // store receiver position in radian
    CUSTOMREAL rec_lon = rec.lon*DEG2RAD;
    CUSTOMREAL rec_lat = rec.lat*DEG2RAD;
    CUSTOMREAL rec_r = depth2radius(rec.dep); // r in km

    // check if the receiver is in this subdomain
    bool is_in_subdomain = false;
    if (grid.get_lon_min_loc() <= rec_lon && rec_lon < grid.get_lon_max_loc() && \
        grid.get_lat_min_loc() <= rec_lat && rec_lat < grid.get_lat_max_loc() && \
        grid.get_r_min_loc()   <= rec_r   && rec_r   < grid.get_r_max_loc()   ) {
        is_in_subdomain = true;
    }

    // check the rank where the source is located
    int rec_rank = -1;
    bool* rec_flags = new bool[nprocs];
    allgather_bool_single(&is_in_subdomain, rec_flags);
    for (int i = 0; i < nprocs; i++) {
        if (rec_flags[i]) {
            rec_rank = i;
            break; // this break means that the first subdomain is used if the receiver is in multiple subdomains (ghost layer)
        }
    }
    delete[] rec_flags;
     // check if the receiver is in the global domain
    if (rec_rank == -1) {
        std::cout << "Error: the receiver is not in the global domain" << std::endl;
        // print rec
        std::cout << " id_rec: " << rec.id << " name: " << rec.name << " depth: " << rec.dep << " lat: " << rec.lat << " lon: " << rec.lon << std::endl;
        std::cout << "lon+bound min max rec: " << (grid.get_lon_min_loc() - delta_lon)*RAD2DEG     << " " << (grid.get_lon_max_loc() + delta_lon)*RAD2DEG     << " " << rec_lon*RAD2DEG    << std::endl;
        std::cout << "lat+bound min max rec: " << (grid.get_lat_min_loc() - delta_lat)*RAD2DEG     << " " << (grid.get_lat_max_loc() + delta_lat)*RAD2DEG     << " " << rec_lat*RAD2DEG    << std::endl;
        std::cout << "r+bound min max rec: "   << radius2depth(grid.get_r_min_loc()   - delta_r  ) << " " << radius2depth(grid.get_r_max_loc()   + delta_r  ) << " " << radius2depth(rec_r)<< std::endl;
        exit(1);
    }

    CUSTOMREAL DTi = 0.0, DTj = 0.0, DTk = 0.0;

    if (is_in_subdomain) {
        // calculate the interpolated travel time and broadcast it

        // descretize source position (LOCAL) ID)
        int i_rec = std::floor((rec_lon - grid.get_lon_min_loc()) / delta_lon);
        int j_rec = std::floor((rec_lat - grid.get_lat_min_loc()) / delta_lat);
        int k_rec = std::floor((rec_r   - grid.get_r_min_loc())   / delta_r);

        // std::cout << "lon: " << rec.lon << ", lat: " << rec.lat << ", dep: " << rec.dep << std::endl;
        // std::cout << "i_rec: " << i_rec << ", j_rec: " << j_rec << ", k_rec: " << k_rec << std::endl;

        // discretized receiver position
        CUSTOMREAL dis_rec_lon = grid.p_loc_1d[i_rec];
        CUSTOMREAL dis_rec_lat = grid.t_loc_1d[j_rec];
        CUSTOMREAL dis_rec_r   = grid.r_loc_1d[k_rec];

        // relative position errors
        CUSTOMREAL e_lon = std::min({_1_CR,(rec_lon - dis_rec_lon)/delta_lon});
        CUSTOMREAL e_lat = std::min({_1_CR,(rec_lat - dis_rec_lat)/delta_lat});
        CUSTOMREAL e_r   = std::min({_1_CR,(rec_r   - dis_rec_r)  /delta_r});

        int i_rec_p1 = i_rec + 1;
        int j_rec_p1 = j_rec + 1;
        int k_rec_p1 = k_rec + 1;


// TODO: MNMN : temporally comment out because this condition hits on the boundary of subdomains
/*        // exclude the points if they are adjancet to the domain boundary
        if (i_rec_p1 > loc_I-2 || j_rec_p1 > loc_J-2 || k_rec_p1 > loc_K-2 || \
            i_rec    < 1       || j_rec    < 1       || k_rec    < 1 ) {
            // exit(1) as the source is out of the domain
            std::cout << "Error: the location is too close to the domain boundary" << std::endl;
            std::cout << " id_rec: " << rec.id << " name: " << rec.name << " depth: " << rec.dep << " lat: " << rec.lat << " lon: " << rec.lon << std::endl;
            std::cout << "lon min max rec: " << grid.get_lon_min_loc()*RAD2DEG << " " << grid.get_lon_max_loc()*RAD2DEG << " " << rec_lon*RAD2DEG << std::endl;
            std::cout << "lat min max rec: " << grid.get_lat_min_loc()*RAD2DEG << " " << grid.get_lat_max_loc()*RAD2DEG << " " << rec_lat*RAD2DEG << std::endl;
            std::cout << "r min max rec: " << radius2depth(grid.get_r_min_loc()) << " " << radius2depth(grid.get_r_max_loc()) << " " << radius2depth(rec_r) << std::endl;
            std::cout << "i_rec: " << i_rec << " j_rec: " << j_rec << " k_rec: " << k_rec << std::endl;
            std::cout << "i_rec_p1: " << i_rec_p1 << " j_rec_p1: " << j_rec_p1 << " k_rec_p1: " << k_rec_p1 << std::endl;
            std::cout << "loc_I-1: " << loc_I-1 << " loc_J-1: " << loc_J-1 << " loc_K-1: " << loc_K-1 << std::endl;
            //finalize_mpi();
            exit(1);
        }
*/
        // CUSTOMREAL Ti, Tip, Tj, Tjp, Tk, Tkp;
        // Tk =      (_1_CR - e_lon) * (_1_CR - e_lat) * _1_CR         * grid.T_loc[I2V(i_rec,   j_rec,   k_rec)]
        //         +          e_lon  * (_1_CR - e_lat) * _1_CR         * grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec)]
        //         + (_1_CR - e_lon) *          e_lat  * _1_CR         * grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec)]
        //         +          e_lon  *          e_lat  * _1_CR         * grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec)];

        // Tkp =     (_1_CR - e_lon) * (_1_CR - e_lat) * _1_CR         * grid.T_loc[I2V(i_rec,   j_rec,   k_rec_p1)]
        //         +          e_lon  * (_1_CR - e_lat) * _1_CR         * grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec_p1)]
        //         + (_1_CR - e_lon) *          e_lat  * _1_CR         * grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec_p1)]
        //         +          e_lon  *          e_lat  * _1_CR         * grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec_p1)];

        // Tj =      (_1_CR - e_lon) * _1_CR           * (_1_CR - e_r) * grid.T_loc[I2V(i_rec,   j_rec,   k_rec)]
        //         +          e_lon  * _1_CR           * (_1_CR - e_r) * grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec)]
        //         + (_1_CR - e_lon) * _1_CR           *          e_r  * grid.T_loc[I2V(i_rec,   j_rec,   k_rec_p1)]
        //         +          e_lon  * _1_CR           *          e_r  * grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec_p1)];

        // Tjp =     (_1_CR - e_lon) * _1_CR           * (_1_CR - e_r) * grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec)]
        //         +          e_lon  * _1_CR           * (_1_CR - e_r) * grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec)]
        //         + (_1_CR - e_lon) * _1_CR           *          e_r  * grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec_p1)]
        //         +          e_lon  * _1_CR           *          e_r  * grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec_p1)];

        // Ti =      _1_CR           * (_1_CR - e_lat) * (_1_CR - e_r) * grid.T_loc[I2V(i_rec,   j_rec,   k_rec)]
        //         + _1_CR           *          e_lat  * (_1_CR - e_r) * grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec)]
        //         + _1_CR           * (_1_CR - e_lat) *          e_r  * grid.T_loc[I2V(i_rec,   j_rec,   k_rec_p1)]
        //         + _1_CR           *          e_lat  *          e_r  * grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec_p1)];

        // Tip =     _1_CR           * (_1_CR - e_lat) * (_1_CR - e_r) * grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec)]
        //         + _1_CR           *          e_lat  * (_1_CR - e_r) * grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec)]
        //         + _1_CR           * (_1_CR - e_lat) *          e_r  * grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec_p1)]
        //         + _1_CR           *          e_lat  *          e_r  * grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec_p1)];

        // DTk = (Tkp - Tk) / delta_r;
        // DTj = (Tjp - Tj) / delta_lat;
        // DTi = (Tip - Ti) / delta_lon;

        CUSTOMREAL DT1, DT2, DT3, DT4, DT5, DT6, DT7, DT8;
        DT1 = (grid.T_loc[I2V(i_rec,     j_rec,     k_rec    + 1)] - grid.T_loc[I2V(i_rec,     j_rec,     k_rec    - 1)]) / _2_CR / delta_r;
        DT2 = (grid.T_loc[I2V(i_rec_p1,  j_rec,     k_rec    + 1)] - grid.T_loc[I2V(i_rec_p1,  j_rec,     k_rec    - 1)]) / _2_CR / delta_r;
        DT3 = (grid.T_loc[I2V(i_rec,     j_rec_p1,  k_rec    + 1)] - grid.T_loc[I2V(i_rec,     j_rec_p1,  k_rec    - 1)]) / _2_CR / delta_r;
        DT4 = (grid.T_loc[I2V(i_rec_p1,  j_rec_p1,  k_rec    + 1)] - grid.T_loc[I2V(i_rec_p1,  j_rec_p1,  k_rec    - 1)]) / _2_CR / delta_r;
        DT5 = (grid.T_loc[I2V(i_rec,     j_rec,     k_rec_p1 + 1)] - grid.T_loc[I2V(i_rec,     j_rec,     k_rec_p1 - 1)]) / _2_CR / delta_r;
        DT6 = (grid.T_loc[I2V(i_rec_p1,  j_rec,     k_rec_p1 + 1)] - grid.T_loc[I2V(i_rec_p1,  j_rec,     k_rec_p1 - 1)]) / _2_CR / delta_r;
        DT7 = (grid.T_loc[I2V(i_rec,     j_rec_p1,  k_rec_p1 + 1)] - grid.T_loc[I2V(i_rec,     j_rec_p1,  k_rec_p1 - 1)]) / _2_CR / delta_r;
        DT8 = (grid.T_loc[I2V(i_rec_p1,  j_rec_p1,  k_rec_p1 + 1)] - grid.T_loc[I2V(i_rec_p1,  j_rec_p1,  k_rec_p1 - 1)]) / _2_CR / delta_r;

        DTk =   (_1_CR - e_lon) * (_1_CR - e_lat) * (_1_CR - e_r) * DT1
              +          e_lon  * (_1_CR - e_lat) * (_1_CR - e_r) * DT2
              + (_1_CR - e_lon) *          e_lat  * (_1_CR - e_r) * DT3
              +          e_lon  *          e_lat  * (_1_CR - e_r) * DT4
              + (_1_CR - e_lon) * (_1_CR - e_lat) *          e_r  * DT5
              +          e_lon  * (_1_CR - e_lat) *          e_r  * DT6
              + (_1_CR - e_lon) *          e_lat  *          e_r  * DT7
              +          e_lon  *          e_lat  *          e_r  * DT8;

        DT1 = (grid.T_loc[I2V(i_rec,     j_rec    + 1,  k_rec   )] - grid.T_loc[I2V(i_rec,     j_rec    - 1,  k_rec   )]) / _2_CR / delta_lat;
        DT2 = (grid.T_loc[I2V(i_rec_p1,  j_rec    + 1,  k_rec   )] - grid.T_loc[I2V(i_rec_p1,  j_rec    - 1,  k_rec   )]) / _2_CR / delta_lat;
        DT3 = (grid.T_loc[I2V(i_rec,     j_rec_p1 + 1,  k_rec   )] - grid.T_loc[I2V(i_rec,     j_rec_p1 - 1,  k_rec   )]) / _2_CR / delta_lat;
        DT4 = (grid.T_loc[I2V(i_rec_p1,  j_rec_p1 + 1,  k_rec   )] - grid.T_loc[I2V(i_rec_p1,  j_rec_p1 - 1,  k_rec   )]) / _2_CR / delta_lat;
        DT5 = (grid.T_loc[I2V(i_rec,     j_rec    + 1,  k_rec_p1)] - grid.T_loc[I2V(i_rec,     j_rec    - 1,  k_rec_p1)]) / _2_CR / delta_lat;
        DT6 = (grid.T_loc[I2V(i_rec_p1,  j_rec    + 1,  k_rec_p1)] - grid.T_loc[I2V(i_rec_p1,  j_rec    - 1,  k_rec_p1)]) / _2_CR / delta_lat;
        DT7 = (grid.T_loc[I2V(i_rec,     j_rec_p1 + 1,  k_rec_p1)] - grid.T_loc[I2V(i_rec,     j_rec_p1 - 1,  k_rec_p1)]) / _2_CR / delta_lat;
        DT8 = (grid.T_loc[I2V(i_rec_p1,  j_rec_p1 + 1,  k_rec_p1)] - grid.T_loc[I2V(i_rec_p1,  j_rec_p1 - 1,  k_rec_p1)]) / _2_CR / delta_lat;

        DTj =   (_1_CR - e_lon) * (_1_CR - e_lat) * (_1_CR - e_r) * DT1
              +          e_lon  * (_1_CR - e_lat) * (_1_CR - e_r) * DT2
              + (_1_CR - e_lon) *          e_lat  * (_1_CR - e_r) * DT3
              +          e_lon  *          e_lat  * (_1_CR - e_r) * DT4
              + (_1_CR - e_lon) * (_1_CR - e_lat) *          e_r  * DT5
              +          e_lon  * (_1_CR - e_lat) *          e_r  * DT6
              + (_1_CR - e_lon) *          e_lat  *          e_r  * DT7
              +          e_lon  *          e_lat  *          e_r  * DT8;

        DT1 = (grid.T_loc[I2V(i_rec    + 1,  j_rec   ,  k_rec   )] - grid.T_loc[I2V(i_rec    - 1,  j_rec   ,  k_rec   )]) / _2_CR / delta_lon;
        DT2 = (grid.T_loc[I2V(i_rec_p1 + 1,  j_rec   ,  k_rec   )] - grid.T_loc[I2V(i_rec_p1 - 1,  j_rec   ,  k_rec   )]) / _2_CR / delta_lon;
        DT3 = (grid.T_loc[I2V(i_rec    + 1,  j_rec_p1,  k_rec   )] - grid.T_loc[I2V(i_rec    - 1,  j_rec_p1,  k_rec   )]) / _2_CR / delta_lon;
        DT4 = (grid.T_loc[I2V(i_rec_p1 + 1,  j_rec_p1,  k_rec   )] - grid.T_loc[I2V(i_rec_p1 - 1,  j_rec_p1,  k_rec   )]) / _2_CR / delta_lon;
        DT5 = (grid.T_loc[I2V(i_rec    + 1,  j_rec   ,  k_rec_p1)] - grid.T_loc[I2V(i_rec    - 1,  j_rec   ,  k_rec_p1)]) / _2_CR / delta_lon;
        DT6 = (grid.T_loc[I2V(i_rec_p1 + 1,  j_rec   ,  k_rec_p1)] - grid.T_loc[I2V(i_rec_p1 - 1,  j_rec   ,  k_rec_p1)]) / _2_CR / delta_lon;
        DT7 = (grid.T_loc[I2V(i_rec    + 1,  j_rec_p1,  k_rec_p1)] - grid.T_loc[I2V(i_rec    - 1,  j_rec_p1,  k_rec_p1)]) / _2_CR / delta_lon;
        DT8 = (grid.T_loc[I2V(i_rec_p1 + 1,  j_rec_p1,  k_rec_p1)] - grid.T_loc[I2V(i_rec_p1 - 1,  j_rec_p1,  k_rec_p1)]) / _2_CR / delta_lon;

        DTi =   (_1_CR - e_lon) * (_1_CR - e_lat) * (_1_CR - e_r) * DT1
              +          e_lon  * (_1_CR - e_lat) * (_1_CR - e_r) * DT2
              + (_1_CR - e_lon) *          e_lat  * (_1_CR - e_r) * DT3
              +          e_lon  *          e_lat  * (_1_CR - e_r) * DT4
              + (_1_CR - e_lon) * (_1_CR - e_lat) *          e_r  * DT5
              +          e_lon  * (_1_CR - e_lat) *          e_r  * DT6
              + (_1_CR - e_lon) *          e_lat  *          e_r  * DT7
              +          e_lon  *          e_lat  *          e_r  * DT8;
        // std::cout << grid.T_loc[I2V(i_rec    + 1,  j_rec   ,  k_rec   )] << ", " << grid.T_loc[I2V(i_rec    - 1,  j_rec   ,  k_rec   )]
        //           << ", " << _2_CR * delta_lon << std::endl;
        // std::cout << DT1 << "," << DT2 << "," << DT3 << "," << DT4 << "," << DT5 << "," << DT6 << "," << DT7 << "," << DT8 <<std::endl;
        // std::cout << "DTi: " << DTi << std::endl;

        // broadcast the gradient components to all processes
        broadcast_cr_single(DTk, rec_rank);
        broadcast_cr_single(DTj, rec_rank);
        broadcast_cr_single(DTi, rec_rank);

    } else {
        // receive the gradient components
        broadcast_cr_single(DTk, rec_rank);
        broadcast_cr_single(DTj, rec_rank);
        broadcast_cr_single(DTi, rec_rank);
    }

    // store the calculated travel time
    rec.DTk = DTk;
    rec.DTj = DTj;
    rec.DTi = DTi;

}


// approximated optimal origin time
void Receiver::calculate_optimal_origin_time(InputParams& IP, const std::string& name_sim_src){
    if (subdom_main) {

        // calculate gradient of travel time at each receiver (swapped source)
        for (auto iter = IP.data_map[name_sim_src].begin(); iter != IP.data_map[name_sim_src].end(); iter++) {
            std::string name_rec = iter->first;
            DataInfo&   data     = iter->second;
            const CUSTOMREAL& weight = data.weight;
            const CUSTOMREAL& misfit = data.travel_time_obs - data.travel_time;
            IP.rec_map[name_rec].tau_opt    += misfit * weight;
            IP.rec_map[name_rec].sum_weight += weight;
        }

    }
}

void Receiver::divide_optimal_origin_time_by_summed_weight(InputParams& IP) {
    if (subdom_main) {

        for (auto iter = IP.rec_map.begin(); iter != IP.rec_map.end();  iter++) {
            allreduce_cr_sim_single_inplace(iter->second.tau_opt);
            allreduce_cr_sim_single_inplace(iter->second.sum_weight);
            iter->second.tau_opt /= iter->second.sum_weight;
            // std::cout << "id_sim" << id_sim << ", name: " << iter->first << ", ortime: " << iter->second.tau_opt <<std::endl;
        }
    }
    //synchronize_all_world(); // not necessary because allreduce is already synchronizing communication
}

void Receiver::calculate_obj_reloc(InputParams& IP, int i_iter){

    for (auto it_src = IP.data_map.begin(); it_src != IP.data_map.end(); it_src++) {
        std::string name_src = it_src->first;
        for (auto it_rec = IP.data_map[name_src].begin(); it_rec != IP.data_map[name_src].end(); it_rec++) {
            std::string name_rec = it_rec->first;
            DataInfo&   data     = it_rec->second;
            const CUSTOMREAL& weight = data.weight;
            const CUSTOMREAL& misfit = data.travel_time - data.travel_time_obs;

            IP.rec_map[name_rec].vobj_src_reloc += weight / _2_CR * my_square(misfit+IP.rec_map[name_rec].tau_opt);
        }
    }

    // sum the obj from all sources (swapped receivers)
    if (subdom_main) {
        for (auto iter = IP.rec_map.begin(); iter != IP.rec_map.end(); iter++){
            allreduce_cr_sim_single_inplace(iter->second.vobj_src_reloc);
        }
    }

    //synchronize_all_world(); // not necessary because allreduce is already synchronizing communication

    for (auto iter = IP.rec_map.begin(); iter != IP.rec_map.end(); iter++){
        CUSTOMREAL obj = iter->second.vobj_src_reloc;
        CUSTOMREAL old_obj = iter->second.vobj_src_reloc_old;
        if (i_iter != 0 && old_obj < obj){    // if obj increase, decrease the step length of this (swapped) source
            // std::cout << "before, step_length_max: " << iter->second.step_length_max << "step_size_decay: " << step_size_decay << std::endl;
            iter->second.step_length_max *= step_length_decay;
            // std::cout << "after, step_length_max: " << iter->second.step_length_max << "step_size_decay: " << step_size_decay << std::endl;
        }
        // std::cout << "id_sim: " << id_sim << ", name: " << iter->first << ", obj: " << obj << ", old obj: " << old_obj << ", step_length_max: " << iter->second.step_length_max
        //           << ", step_size_decay: " << step_size_decay
        //           << std::endl;
    }

}


// calculate the gradient of the objective function
void Receiver::calculate_grad_obj_src_reloc(InputParams& IP, const std::string& name_sim_src) {

    if(subdom_main){
        // get list of receivers from input parameters
        // std::vector<SrcRec>& receivers = IP.get_rec_names(id_sim_src);

        // get the source weight
        // CUSTOMREAL src_weight = IP.get_src_point(id_sim_src).weight;

        // calculate gradient of travel time at each receiver (swapped source)
        for (auto it_rec = IP.data_map[name_sim_src].begin(); it_rec != IP.data_map[name_sim_src].end(); it_rec++){
            std::string name_rec = it_rec->first;
            DataInfo&   data     = it_rec->second;
            const CUSTOMREAL& weight = data.weight;
            const CUSTOMREAL& misfit = data.travel_time - data.travel_time_obs;
            IP.rec_map[name_rec].grad_chi_k += (misfit + IP.rec_map[name_rec].tau_opt) * IP.rec_map[name_rec].DTk * weight;
            IP.rec_map[name_rec].grad_chi_j += (misfit + IP.rec_map[name_rec].tau_opt) * IP.rec_map[name_rec].DTj * weight;
            IP.rec_map[name_rec].grad_chi_i += (misfit + IP.rec_map[name_rec].tau_opt) * IP.rec_map[name_rec].DTi * weight;

            // std::cout << "misfit: " << misfit << "tau_opt: " << IP.rec_map_nv[name_rec].tau_opt << std::endl;
            // std::cout << "DTk, DTj, DTi = " << IP.rec_map_nv[name_rec].DTk << ", " << IP.rec_map_nv[name_rec].DTj
            //           << ", " << IP.rec_map_nv[name_rec].DTi <<std::endl;
        }

        // for(auto iter = IP.rec_map_nv.begin(); iter != IP.rec_map_nv.end(); iter++){
            // std::cout << "grad_chi_i (lon,rad) is: " << iter->second.grad_chi_i << std::endl;
            // std::cout << "grad_chi_j (lat,rad) is: " << iter->second.grad_chi_j << std::endl;
            // std::cout << "grad_chi_k (radius,km) is: " << iter->second.grad_chi_k << std::endl << std::endl;
        // }
    }

    // for(auto iter = IP.rec_map_nv.begin(); iter != IP.rec_map_nv.end(); iter++){
    //     std::cout << "ckp1, id_sim" << id_sim
    //     << ", src name: " << iter->first
    //     << ", kernel k: " << iter->second.grad_chi_k
    //     << ", kernel j: " << iter->second.grad_chi_j
    //     << ", kernel i: " << iter->second.grad_chi_i
    //     << std::endl;
    // }
}

void Receiver::update_source_location(InputParams& IP, Grid& grid) {

    if (subdom_main) {
        // get list of receivers from input parameters
        // std::vector<SrcRec>& receivers = IP.get_rec_names(id_sim_src);

        for(auto iter = IP.rec_map.begin(); iter != IP.rec_map.end(); iter++){
            CUSTOMREAL grad_dep_km = - iter->second.grad_chi_k;
            CUSTOMREAL grad_lat_km = iter->second.grad_chi_j/(R_earth);
            CUSTOMREAL grad_lon_km = iter->second.grad_chi_i/(R_earth * cos(iter->second.lat * DEG2RAD));
            CUSTOMREAL norm_grad   = std::sqrt(my_square(grad_dep_km) + my_square(grad_lat_km) + my_square(grad_lon_km));
            CUSTOMREAL step_length;
            if (iter->second.is_stop || isZero(norm_grad)){
                step_length = 0;
            } else {
                step_length = 0.5 * iter->second.vobj_src_reloc/my_square(norm_grad); // becomes infinity if norm_grad is zero
            }

            // make 0 if step_length is nan or inf
//            if (std::isnan(step_length) || std::isinf(step_length)){
//                step_length = 0;
//            }

//            // check if step_length is nan
//            if (std::isnan(step_length) || std::isinf(step_length)){
//                std::cout << "step_length is nan, id_sim: " << id_sim << ", name: " << iter->first << std::endl;
//                std::cout << "grad_dep_km: " << grad_dep_km << ", grad_lat_km: " << grad_lat_km << ", grad_lon_km: " << grad_lon_km << std::endl;
//                std::cout << "norm_grad: " << norm_grad << ", vobj_src_reloc: " << iter->second.vobj_src_reloc << std::endl;
//                std::cout << "grad_chi_k: " << iter->second.grad_chi_k << ", grad_chi_j: " << iter->second.grad_chi_j << ", grad_chi_i: " << iter->second.grad_chi_i << std::endl;
//                std::cout << "lat: " << iter->second.lat << ", lon: " << iter->second.lon << ", dep: " << iter->second.dep << std::endl;
//                std::cout << "DTk: " << iter->second.DTk << ", DTj: " << iter->second.DTj << ", DTi: " << iter->second.DTi << std::endl;
//                std::cout << "tau_opt: " << iter->second.tau_opt << std::endl;
//                std::cout << "is_stop: " << iter->second.is_stop << std::endl;
//                std::cout << "vobj_src_reloc: " << iter->second.vobj_src_reloc << std::endl;
//                std::cout << "step_length: " << step_length << std::endl;
//            }
//

            // calculated grad are too small or even 0

            //std::cout << "ckp2, id_sim" << id_sim
            //        << ", src name: " << iter->first
            //        << ", kernel k: " << iter->second.grad_chi_k
            //        << ", kernel j: " << iter->second.grad_chi_j
            //        << ", kernel i: " << iter->second.grad_chi_i
            //        << std::endl;

            CUSTOMREAL update_max = -1.0;
            CUSTOMREAL update_dep = step_length * grad_dep_km;
            update_max = std::max(update_max,abs(update_dep));
            CUSTOMREAL update_lat = step_length * grad_lat_km;
            update_max = std::max(update_max,abs(update_lat));
            CUSTOMREAL update_lon = step_length * grad_lon_km;
            update_max = std::max(update_max,abs(update_lon));

            CUSTOMREAL downscale = 1.0;
            if (update_max > iter->second.step_length_max){
                downscale = iter->second.step_length_max / update_max;
            }

            iter->second.vobj_grad_norm_src_reloc = norm_grad;
            iter->second.dep -= update_dep * downscale;
            iter->second.lat -= update_lat * downscale / (R_earth * DEG2RAD);
            iter->second.lon -= update_lon * downscale / (R_earth * DEG2RAD * cos(iter->second.lat * DEG2RAD));

            // detect nan and inf then exit the program
            if (std::isnan(iter->second.dep) || std::isinf(iter->second.dep) ||
                std::isnan(iter->second.lat) || std::isinf(iter->second.lat) ||
                std::isnan(iter->second.lon) || std::isinf(iter->second.lon)){
                std::cout << "Error: nan or inf detected in source relocation!" << std::endl;
                std::cout << "id_sim: " << id_sim
                          << ", src name: " << iter->first
                          << ", obj: " << iter->second.vobj_src_reloc
                          << ", lat: " << iter->second.lat
                          << ", lon: " << iter->second.lon
                          << ", dep: " << iter->second.dep
                          << ", ortime: " << iter->second.tau_opt
                          << ", is_stop: " << iter->second.is_stop
                          << ", grad_dep_km: " << grad_dep_km
                          << ", grad_lat_km: " << grad_lat_km
                          << ", grad_lon_km: " << grad_lon_km
                          << ", vobj_src_reloc: " << iter->second.vobj_src_reloc
                          << std::endl;

                exit(1);
            }


            if (norm_grad < TOL_SRC_RELOC || iter->second.step_length_max < TOL_STEP_SIZE){
                iter->second.is_stop = true;
            }

            // if(id_sim == 0 && myrank == 0 && iter->first == "s1"){
            //     std::cout << "id_sim: " << id_sim
            //           << ", src name: " << iter->first
            //           << ", obj: " << iter->second.vobj_src_reloc
            //           << ", lat: " << iter->second.lat
            //           << ", lon: " << iter->second.lon
            //           << ", dep: " << iter->second.dep
            //           << ", ortime: " << iter->second.tau_opt
            //           << ", is_stop: " << iter->second.is_stop
            //           << std::endl;
            //     std::cout << "time, r1-s1: " << IP.syn_time_map_sr["r1"]["s1"]
            //               << ", src_lat: " << IP.rec_map_nv["s1"].lat
            //               << ", src_lon: " << IP.rec_map_nv["s1"].lon
            //               << ", src_dep: " << IP.rec_map_nv["s1"].dep
            //               << ", rec_lat: " << IP.src_map_nv["r1"].lat
            //               << ", rec_lon: " << IP.src_map_nv["r1"].lon
            //               << ", rec_dep: " << IP.src_map_nv["r1"].dep
            //               << std::endl;
            // }

            // std::cout << "update_dep: " << update_dep
            //           << ", update_lat: " << update_lat
            //           << ", update_lon: " << update_lon
            //           << ", update_max: " << update_max
            //           << ", step_length_max: " << iter->second.step_length_max
            //           << ", downscale: " << downscale
            //           << ", is_stop: "   << iter->second.is_stop
            //           << std::endl;


            // std::cout << "ckp3,id_sim" << id_sim
            //           << ", src name: " << iter->first
            //           << ", kernel k(dep,km): " << grad_dep_km
            //           << ", kernel j(lon,km): " << grad_lat_km
            //           << ", kernel i(lat,km): " << grad_lon_km
            //           << std::endl;


            // std::cout << "dep update km: " << update_dep * downscale
            //           << ", lat update km: " << update_lat * downscale
            //           << ", lon update km: " << update_lon * downscale
            //           << std::endl;

            // check if the new receiver position is within the domain
            // if not then set the receiver position to the closest point on the domain

            // grid size + 1% mergin to avoid the receiver position is exactly on the boundary
            CUSTOMREAL mergin_lon = 1.01 * grid.get_delta_lon();
            CUSTOMREAL mergin_lat = 1.01 * grid.get_delta_lat();
            CUSTOMREAL mergin_r   = 1.01 * grid.get_delta_r();

            if (iter->second.lon < IP.get_min_lon()*RAD2DEG)
                iter->second.lon = IP.get_min_lon()*RAD2DEG + mergin_lon;
            if (iter->second.lon > IP.get_max_lon()*RAD2DEG)
                iter->second.lon = IP.get_max_lon()*RAD2DEG - mergin_lon;
            if (iter->second.lat < IP.get_min_lat()*RAD2DEG)
                iter->second.lat = IP.get_min_lat()*RAD2DEG + mergin_lat;
            if (iter->second.lat > IP.get_max_lat()*RAD2DEG)
                iter->second.lat = IP.get_max_lat()*RAD2DEG - mergin_lat;
            if (iter->second.dep < IP.get_min_dep())
                iter->second.dep = IP.get_min_dep() + mergin_r;
            if (iter->second.dep > IP.get_max_dep())
                iter->second.dep = IP.get_max_dep() - mergin_r;
        }
    }
}
