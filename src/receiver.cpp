#include "receiver.h"


Receiver::Receiver() {
}


Receiver::~Receiver() {
}


void Receiver::calculate_arrival_time(InputParams& IP, Grid& grid) {
    if(subdom_main){
        // get receiver positions from input parameters
        std::vector<SrcRec>& receivers = IP.get_rec_points(id_sim_src);

        // calculate the travel time of the receiver by interpolation
        for (auto& rec: receivers) {
            if(!rec.is_rec_pair)
                interpolate_travel_time(grid, rec);      // calculate the travel time of the receiver by interpolation
            else
                interpolate_differential_travel_time(grid, rec);   // calculate the differential travel time of the receiver pair by interpolation
        }
    }
}

std::vector<CUSTOMREAL> Receiver::calculate_adjoint_source(InputParams& IP) {

    CUSTOMREAL allsum_obj = 0;
    CUSTOMREAL allsum_misfit = 0;

    std::vector<CUSTOMREAL> allsum = std::vector<CUSTOMREAL>(2);

    if(subdom_main){
        // get receiver positions from input parameters
        std::vector<SrcRec>& receivers = IP.get_rec_points(id_sim_src);

        // get the source weight
        CUSTOMREAL src_weight = IP.get_src_point(id_sim_src).weight;

        // store sum of adjoint sources for the current source as an objective function value
        CUSTOMREAL sum_adj_src = 0;

        // calculate the adjoint source of the receiver by interpolation
        for (auto& rec: receivers) {
            if(!rec.is_rec_pair){
                rec.t_adj = (rec.arr_time - rec.arr_time_ori) * (rec.weight * src_weight);
                sum_adj_src += my_square(rec.t_adj); // multiply by weight and
                allsum_misfit += my_square(rec.t_adj);
            } else {
                if(!IP.get_src_point(id_sim_src).is_teleseismic){
                    // differential traveltimes of local earthquake, do not consider station correction
                    rec.ddt_adj_pair[0] = (rec.dif_arr_time - rec.dif_arr_time_ori) * (rec.weight * src_weight);
                    rec.ddt_adj_pair[1] = -rec.ddt_adj_pair[0];
                    sum_adj_src += my_square(rec.ddt_adj_pair[0]);
                    allsum_misfit += my_square(rec.ddt_adj_pair[0]);
                } else {
                    // differential traveltimes of teleseismic earthquake, consider station correction
                    rec.ddt_adj_pair[0] = (rec.dif_arr_time + (rec.station_correction_pair[0] - rec.station_correction_pair[1]) - rec.dif_arr_time_ori) * (rec.weight * src_weight);
                    rec.ddt_adj_pair[1] = -rec.ddt_adj_pair[0];
                    sum_adj_src += my_square(rec.ddt_adj_pair[0]);
                    allsum_misfit += my_square(rec.ddt_adj_pair[0]);
                    // std::cout << "dif_arr_time: " << rec.dif_arr_time << ", correction_pair[0]: " << rec.station_correction_pair[0]
                    //           << ", correction_pair[1]: " << rec.station_correction_pair[1]
                    //           << ", dif_arr_time_ori: " << rec.dif_arr_time_ori << std::endl;
                }

            }
        }

        allsum_obj = sum_adj_src;
    }

    // share the calculated objective function value to other processes
    broadcast_cr_single_sub(allsum_obj,0);
    broadcast_cr_single_sub(allsum_misfit,0);

    allsum[0] = allsum_obj;
    allsum[1] = allsum_misfit;

    return allsum;
}


void Receiver::interpolate_travel_time(Grid& grid, SrcRec& rec) {
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
        std::cout << "id_src: " << rec.id_src << " id_rec: " << rec.id_rec << " n_rec: " << rec.n_rec << " depth: " << rec.dep << " lat: " << rec.lat << " lon: " << rec.lon << " arr_time: " << rec.arr_time << std::endl;
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
            std::cout << "id_src: " << rec.id_src << " id_rec: " << rec.id_rec << " n_rec: " << rec.n_rec << " depth: " << rec.dep << " lat: " << rec.lat << " lon: " << rec.lon << " arr_time: " << rec.arr_time << std::endl;
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

        // broadcast interpolated travel time
        broadcast_cr_single(vinterp, rec_rank);

    } else {
        // receive the calculated traveltime
        broadcast_cr_single(vinterp, rec_rank);
    }

    // store the calculated travel time
    rec.arr_time = vinterp;

}


void Receiver::interpolate_differential_travel_time(Grid& grid, SrcRec& rec) {
    // calculate the differential travel time of the receiver by 3d linear interpolation

    // copy some parameters
    CUSTOMREAL delta_lon = grid.get_delta_lon();
    CUSTOMREAL delta_lat = grid.get_delta_lat();
    CUSTOMREAL delta_r   = grid.get_delta_r();

    // store receiver position in radian

    std::vector<CUSTOMREAL> vinterp(2);

    for(int index_rec = 0; index_rec < 2; index_rec ++){
        // store receiver position in radian
        CUSTOMREAL rec_lon = rec.lon_pair[index_rec]*DEG2RAD;
        CUSTOMREAL rec_lat = rec.lat_pair[index_rec]*DEG2RAD;
        CUSTOMREAL rec_r   = depth2radius(rec.dep_pair[index_rec]);
        CUSTOMREAL rec_id  = rec.id_rec_pair[index_rec];

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
            std::cout << "id_src: " << rec.id_src << " id_rec: " << rec_id << " n_rec_pair: " << rec.n_rec_pair << " depth: " << radius2depth(rec_r) \
                      << " lat: " << rec_lat << " lon: " << rec_lon << std::endl;
            // print boundary
            //std::cout << "lon min max rec: " << grid.get_lon_min_loc() << " " << grid.get_lon_max_loc() << " " << rec_lon << std::endl;
            //std::cout << "lat min max rec: " << grid.get_lat_min_loc() << " " << grid.get_lat_max_loc() << " " << rec_lat << std::endl;
            //std::cout << "r min max rec: " << grid.get_r_min_loc() << " " << grid.get_r_max_loc() << " " << rec_r << std::endl;

            std::cout << "lon+bound min max rec: " << (grid.get_lon_min_loc() - delta_lon)*RAD2DEG     << " " << (grid.get_lon_max_loc() + delta_lon)*RAD2DEG     << " " << rec_lon*RAD2DEG    << std::endl;
            std::cout << "lat+bound min max rec: " << (grid.get_lat_min_loc() - delta_lat)*RAD2DEG     << " " << (grid.get_lat_max_loc() + delta_lat)*RAD2DEG     << " " << rec_lat*RAD2DEG    << std::endl;
            std::cout << "r+bound min max rec: "   << radius2depth(grid.get_r_min_loc()   - delta_r  ) << " " << radius2depth(grid.get_r_max_loc()   + delta_r  ) << " " << radius2depth(rec_r)<< std::endl;
            exit(1);
        }

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

            if(if_verbose){
                std::cout << "(rec_lon - dis_rec_lon)/dlon: " << (rec_lon - dis_rec_lon)/delta_lon << std::endl;
                std::cout << "(rec_lat - dis_rec_lat)/dlat: " << (rec_lat - dis_rec_lat)/delta_lat << std::endl;
                std::cout << "(rec_r   - dis_rec_r  )/dr: "   << (rec_r   - dis_rec_r  )/delta_r << std::endl;
                std::cout << "rec_lon, dis_rec_lon, delta_lon: " << rec_lon << " " << dis_rec_lon << " " << delta_lon << std::endl;
                std::cout << "rec_lat, dis_rec_lat, delta_lat: " << rec_lat << " " << dis_rec_lat << " " << delta_lat << std::endl;
                std::cout << "rec_r, dis_rec_r, delta_r  : " << rec_r << ", " << dis_rec_r << ", " <<  delta_r  << std::endl;
                std::cout << "loc_K. k_rec: " << loc_K << "," << k_rec << std::endl;
                std::cout << "r_loc_1d[loc_K-1]: " << grid.r_loc_1d[loc_K-1] << std::endl;
            }

            int i_rec_p1 = i_rec + 1;
            int j_rec_p1 = j_rec + 1;
            int k_rec_p1 = k_rec + 1;

            // exclude the points if they are out of the domain
            if (i_rec_p1 > loc_I-1 \
            || j_rec_p1 > loc_J-1 \
            || k_rec_p1 > loc_K-1) {
                // exit(1) as the source is out of the domain
                std::cout << "Error: the receiver is out of the domain" << std::endl;
                std::cout << "id_src: " << rec.id_src << " id_rec: " << rec_id << " n_rec_pair: " << rec.n_rec_pair << " depth: " << radius2depth(rec_r) \
                          << " lat: " << rec_lat << " lon: " << rec_lon << std::endl;
                std::cout << "lon min max rec: " << grid.get_lon_min_loc()*RAD2DEG << " " << grid.get_lon_max_loc()*RAD2DEG << " " << rec_lon*RAD2DEG << std::endl;
                std::cout << "lat min max rec: " << grid.get_lat_min_loc()*RAD2DEG << " " << grid.get_lat_max_loc()*RAD2DEG << " " << rec_lat*RAD2DEG << std::endl;
                std::cout << "r min max rec: " << radius2depth(grid.get_r_min_loc()) << " " << radius2depth(grid.get_r_max_loc()) << " " << radius2depth(rec_r) << std::endl;
                std::cout << "i_rec: " << i_rec << " j_rec: " << j_rec << " k_rec: " << k_rec << std::endl;
                std::cout << "i_rec_p1: " << i_rec_p1 << " j_rec_p1: " << j_rec_p1 << " k_rec_p1: " << k_rec_p1 << std::endl;
                std::cout << "loc_I-1: " << loc_I-1 << " loc_J-1: " << loc_J-1 << " loc_K-1: " << loc_K-1 << std::endl;
                //finalize_mpi();
                exit(1);
            }


            vinterp[index_rec] =  (_1_CR - e_lon) * (_1_CR - e_lat) * (_1_CR - e_r) * grid.T_loc[I2V(i_rec,   j_rec,   k_rec)]   \
                                +          e_lon  * (_1_CR - e_lat) * (_1_CR - e_r) * grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec)]   \
                                + (_1_CR - e_lon) *          e_lat  * (_1_CR - e_r) * grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec)]   \
                                + (_1_CR - e_lon) * (_1_CR - e_lat) *          e_r  * grid.T_loc[I2V(i_rec,   j_rec,   k_rec_p1)] \
                                +          e_lon  *          e_lat  * (_1_CR - e_r) * grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec)]   \
                                +          e_lon  * (_1_CR - e_lat) *          e_r  * grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec_p1)] \
                                + (_1_CR - e_lon) *          e_lat  *          e_r  * grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec_p1)] \
                                +          e_lon  *          e_lat  *          e_r  * grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec_p1)];

            //std::cout << "DEBUG near and vinterp : " << grid.T_loc[I2V(i_rec,j_rec,k_rec)] << ", " << vinterp << std::endl;

            // broadcast interpolated travel time
            broadcast_cr_single(vinterp[index_rec], rec_rank);

        } else {
            // receive the calculated traveltime
            broadcast_cr_single(vinterp[index_rec], rec_rank);
        }

    }
    // store the calculated travel time
    rec.dif_arr_time = vinterp[0] - vinterp[1];
}


void Receiver::init_vars_src_reloc(InputParams& IP,std::vector<SrcRec>& unique_rec_list){
    if (subdom_main) {
        // get list of receivers from input parameters
        std::vector<SrcRec>& receivers = IP.get_rec_points(id_sim_src);

        // calculate gradient of travel time at each receiver (swapped source)
        for (auto& rec: receivers) {
            rec.tau_opt    = _0_CR;
            rec.grad_chi_i = _0_CR;
            rec.grad_chi_j = _0_CR;
            rec.grad_chi_k = _0_CR;
            rec.sum_weight = _0_CR;

            // initialize on unique list
            int target_id = rec.id_unique_list;
            unique_rec_list[target_id].tau_opt                  = _0_CR;
            unique_rec_list[target_id].grad_chi_i               = _0_CR;
            unique_rec_list[target_id].grad_chi_j               = _0_CR;
            unique_rec_list[target_id].grad_chi_k               = _0_CR;
            unique_rec_list[target_id].sum_weight               = _0_CR;
            unique_rec_list[target_id].vobj_src_reloc           = _0_CR;
            unique_rec_list[target_id].vobj_grad_norm_src_reloc = _0_CR;
        }

    }
}


void Receiver::calculate_T_gradient(InputParams& IP, Grid& grid){

    if(subdom_main){
        // get list of receivers from input parameters
        std::vector<SrcRec>& receivers = IP.get_rec_points(id_sim_src);

        // calculate gradient of travel time at each receiver (swapped source)
        for (auto& rec: receivers) {
            calculate_T_gradient_one_rec(grid, rec);
        }

    }

}


void Receiver::calculate_T_gradient_one_rec(Grid& grid, SrcRec& rec){

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
        std::cout << "id_src: " << rec.id_src << " id_rec: " << rec.id_rec << " n_rec: " << rec.n_rec << " depth: " << rec.dep << " lat: " << rec.lat << " lon: " << rec.lon << " arr_time: " << rec.arr_time << std::endl;
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

        // exclude the points if they are out of the domain
        if (i_rec_p1 > loc_I-1 \
         || j_rec_p1 > loc_J-1 \
         || k_rec_p1 > loc_K-1) {
            // exit(1) as the source is out of the domain
            std::cout << "Error: the receiver is out of the domain" << std::endl;
            std::cout << "id_src: " << rec.id_src << " id_rec: " << rec.id_rec << " n_rec: " << rec.n_rec << " depth: " << rec.dep << " lat: " << rec.lat << " lon: " << rec.lon << " arr_time: " << rec.arr_time << std::endl;
            std::cout << "lon min max rec: " << grid.get_lon_min_loc()*RAD2DEG << " " << grid.get_lon_max_loc()*RAD2DEG << " " << rec_lon*RAD2DEG << std::endl;
            std::cout << "lat min max rec: " << grid.get_lat_min_loc()*RAD2DEG << " " << grid.get_lat_max_loc()*RAD2DEG << " " << rec_lat*RAD2DEG << std::endl;
            std::cout << "r min max rec: " << radius2depth(grid.get_r_min_loc()) << " " << radius2depth(grid.get_r_max_loc()) << " " << radius2depth(rec_r) << std::endl;
            std::cout << "i_rec: " << i_rec << " j_rec: " << j_rec << " k_rec: " << k_rec << std::endl;
            std::cout << "i_rec_p1: " << i_rec_p1 << " j_rec_p1: " << j_rec_p1 << " k_rec_p1: " << k_rec_p1 << std::endl;
            std::cout << "loc_I-1: " << loc_I-1 << " loc_J-1: " << loc_J-1 << " loc_K-1: " << loc_K-1 << std::endl;
            //finalize_mpi();
            exit(1);
        }

        CUSTOMREAL Ti, Tip, Tj, Tjp, Tk, Tkp;
        Tk =      (_1_CR - e_lon) * (_1_CR - e_lat) * _1_CR         * grid.T_loc[I2V(i_rec,   j_rec,   k_rec)]   \
                +          e_lon  * (_1_CR - e_lat) * _1_CR         * grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec)]   \
                + (_1_CR - e_lon) *          e_lat  * _1_CR         * grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec)]   \
                +          e_lon  *          e_lat  * _1_CR         * grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec)];

        Tkp =     (_1_CR - e_lon) * (_1_CR - e_lat) * _1_CR         * grid.T_loc[I2V(i_rec,   j_rec,   k_rec_p1)] \
                +          e_lon  * (_1_CR - e_lat) * _1_CR         * grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec_p1)] \
                + (_1_CR - e_lon) *          e_lat  * _1_CR         * grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec_p1)] \
                +          e_lon  *          e_lat  * _1_CR         * grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec_p1)];

        Tj =      (_1_CR - e_lon) * _1_CR           * (_1_CR - e_r) * grid.T_loc[I2V(i_rec,   j_rec,   k_rec)]    \
                +          e_lon  * _1_CR           * (_1_CR - e_r) * grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec)]    \
                + (_1_CR - e_lon) * _1_CR           *          e_r  * grid.T_loc[I2V(i_rec,   j_rec,   k_rec_p1)] \
                +          e_lon  * _1_CR           *          e_r  * grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec_p1)];

        Tjp =     (_1_CR - e_lon) * _1_CR           * (_1_CR - e_r) * grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec)]    \
                +          e_lon  * _1_CR           * (_1_CR - e_r) * grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec)]    \
                + (_1_CR - e_lon) * _1_CR           *          e_r  * grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec_p1)] \
                +          e_lon  * _1_CR           *          e_r  * grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec_p1)];

        Ti =      _1_CR           * (_1_CR - e_lat) * (_1_CR - e_r) * grid.T_loc[I2V(i_rec,   j_rec,   k_rec)]    \
                + _1_CR           *          e_lat  * (_1_CR - e_r) * grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec)]    \
                + _1_CR           * (_1_CR - e_lat) *          e_r  * grid.T_loc[I2V(i_rec,   j_rec,   k_rec_p1)] \
                + _1_CR           *          e_lat  *          e_r  * grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec_p1)];

        Tip =     _1_CR           * (_1_CR - e_lat) * (_1_CR - e_r) * grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec)]    \
                + _1_CR           *          e_lat  * (_1_CR - e_r) * grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec)]    \
                + _1_CR           * (_1_CR - e_lat) *          e_r  * grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec_p1)] \
                + _1_CR           *          e_lat  *          e_r  * grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec_p1)];

        DTk = (Tkp - Tk) / delta_r;
        DTj = (Tjp - Tj) / delta_lat;
        DTi = (Tip - Ti) / delta_lon;

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
void Receiver::calculate_optimal_origin_time(InputParams& IP, std::vector<SrcRec>& unique_rec_list){
    if (subdom_main) {
        // get list of receivers from input parameters
        std::vector<SrcRec>& receivers = IP.get_rec_points(id_sim_src);

        // get the source weight
        CUSTOMREAL src_weight = IP.get_src_point(id_sim_src).weight;

        // calculate gradient of travel time at each receiver (swapped source)
        for (auto& rec: receivers) {
            SrcRec& rec_unique = unique_rec_list[rec.id_unique_list];
            rec_unique.tau_opt    += (rec.arr_time_ori - rec.arr_time) * rec.weight;
            rec_unique.sum_weight += rec.weight * src_weight;

            // calculate objective function value
            rec_unique.vobj_src_reloc += rec.weight / _2_CR * my_square(rec.arr_time - rec.arr_time_ori);

            // calculate grad norm of objective function value
            rec_unique.vobj_grad_norm_src_reloc += std::sqrt(rec.DTk*rec.DTk + rec.DTj*rec.DTj + rec.DTi*rec.DTi);
        }
    }
}


void Receiver::divide_optimal_origin_time_by_summed_weight(InputParams& IP, std::vector<SrcRec>& unique_rec_list) {
    if (subdom_main) {
        for (auto& rec_unique: unique_rec_list) {
            rec_unique.tau_opt /= rec_unique.sum_weight;
        }
    }
}


// calculate the gradient of the objective function
void Receiver::calculate_grad_obj_src_reloc(InputParams& IP, std::vector<SrcRec>& unique_rec_list) {

    if(subdom_main){
        // get list of receivers from input parameters
        std::vector<SrcRec>& receivers = IP.get_rec_points(id_sim_src);

        // get the source weight
        CUSTOMREAL src_weight = IP.get_src_point(id_sim_src).weight;

        // calculate gradient of travel time at each receiver (swapped source)
        for (auto& rec: receivers) {
            SrcRec& rec_unique = unique_rec_list[rec.id_unique_list];

            rec_unique.grad_chi_k += (rec.arr_time - rec.arr_time_ori + rec_unique.tau_opt) * rec.DTk * rec.weight * src_weight;
            rec_unique.grad_chi_j += (rec.arr_time - rec.arr_time_ori + rec_unique.tau_opt) * rec.DTj * rec.weight * src_weight;
            rec_unique.grad_chi_i += (rec.arr_time - rec.arr_time_ori + rec_unique.tau_opt) * rec.DTi * rec.weight * src_weight;

            //// debug
            //if (myrank == 0){
            //    std::cout << "DTijk " << rec.DTi << " " << rec.DTj << " " << rec.DTk << std::endl;
            //    std::cout << "tau_opt, grad_chi_i, grad_chi_j, grad_chi_k: " << rec_unique.tau_opt << ", " << rec_unique.grad_chi_i << ", " << rec_unique.grad_chi_j << ", " << rec_unique.grad_chi_k << std::endl;

            //}
        }
    }

}


void Receiver::update_source_location(InputParams& IP, Grid& grid, std::vector<SrcRec>& unique_rec_list) {

    if (subdom_main) {
        // get list of receivers from input parameters
        std::vector<SrcRec>& receivers = IP.get_rec_points(id_sim_src);

        for (auto& rec: receivers) {
            SrcRec& rec_unique = unique_rec_list[rec.id_unique_list];
            rec.dep -= rec_unique.grad_chi_k * step_length_src_reloc;
            rec.lat -= rec_unique.grad_chi_j * step_length_src_reloc;
            rec.lon -= rec_unique.grad_chi_i * step_length_src_reloc;

            // check if the new receiver position is within the domain
            // if not then set the receiver position to the closest point on the domain

            CUSTOMREAL mergin_lon = 1.01 * grid.get_delta_lon();
            CUSTOMREAL mergin_lat = 1.01 * grid.get_delta_lat();
            CUSTOMREAL mergin_r   = 1.01 * grid.get_delta_r();

            if (rec.lon < IP.get_min_lon()*RAD2DEG)
                rec.lon = IP.get_min_lon()*RAD2DEG + mergin_lon;
            if (rec.lon > IP.get_max_lon()*RAD2DEG)
                rec.lon = IP.get_max_lon()*RAD2DEG - mergin_lon;
            if (rec.lat < IP.get_min_lat()*RAD2DEG)
                rec.lat = IP.get_min_lat()*RAD2DEG + mergin_lat;
            if (rec.lat > IP.get_max_lat()*RAD2DEG)
                rec.lat = IP.get_max_lat()*RAD2DEG - mergin_lat;
            if (rec.dep < IP.get_min_dep())
                rec.dep = IP.get_min_dep() + mergin_r;
            if (rec.dep > IP.get_max_dep())
                rec.dep = IP.get_max_dep() - mergin_r;
       }
    }
}

