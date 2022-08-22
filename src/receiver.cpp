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
            // calculate the travel time of the receiver by interpolation
            interpolate_travel_time(grid, rec);
        }
    }
}


CUSTOMREAL Receiver::calculate_adjoint_source(InputParams& IP) {

    CUSTOMREAL allsum_obj = 0;

    if(subdom_main){
        // get receiver positions from input parameters
        std::vector<SrcRec>& receivers = IP.get_rec_points(id_sim_src);

        // store sum of adjoint sources for the current source as an objective function value
        CUSTOMREAL sum_adj_src = 0;

        // calculate the adjoint source of the receiver by interpolation
        for (auto& rec: receivers) {
            rec.t_adj = rec.arr_time - rec.arr_time_ori;
            sum_adj_src += std::pow(rec.t_adj,2);
        }

        allsum_obj = sum_adj_src;
    }

    // share the calculated objective function value to other processes
    broadcast_cr_single_sub(allsum_obj,0);

    return allsum_obj;
}


CUSTOMREAL Receiver::calculate_adjoint_source_teleseismic(InputParams& IP) {

    CUSTOMREAL allsum_obj = 0;

    if(subdom_main){
        // get receiver positions from input parameters
        std::vector<SrcRec>& receivers = IP.get_rec_points(id_sim_src);

        // store sum of adjoint sources for the current source as an objective function value
        CUSTOMREAL sum_adj_src = 0;

        // initialize adjoint sources
        for (auto& rec: receivers) {
            rec.t_adj = 0;
        }

        // calculate double-difference travel times
        int nrec = receivers.size();
        for (int irec=0; irec<nrec-1; irec++) {
            SrcRec& rec_i = receivers[irec];
            for (int jrec=irec+1; jrec<nrec; jrec++) {

                SrcRec& rec_j = receivers[jrec];

                CUSTOMREAL deg_diff;
                Epicentral_distance_sphere(rec_i.lat*DEG2RAD, rec_i.lon*DEG2RAD, rec_j.lat*DEG2RAD, rec_j.lon*DEG2RAD, deg_diff);

                if (std::abs(deg_diff) < DIST_SRC_DDT){
                    CUSTOMREAL DDT = (rec_i.arr_time - rec_j.arr_time) - (rec_i.arr_time_ori - rec_j.arr_time_ori);
                    rec_i.t_adj += DDT;
                    rec_j.t_adj -= DDT;
                    sum_adj_src += std::pow(DDT,2);
                }
            }
        }

        allsum_obj = sum_adj_src;
    }

    // share the calculated objective function value to other processes
    broadcast_cr_single_sub(allsum_obj,0);

    return allsum_obj;
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
    if (grid.get_lon_min_loc() <= rec_lon && rec_lon <= grid.get_lon_max_loc() && \
        grid.get_lat_min_loc() <= rec_lat && rec_lat <= grid.get_lat_max_loc() && \
        grid.get_r_min_loc()   <= rec_r   && rec_r   <= grid.get_r_max_loc()   ) {
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


