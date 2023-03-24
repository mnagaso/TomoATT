#include "receiver.h"


Receiver::Receiver() {
}


Receiver::~Receiver() {
}


// void Receiver::calculate_arrival_time(InputParams& IP, Grid& grid) {
//     if(subdom_main){
//         // get receiver positions from input parameters
//         std::vector<SrcRec>& receivers = IP.get_rec_points(id_sim_src);

//         // calculate the travel time of the receiver by interpolation
//         for (auto& rec: receivers) {
//             if(!rec.is_rec_pair)
//                 interpolate_travel_time(grid, rec);      // calculate the travel time of the receiver by interpolation
//             else
//                 interpolate_differential_travel_time(grid, rec);   // calculate the differential travel time of the receiver pair by interpolation
//         }
//     }
// }

void Receiver::store_arrival_time(InputParams& IP, Grid& grid, std::string name_src) {
    if(subdom_main){
        // get receiver positions from input parameters
        std::vector<std::string> name_receivers = IP.get_rec_points_nv(name_sim_src);

        // calculate the travel time of the receiver by interpolation
        for(std::string name_rec: name_receivers) {
            interpolate_travel_time_nv(grid, IP, name_src, name_rec);
        }
    }
}

// std::vector<CUSTOMREAL> Receiver::calculate_adjoint_source(InputParams& IP) {

//     CUSTOMREAL allsum_obj = 0;
//     CUSTOMREAL allsum_misfit = 0;

//     std::vector<CUSTOMREAL> allsum = std::vector<CUSTOMREAL>(2);

//     if(subdom_main){
//         // get receiver positions from input parameters
//         std::vector<SrcRec>& receivers = IP.get_rec_points(id_sim_src);

//         // get the source weight
//         CUSTOMREAL src_weight = IP.get_src_point(id_sim_src).weight;

//         // store sum of adjoint sources for the current source as an objective function value
//         CUSTOMREAL sum_adj_src = 0;

//         // calculate the adjoint source of the receiver by interpolation
//         for (auto& rec: receivers) {
//             if(!rec.is_rec_pair){
//                 rec.t_adj = rec.arr_time - rec.arr_time_ori;
//                 sum_adj_src += my_square(rec.t_adj) * rec.weight * src_weight; // multiply by weight and
//                 allsum_misfit += my_square(rec.t_adj);
//             } else {
//                 if(!IP.get_src_point(id_sim_src).is_teleseismic){        // differential traveltimes of local earthquake, do not consider station correction
//                     rec.ddt_adj_pair[0] = rec.dif_arr_time - rec.dif_arr_time_ori;
//                     rec.ddt_adj_pair[1] = -rec.ddt_adj_pair[0];
//                     sum_adj_src += my_square(rec.ddt_adj_pair[0]) * rec.weight * src_weight;
//                     allsum_misfit += my_square(rec.ddt_adj_pair[0]);
//                 } else {                        // differential traveltimes of teleseismic earthquake, consider station correction
//                     rec.ddt_adj_pair[0] = rec.dif_arr_time + (rec.station_correction_pair[0] - rec.station_correction_pair[1]) - rec.dif_arr_time_ori;
//                     rec.ddt_adj_pair[1] = -rec.ddt_adj_pair[0];
//                     sum_adj_src += my_square(rec.ddt_adj_pair[0]) * rec.weight * src_weight;
//                     allsum_misfit += my_square(rec.ddt_adj_pair[0]);
//                     // std::cout << "dif_arr_time: " << rec.dif_arr_time << ", correction_pair[0]: " << rec.station_correction_pair[0] 
//                     //           << ", correction_pair[1]: " << rec.station_correction_pair[1] 
//                     //           << ", dif_arr_time_ori: " << rec.dif_arr_time_ori << std::endl; 
//                 }

//             }
//         }

//         allsum_obj = sum_adj_src;
//     }

//     // share the calculated objective function value to other processes
//     broadcast_cr_single_sub(allsum_obj,0);
//     broadcast_cr_single_sub(allsum_misfit,0);

//     allsum[0] = allsum_obj;
//     allsum[1] = allsum_misfit;

//     return allsum;
// }

std::vector<CUSTOMREAL> Receiver::calculate_adjoint_source_nv(InputParams& IP) {

    CUSTOMREAL obj = 0;
    CUSTOMREAL obj_abs = 0;
    CUSTOMREAL obj_cs_dif = 0;
    CUSTOMREAL obj_cr_dif = 0;
    CUSTOMREAL obj_tele = 0;
    CUSTOMREAL misfit = 0;
    CUSTOMREAL misfit_abs = 0;
    CUSTOMREAL misfit_cs_dif = 0;
    CUSTOMREAL misfit_cr_dif = 0;
    CUSTOMREAL misfit_tele = 0;
     
    // obj, obj_abs, obj_cs_dif, obj_cr_dif, obj_tele, misfit, misfit_abs, misfit_cs_dif, misfit_cr_dif, misfit_tele
    std::vector<CUSTOMREAL> allsum = std::vector<CUSTOMREAL>(20);

    // rec.adjoint_source = 0
    IP.initialize_adjoint_source();

    int count = 0;
    
    if(subdom_main){
        for (auto data : IP.data_info_smap[name_sim_src]){     // loop all data related to this source
            if (data.is_src_rec){   // absolute traveltime
                std::string name_src = data.name_src;
                std::string name_rec  = data.name_rec;
                CUSTOMREAL syn_time = IP.syn_time_list_sr[name_src][data.name_rec];
                CUSTOMREAL obs_time = data.travel_time_obs;
                CUSTOMREAL adjoint_source = IP.get_rec_list_nv(name_rec).adjoint_source + (syn_time - obs_time)*data.weight;
                IP.set_adjoint_source(name_rec, adjoint_source);

                // contribute misfit
                obj     += 1.0 * my_square(syn_time - obs_time)*data.weight;
                misfit  += 1.0 * my_square(syn_time - obs_time);
                if (IP.get_src_list_nv(name_src).is_out_of_region || IP.get_rec_list_nv(name_rec).is_out_of_region){
                    obj_tele        +=  1.0 * my_square(syn_time - obs_time)*data.weight;
                    misfit_tele     +=  1.0 * my_square(syn_time - obs_time);
                } else{
                    obj_abs         +=  1.0 * my_square(syn_time - obs_time)*data.weight;
                    misfit_abs      +=  1.0 * my_square(syn_time - obs_time);
                }
                
            } else if (data.is_src_pair) {  // common receiver differential traveltime 
                std::string name_src1 = data.name_src_pair[0];
                std::string name_src2 = data.name_src_pair[1];
                std::string name_rec  = data.name_rec_single;
                if (name_sim_src == name_src1){
                    CUSTOMREAL syn_dif_time = IP.syn_time_list_sr[name_src1][name_rec] - IP.syn_time_list_sr[name_src2][name_rec];
                    CUSTOMREAL obs_dif_time = data.cr_dif_travel_time_obs;
                    CUSTOMREAL adjoint_source = IP.get_rec_list_nv(name_rec).adjoint_source + (syn_dif_time - obs_dif_time)*data.weight;
                    IP.set_adjoint_source(name_rec, adjoint_source);
                    
                    // contribute misfit
                    obj     += 0.5 * my_square(syn_dif_time - obs_dif_time)*data.weight;
                    misfit  += 0.5 * my_square(syn_dif_time - obs_dif_time);
                    if (IP.get_src_list_nv(name_src1).is_out_of_region || \
                        IP.get_src_list_nv(name_src2).is_out_of_region || \
                        IP.get_rec_list_nv(name_rec).is_out_of_region){
                        obj_tele        += 0.5 * my_square(syn_dif_time - obs_dif_time)*data.weight;   
                        misfit_tele     += 0.5  * my_square(syn_dif_time - obs_dif_time);   // because a pair sf sources are counted twice, thus * 0.5
                    } else{
                        obj_cr_dif      += 0.5 * my_square(syn_dif_time - obs_dif_time)*data.weight;
                        misfit_cr_dif   += 0.5  * my_square(syn_dif_time - obs_dif_time);
                    }

                } else if (name_sim_src == name_src2) {
                    CUSTOMREAL syn_dif_time = IP.syn_time_list_sr[name_src2][name_rec] - IP.syn_time_list_sr[name_src1][name_rec];
                    CUSTOMREAL obs_dif_time = - data.cr_dif_travel_time_obs;
                    CUSTOMREAL adjoint_source = IP.get_rec_list_nv(name_rec).adjoint_source + (syn_dif_time - obs_dif_time)*data.weight;
                    IP.set_adjoint_source(name_rec, adjoint_source);
                    
                    // contribute misfit
                    obj     += 0.5 * my_square(syn_dif_time - obs_dif_time)*data.weight;
                    misfit  += 0.5  * my_square(syn_dif_time - obs_dif_time);
                    if (IP.get_src_list_nv(name_src1).is_out_of_region || \
                        IP.get_src_list_nv(name_src2).is_out_of_region || \
                        IP.get_rec_list_nv(name_rec).is_out_of_region){
                        obj_tele        +=  0.5 * my_square(syn_dif_time - obs_dif_time)*data.weight;
                        misfit_tele     +=  0.5  * my_square(syn_dif_time - obs_dif_time);
                    } else{
                        obj_cr_dif      +=  0.5 * my_square(syn_dif_time - obs_dif_time)*data.weight;
                        misfit_cr_dif   +=  0.5  * my_square(syn_dif_time - obs_dif_time);
                    }

                } else {
                    std::cout << "error match of data in function: calculate_adjoint_source_nv() " << std::endl;
                }

            } else if (data.is_rec_pair) {  // common source differential traveltime
                std::string name_src = data.name_src_single;
                std::string name_rec1 = data.name_rec_pair[0];
                std::string name_rec2  = data.name_rec_pair[1];
                CUSTOMREAL syn_dif_time = IP.syn_time_list_sr[name_src][name_rec1] - IP.syn_time_list_sr[name_src][name_rec2];
                CUSTOMREAL obs_dif_time = data.cs_dif_travel_time_obs;

                CUSTOMREAL adjoint_source = IP.get_rec_list_nv(name_rec1).adjoint_source + (syn_dif_time - obs_dif_time)*data.weight;
                IP.set_adjoint_source(name_rec1, adjoint_source);
                adjoint_source = IP.get_rec_list_nv(name_rec2).adjoint_source - (syn_dif_time - obs_dif_time)*data.weight;
                IP.set_adjoint_source(name_rec2, adjoint_source);

                // contribute misfit
                obj     += 1.0 * my_square(syn_dif_time - obs_dif_time)*data.weight;
                misfit  += 1.0 * my_square(syn_dif_time - obs_dif_time);
                if (IP.get_src_list_nv(name_src).is_out_of_region || \
                    IP.get_rec_list_nv(name_rec1).is_out_of_region || \
                    IP.get_rec_list_nv(name_rec2).is_out_of_region){
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
    // broadcast_cr_sub(allsum,2,0);

    return allsum;
    
}


// void Receiver::interpolate_travel_time(Grid& grid, SrcRec& rec) {
//     // calculate the travel time of the receiver by 3d linear interpolation

//     // copy some parameters
//     CUSTOMREAL delta_lon = grid.get_delta_lon();
//     CUSTOMREAL delta_lat = grid.get_delta_lat();
//     CUSTOMREAL delta_r   = grid.get_delta_r();

//     // store receiver position in radian
//     CUSTOMREAL rec_lon = rec.lon*DEG2RAD;
//     CUSTOMREAL rec_lat = rec.lat*DEG2RAD;
//     CUSTOMREAL rec_r = depth2radius(rec.dep); // r in km

//     // check if the receiver is in this subdomain
//     bool is_in_subdomain = false;
//     if (grid.get_lon_min_loc() <= rec_lon && rec_lon < grid.get_lon_max_loc() && 
//         grid.get_lat_min_loc() <= rec_lat && rec_lat < grid.get_lat_max_loc() && 
//         grid.get_r_min_loc()   <= rec_r   && rec_r   < grid.get_r_max_loc()   ) {
//         is_in_subdomain = true;
//     }

//     // check the rank where the source is located
//     int rec_rank = -1;
//     bool* rec_flags = new bool[nprocs];
//     allgather_bool_single(&is_in_subdomain, rec_flags);
//     for (int i = 0; i < nprocs; i++) {
//         if (rec_flags[i]) {
//             rec_rank = i;
//             break; // this break means that the first subdomain is used if the receiver is in multiple subdomains (ghost layer)
//         }
//     }
//     delete[] rec_flags;
//      // check if the receiver is in the global domain
//     if (rec_rank == -1) {
//         std::cout << "Error: the receiver is not in the global domain" << std::endl;
//         // print rec
//         std::cout << "id_src: " << rec.id_src << " id_rec: " << rec.id_rec << " n_rec: " << rec.n_rec << " depth: " << rec.dep << " lat: " << rec.lat << " lon: " << rec.lon << " arr_time: " << rec.arr_time << std::endl;
//         // print boundary
//         //std::cout << "lon min max rec: " << grid.get_lon_min_loc() << " " << grid.get_lon_max_loc() << " " << rec_lon << std::endl;
//         //std::cout << "lat min max rec: " << grid.get_lat_min_loc() << " " << grid.get_lat_max_loc() << " " << rec_lat << std::endl;
//         //std::cout << "r min max rec: " << grid.get_r_min_loc() << " " << grid.get_r_max_loc() << " " << rec_r << std::endl;

//         std::cout << "lon+bound min max rec: " << (grid.get_lon_min_loc() - delta_lon)*RAD2DEG     << " " << (grid.get_lon_max_loc() + delta_lon)*RAD2DEG     << " " << rec_lon*RAD2DEG    << std::endl;
//         std::cout << "lat+bound min max rec: " << (grid.get_lat_min_loc() - delta_lat)*RAD2DEG     << " " << (grid.get_lat_max_loc() + delta_lat)*RAD2DEG     << " " << rec_lat*RAD2DEG    << std::endl;
//         std::cout << "r+bound min max rec: "   << radius2depth(grid.get_r_min_loc()   - delta_r  ) << " " << radius2depth(grid.get_r_max_loc()   + delta_r  ) << " " << radius2depth(rec_r)<< std::endl;
//         exit(1);
//     }

//     CUSTOMREAL vinterp = 0.0;

//     if (is_in_subdomain) {
//         // calculate the interpolated travel time and broadcast it

//         // descretize source position (LOCAL) ID)
//         int i_rec = std::floor((rec_lon - grid.get_lon_min_loc()) / delta_lon);
//         int j_rec = std::floor((rec_lat - grid.get_lat_min_loc()) / delta_lat);
//         int k_rec = std::floor((rec_r   - grid.get_r_min_loc())   / delta_r);

//         // discretized receiver position
//         CUSTOMREAL dis_rec_lon = grid.p_loc_1d[i_rec];
//         CUSTOMREAL dis_rec_lat = grid.t_loc_1d[j_rec];
//         CUSTOMREAL dis_rec_r   = grid.r_loc_1d[k_rec];

//         // relative position errors
//         CUSTOMREAL e_lon = std::min({_1_CR,(rec_lon - dis_rec_lon)/delta_lon});
//         CUSTOMREAL e_lat = std::min({_1_CR,(rec_lat - dis_rec_lat)/delta_lat});
//         CUSTOMREAL e_r   = std::min({_1_CR,(rec_r   - dis_rec_r)  /delta_r});

// //        if(if_verbose){
// //            std::cout << "(rec_lon - dis_rec_lon)/dlon: " << (rec_lon - dis_rec_lon)/delta_lon << std::endl;
// //            std::cout << "(rec_lat - dis_rec_lat)/dlat: " << (rec_lat - dis_rec_lat)/delta_lat << std::endl;
// //            std::cout << "(rec_r   - dis_rec_r  )/dr: "   << (rec_r   - dis_rec_r  )/delta_r << std::endl;
// //            std::cout << "rec_lon, dis_rec_lon, delta_lon: " << rec_lon << " " << dis_rec_lon << " " << delta_lon << std::endl;
// //            std::cout << "rec_lat, dis_rec_lat, delta_lat: " << rec_lat << " " << dis_rec_lat << " " << delta_lat << std::endl;
// //            std::cout << "rec_r, dis_rec_r, delta_r  : " << rec_r << ", " << dis_rec_r << ", " <<  delta_r  << std::endl;
// //            std::cout << "loc_K. k_rec: " << loc_K << "," << k_rec << std::endl;
// //            std::cout << "r_loc_1d[loc_K-1]: " << grid.r_loc_1d[loc_K-1] << std::endl;
// //        }

//         int i_rec_p1 = i_rec + 1;
//         int j_rec_p1 = j_rec + 1;
//         int k_rec_p1 = k_rec + 1;

//         // exclude the points if they are out of the domain
//         if (i_rec_p1 > loc_I-1 
//          || j_rec_p1 > loc_J-1 
//          || k_rec_p1 > loc_K-1) {
//             // exit(1) as the source is out of the domain
//             std::cout << "Error: the receiver is out of the domain" << std::endl;
//             std::cout << "id_src: " << rec.id_src << " id_rec: " << rec.id_rec << " n_rec: " << rec.n_rec << " depth: " << rec.dep << " lat: " << rec.lat << " lon: " << rec.lon << " arr_time: " << rec.arr_time << std::endl;
//             std::cout << "lon min max rec: " << grid.get_lon_min_loc()*RAD2DEG << " " << grid.get_lon_max_loc()*RAD2DEG << " " << rec_lon*RAD2DEG << std::endl;
//             std::cout << "lat min max rec: " << grid.get_lat_min_loc()*RAD2DEG << " " << grid.get_lat_max_loc()*RAD2DEG << " " << rec_lat*RAD2DEG << std::endl;
//             std::cout << "r min max rec: " << radius2depth(grid.get_r_min_loc()) << " " << radius2depth(grid.get_r_max_loc()) << " " << radius2depth(rec_r) << std::endl;
//             std::cout << "i_rec: " << i_rec << " j_rec: " << j_rec << " k_rec: " << k_rec << std::endl;
//             std::cout << "i_rec_p1: " << i_rec_p1 << " j_rec_p1: " << j_rec_p1 << " k_rec_p1: " << k_rec_p1 << std::endl;
//             std::cout << "loc_I-1: " << loc_I-1 << " loc_J-1: " << loc_J-1 << " loc_K-1: " << loc_K-1 << std::endl;
//             //finalize_mpi();
//             exit(1);
//          }


//         vinterp = (_1_CR - e_lon) * (_1_CR - e_lat) * (_1_CR - e_r) * grid.T_loc[I2V(i_rec,   j_rec,   k_rec)]   
//                 +          e_lon  * (_1_CR - e_lat) * (_1_CR - e_r) * grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec)]   
//                 + (_1_CR - e_lon) *          e_lat  * (_1_CR - e_r) * grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec)]   
//                 + (_1_CR - e_lon) * (_1_CR - e_lat) *          e_r  * grid.T_loc[I2V(i_rec,   j_rec,   k_rec_p1)] 
//                 +          e_lon  *          e_lat  * (_1_CR - e_r) * grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec)]   
//                 +          e_lon  * (_1_CR - e_lat) *          e_r  * grid.T_loc[I2V(i_rec_p1,j_rec,   k_rec_p1)] 
//                 + (_1_CR - e_lon) *          e_lat  *          e_r  * grid.T_loc[I2V(i_rec,   j_rec_p1,k_rec_p1)] 
//                 +          e_lon  *          e_lat  *          e_r  * grid.T_loc[I2V(i_rec_p1,j_rec_p1,k_rec_p1)];

//         //std::cout << "DEBUG near and vinterp : " << grid.T_loc[I2V(i_rec,j_rec,k_rec)] << ", " << vinterp << std::endl;

//         // broadcast interpolated travel time
//         broadcast_cr_single(vinterp, rec_rank);

//     } else {
//         // receive the calculated traveltime
//         broadcast_cr_single(vinterp, rec_rank);
//     }

//     // store the calculated travel time
//     rec.arr_time = vinterp;

// }

void Receiver::interpolate_travel_time_nv(Grid& grid, InputParams& IP, std::string name_src, std::string name_rec) {
    // calculate the travel time of the receiver by 3d linear interpolation

    SrcRecInfo rec = IP.get_rec_point_nv(name_rec);

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

    // store the calculated travel time
    IP.syn_time_list_sr[name_src][name_rec] = vinterp;
}



void Receiver::init_vars_src_reloc(InputParams& IP){
    if (subdom_main) {
   
        // // calculate gradient of travel time at each receiver (swapped source)
        // for (auto iter = IP.rec_list_nv.begin(); iter != IP.rec_list_nv.end(); iter++) {
        //     iter->second.tau_opt                    = _0_CR;
        //     iter->second.grad_chi_i                 = _0_CR;
        //     iter->second.grad_chi_j                 = _0_CR;
        //     iter->second.grad_chi_k                 = _0_CR;
        //     iter->second.sum_weight                 = _0_CR;
        //     iter->second.vobj_src_reloc_old         = iter->second.vobj_src_reloc;
        //     iter->second.vobj_src_reloc             = _0_CR;
        //     iter->second.vobj_grad_norm_src_reloc   = _0_CR;
        // }

        // initialize the kernel of unrelocated source
        for (int i = 0; i < (int)IP.name_for_reloc.size(); i++){
            std::string name_rec = IP.name_for_reloc[i];
            if (IP.is_ortime_local_search == 0) // origin time global search
                IP.rec_list_nv[name_rec].tau_opt                    = _0_CR;
            else    // origin time local search
                IP.rec_list_nv[name_rec].grad_tau                   = _0_CR;

            IP.rec_list_nv[name_rec].grad_chi_i                 = _0_CR;
            IP.rec_list_nv[name_rec].grad_chi_j                 = _0_CR;
            IP.rec_list_nv[name_rec].grad_chi_k                 = _0_CR;
            IP.rec_list_nv[name_rec].sum_weight                 = _0_CR;
            IP.rec_list_nv[name_rec].vobj_src_reloc_old         = IP.rec_list_nv[name_rec].vobj_src_reloc;
            IP.rec_list_nv[name_rec].vobj_src_reloc             = _0_CR;
            IP.rec_list_nv[name_rec].vobj_grad_norm_src_reloc   = _0_CR;
        }

    }
}

void Receiver::calculate_T_gradient(InputParams& IP, Grid& grid){

    if(subdom_main){
        // get list of receivers from input parameters
        // std::vector<SrcRec>& receivers = IP.get_rec_points(id_sim_src);

        // calculate gradient of travel time at each receiver (swapped source)
        for (int i = 0; i < (int)IP.data_info_smap_reloc[name_sim_src].size(); i++){
            std::string name_rec = IP.data_info_smap_reloc[name_sim_src][i].name_rec;
            if(!IP.rec_list_nv[name_rec].is_stop)
                calculate_T_gradient_one_rec(grid, IP.rec_list_nv[name_rec]);
        }

        // for(auto iter = IP.rec_list_nv.begin(); iter != IP.rec_list_nv.end(); iter++){
        //     std::cout << "DTi (lon) is: " << iter->second.DTi << std::endl;
        //     std::cout << "DTj (lat) is: " << iter->second.DTj << std::endl;
        //     std::cout << "DTk (dep) is: " << iter->second.DTk << std::endl;
        // }
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

        // // exclude the points if they are out of the domain
        // if (i_rec_p1 > loc_I-1 
        //  || j_rec_p1 > loc_J-1 
        //  || k_rec_p1 > loc_K-1) {
        //     // exit(1) as the source is out of the domain
        //     std::cout << "Error: the receiver is out of the domain" << std::endl;
        //     std::cout << " id_rec: " << rec.id << " name: " << rec.name << " depth: " << rec.dep << " lat: " << rec.lat << " lon: " << rec.lon << std::endl;
        //     std::cout << "lon min max rec: " << grid.get_lon_min_loc()*RAD2DEG << " " << grid.get_lon_max_loc()*RAD2DEG << " " << rec_lon*RAD2DEG << std::endl;
        //     std::cout << "lat min max rec: " << grid.get_lat_min_loc()*RAD2DEG << " " << grid.get_lat_max_loc()*RAD2DEG << " " << rec_lat*RAD2DEG << std::endl;
        //     std::cout << "r min max rec: " << radius2depth(grid.get_r_min_loc()) << " " << radius2depth(grid.get_r_max_loc()) << " " << radius2depth(rec_r) << std::endl;
        //     std::cout << "i_rec: " << i_rec << " j_rec: " << j_rec << " k_rec: " << k_rec << std::endl;
        //     std::cout << "i_rec_p1: " << i_rec_p1 << " j_rec_p1: " << j_rec_p1 << " k_rec_p1: " << k_rec_p1 << std::endl;
        //     std::cout << "loc_I-1: " << loc_I-1 << " loc_J-1: " << loc_J-1 << " loc_K-1: " << loc_K-1 << std::endl;
        //     //finalize_mpi();
        //     exit(1);
        // }

        // exclude the points if they are adjancet to the domain boundary
        if (i_rec_p1 > loc_I-2 || j_rec_p1 > loc_J-2 || k_rec_p1 > loc_K-2 || \
            i_rec    < 1       || j_rec    < 1       || k_rec    < 1 ) {
            // exit(1) as the source is out of the domain
            std::cout << "Error: the location is close to the domain boundary" << std::endl;
            std::cout << " id_rec: " << rec.id << " name: " << rec.name << " depth: " << rec.dep << " lat: " << rec.lat << " lon: " << rec.lon << std::endl;
            std::cout << "lon min max rec: " << grid.get_lon_min_loc()*RAD2DEG << " " << grid.get_lon_max_loc()*RAD2DEG << " " << rec_lon*RAD2DEG << std::endl;
            std::cout << "lat min max rec: " << grid.get_lat_min_loc()*RAD2DEG << " " << grid.get_lat_max_loc()*RAD2DEG << " " << rec_lat*RAD2DEG << std::endl;
            std::cout << "r min max rec: " << radius2depth(grid.get_r_min_loc()) << " " << radius2depth(grid.get_r_max_loc()) << " " << radius2depth(rec_r) << std::endl;
            std::cout << "i_rec: " << i_rec << " j_rec: " << j_rec << " k_rec: " << k_rec << std::endl;
            std::cout << "i_rec_p1: " << i_rec_p1 << " j_rec_p1: " << j_rec_p1 << " k_rec_p1: " << k_rec_p1 << std::endl;
            std::cout << "loc_I-1: " << loc_I-1 << " loc_J-1: " << loc_J-1 << " loc_K-1: " << loc_K-1 << std::endl;
            //finalize_mpi();
            // exit(1);
            rec.is_stop = true;          
        }

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

        if (!rec.is_stop){
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
        } else {
            DTk = 0.0;
            DTj = 0.0;
            DTi = 0.0;
        }

        
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
void Receiver::calculate_optimal_origin_time(InputParams& IP){
    if (subdom_main) {
        // get list of receivers from input parameters
        // std::vector<SrcRec>& receivers = IP.get_rec_points(id_sim_src);
        std::vector<DataInfo> data_info_tmp = IP.data_info_smap_reloc[name_sim_src];

        // calculate gradient of travel time at each receiver (swapped source)
        for (int i = 0; i < (int)data_info_tmp.size(); i++){
            std::string name_rec = data_info_tmp[i].name_rec;
            if (! IP.rec_list_nv[name_rec].is_stop){
                CUSTOMREAL weight = data_info_tmp[i].weight;
                CUSTOMREAL misfit = data_info_tmp[i].travel_time_obs - IP.syn_time_list_sr[name_sim_src][name_rec];

                
                if (IP.is_ortime_local_search == 0){    // global search
                    IP.rec_list_nv[name_rec].tau_opt += misfit * weight;
                    IP.rec_list_nv[name_rec].sum_weight += weight;
                } else {    // local search
                    IP.rec_list_nv[name_rec].grad_tau += weight * 
                        (IP.syn_time_list_sr[name_sim_src][name_rec] + IP.rec_list_nv[name_rec].tau_opt - data_info_tmp[i].travel_time_obs);
                }
                  
            }
        }

    }
}

void Receiver::divide_optimal_origin_time_by_summed_weight(InputParams& IP) {
    if (subdom_main) {  

        for (int i = 0; i < (int)IP.name_for_reloc.size(); i++) {
            std::string name_rec = IP.name_for_reloc[i];
            allreduce_cr_sim_single_inplace(IP.rec_list_nv[name_rec].tau_opt);
            allreduce_cr_sim_single_inplace(IP.rec_list_nv[name_rec].sum_weight);
            IP.rec_list_nv[name_rec].tau_opt /= IP.rec_list_nv[name_rec].sum_weight;
            // std::cout << "id_sim" << id_sim << ", name: " << iter->first << ", ortime: " << iter->second.tau_opt <<std::endl;
        }
    }
    synchronize_all_world();
}

void Receiver::calculate_obj_reloc(InputParams& IP, int i_iter){

    for (long unsigned int i_src = 0; i_src < IP.src_names_this_sim.size(); i_src++) {
        // load the global id of this src
        id_sim_src = IP.src_ids_this_sim[i_src]; // local src id to global src id
        name_sim_src = IP.src_names_this_sim[i_src]; // local src name to global src id

        if (subdom_main) {
            // get list of receivers from input parameters
            std::vector<DataInfo> data_info_tmp = IP.data_info_smap_reloc[name_sim_src];

            // calculate objective function at each receiver (swapped source)
            for (int i = 0; i < (int)data_info_tmp.size(); i++){
                std::string name_rec = data_info_tmp[i].name_rec;
                if (! IP.rec_list_nv[name_rec].is_stop){
                    CUSTOMREAL weight = data_info_tmp[i].weight;
                    CUSTOMREAL misfit = IP.syn_time_list_sr[name_sim_src][name_rec] - data_info_tmp[i].travel_time_obs + IP.rec_list_nv[name_rec].tau_opt;
                    IP.rec_list_nv[name_rec].vobj_src_reloc += weight / _2_CR * my_square(misfit);
               }
                
                 
            }
        }
    }
    
    // sum the obj from all sources (swapped receivers)
    if (subdom_main) {
        for (int i = 0; i < (int)IP.name_for_reloc.size(); i++){
            std::string name_rec = IP.name_for_reloc[i];
            allreduce_cr_sim_single_inplace(IP.rec_list_nv[name_rec].vobj_src_reloc);
        }
    }
    synchronize_all_world();

    for (int i = 0; i < (int)IP.name_for_reloc.size(); i++){
        std::string name_rec = IP.name_for_reloc[i];
        CUSTOMREAL obj = IP.rec_list_nv[name_rec].vobj_src_reloc;
        CUSTOMREAL old_obj = IP.rec_list_nv[name_rec].vobj_src_reloc_old;
        if (i_iter != 0 && old_obj < obj){    // if obj increase, decrease the step length of this (swapped) source
            IP.rec_list_nv[name_rec].step_length_max *= step_length_decay;
        }
    }    
}


// calculate the gradient of the objective function
void Receiver::calculate_grad_obj_src_reloc(InputParams& IP) {

    if(subdom_main){
        // get list of receivers from input parameters
        // std::vector<SrcRec>& receivers = IP.get_rec_points(id_sim_src);

        // get the source weight
        // CUSTOMREAL src_weight = IP.get_src_point(id_sim_src).weight;

        // calculate gradient of travel time at each receiver (swapped source)
        for(int i = 0; i < (int)IP.data_info_smap_reloc[name_sim_src].size(); i++){
            std::string name_rec = IP.data_info_smap_reloc[name_sim_src][i].name_rec;
            if(!IP.rec_list_nv[name_rec].is_stop){
                CUSTOMREAL misfit = IP.syn_time_list_sr[name_sim_src][name_rec] - IP.data_info_smap_reloc[name_sim_src][i].travel_time_obs + IP.rec_list_nv[name_rec].tau_opt;
                CUSTOMREAL weight = IP.data_info_smap_reloc[name_sim_src][i].weight;
                IP.rec_list_nv[name_rec].grad_chi_k += (misfit) * IP.rec_list_nv[name_rec].DTk * weight;
                IP.rec_list_nv[name_rec].grad_chi_j += (misfit) * IP.rec_list_nv[name_rec].DTj * weight;
                IP.rec_list_nv[name_rec].grad_chi_i += (misfit) * IP.rec_list_nv[name_rec].DTi * weight;
            }
            // std::cout << "misfit: " << misfit << "tau_opt: " << IP.rec_list_nv[name_rec].tau_opt << std::endl;   
            // std::cout << "id_sim: " << id_sim << ", src: " << name_sim_src << ", misfit: " << misfit << "tau_opt: " << IP.rec_list_nv[name_rec].tau_opt 
            //               << std::endl;
            // std::cout << "DTk, DTj, DTi = " << IP.rec_list_nv[name_rec].DTk << ", " << IP.rec_list_nv[name_rec].DTj 
            //           << ", " << IP.rec_list_nv[name_rec].DTi << std::endl;            
        }

        // for(auto iter = IP.rec_list_nv.begin(); iter != IP.rec_list_nv.end(); iter++){
            // std::cout << "grad_chi_i (lon,rad) is: " << iter->second.grad_chi_i << std::endl;
            // std::cout << "grad_chi_j (lat,rad) is: " << iter->second.grad_chi_j << std::endl;
            // std::cout << "grad_chi_k (radius,km) is: " << iter->second.grad_chi_k << std::endl << std::endl;
        // }
    }
    
}

void Receiver::update_source_location(InputParams& IP, Grid& grid) {

    if (subdom_main) {
        // get list of receivers from input parameters
        // std::vector<SrcRec>& receivers = IP.get_rec_points(id_sim_src);

        for(int i = 0; i < (int)IP.name_for_reloc.size(); i++){
            std::string name_rec = IP.name_for_reloc[i];
            
            if (IP.rec_list_nv[name_rec].is_stop){
                // do nothing
            } else {
                // rescale gradient
                CUSTOMREAL grad_dep_km;
                if (abs(IP.max_change_dep - abs(IP.rec_list_nv[name_rec].change_dep)) < 0.001 || abs(IP.rec_list_nv[name_rec].dep) < 0.001)
                    grad_dep_km = 0.0;
                else
                    grad_dep_km = - IP.rec_list_nv[name_rec].grad_chi_k;
                CUSTOMREAL grad_lat_km;
                if (abs(IP.max_change_lat - abs(IP.rec_list_nv[name_rec].change_lat)) < 0.001)
                    grad_lat_km = 0.0;
                else
                    grad_lat_km = IP.rec_list_nv[name_rec].grad_chi_j/(R_earth);
                CUSTOMREAL grad_lon_km;
                if (abs(IP.max_change_lon - abs(IP.rec_list_nv[name_rec].change_lon)) < 0.001)
                    grad_lon_km = 0.0;
                else
                    grad_lon_km = IP.rec_list_nv[name_rec].grad_chi_i/(R_earth * cos(IP.rec_list_nv[name_rec].lat * DEG2RAD));
                
                CUSTOMREAL norm_grad;
                CUSTOMREAL grad_ortime = 0.0;
                if (IP.is_ortime_local_search == 0){
                    norm_grad   = std::sqrt(my_square(grad_dep_km) + my_square(grad_lat_km) + my_square(grad_lon_km));
                } else {
                    if (abs(IP.max_change_ortime - abs(IP.rec_list_nv[name_rec].change_ortime)) < 0.001)
                        grad_ortime = 0.0;
                    else
                        grad_ortime = IP.rec_list_nv[name_rec].grad_tau/IP.ref_ortime_change;
                    
                    norm_grad   = std::sqrt(my_square(grad_dep_km) + my_square(grad_lat_km) + my_square(grad_lon_km) + my_square(grad_ortime));
                } 

                if (norm_grad < TOL_SRC_RELOC){
                    IP.rec_list_nv[name_rec].is_stop = true;
                    continue;
                }


                CUSTOMREAL step_length;
                step_length = 0.5 * IP.rec_list_nv[name_rec].vobj_src_reloc/my_square(norm_grad);
                

                // rescale update for dep, lat, lon
                CUSTOMREAL update_max = -1.0;
                CUSTOMREAL update_dep = step_length * grad_dep_km;
                update_max = std::max(update_max,abs(update_dep));
                CUSTOMREAL update_lat = step_length * grad_lat_km;
                update_max = std::max(update_max,abs(update_lat));
                CUSTOMREAL update_lon = step_length * grad_lon_km;
                update_max = std::max(update_max,abs(update_lon));
                
                CUSTOMREAL update_ortime = 0.0;
                if (IP.is_ortime_local_search == 1){
                    update_ortime = step_length * grad_ortime * IP.step_length_ortime_rescale;
                    update_max = std::max(update_max,abs(update_ortime));
                }


                CUSTOMREAL downscale = 1.0;
                if (update_max > IP.rec_list_nv[name_rec].step_length_max){
                    downscale = IP.rec_list_nv[name_rec].step_length_max / update_max;
                }

                // update value for dep (km), lat (degree), lon (degree)
                IP.rec_list_nv[name_rec].vobj_grad_norm_src_reloc = norm_grad;
                update_dep = - update_dep * downscale;
                update_lat = - update_lat * downscale / R_earth * RAD2DEG;
                update_lon = - update_lon * downscale / (R_earth * cos(IP.rec_list_nv[name_rec].lat * DEG2RAD)) * RAD2DEG;

                if (IP.is_ortime_local_search == 1){
                    update_ortime = - update_ortime * downscale * IP.step_length_ortime_rescale;
                }

                // limit the update for dep, lat, lon
                if (abs(IP.rec_list_nv[name_rec].change_dep + update_dep) > IP.max_change_dep){
                    if (IP.rec_list_nv[name_rec].change_dep + update_dep > 0)
                        update_dep = IP.max_change_dep - IP.rec_list_nv[name_rec].change_dep;
                    else 
                        update_dep = - IP.max_change_dep - IP.rec_list_nv[name_rec].change_dep;
                    
                } 
                if (abs(IP.rec_list_nv[name_rec].change_lat + update_lat) > IP.max_change_lat){
                    if (IP.rec_list_nv[name_rec].change_dep + update_dep > 0)
                        update_lat = IP.max_change_lat - IP.rec_list_nv[name_rec].change_lat;
                    else 
                        update_lat = - IP.max_change_lat - IP.rec_list_nv[name_rec].change_lat;
                } 
                if (abs(IP.rec_list_nv[name_rec].change_lon + update_lon) > IP.max_change_lon){
                    if (IP.rec_list_nv[name_rec].change_dep + update_dep > 0)
                        update_lon = IP.max_change_lon - IP.rec_list_nv[name_rec].change_lon;
                    else 
                        update_lon = - IP.max_change_lon - IP.rec_list_nv[name_rec].change_lon;
                } 

                if (IP.is_ortime_local_search == 1){
                    if (abs(IP.rec_list_nv[name_rec].change_ortime + update_ortime) > IP.max_change_ortime){
                        if (IP.rec_list_nv[name_rec].change_ortime + update_ortime > 0)
                            update_ortime = IP.max_change_ortime - IP.rec_list_nv[name_rec].change_ortime;
                        else 
                            update_ortime = - IP.max_change_ortime - IP.rec_list_nv[name_rec].change_ortime;
                    } 
                }

                // earthquake should be below the surface
                if (IP.rec_list_nv[name_rec].dep + update_dep < 0){
                    update_dep = - IP.rec_list_nv[name_rec].dep;
                }

                IP.rec_list_nv[name_rec].dep += update_dep;
                IP.rec_list_nv[name_rec].lat += update_lat;
                IP.rec_list_nv[name_rec].lon += update_lon;

                IP.rec_list_nv[name_rec].change_dep += update_dep;
                IP.rec_list_nv[name_rec].change_lat += update_lat;
                IP.rec_list_nv[name_rec].change_lon += update_lon;

                if (IP.is_ortime_local_search == 1){
                    IP.rec_list_nv[name_rec].tau_opt += update_ortime;
                    IP.rec_list_nv[name_rec].change_ortime += update_ortime;
                }

                if (IP.rec_list_nv[name_rec].step_length_max < TOL_STEP_SIZE){
                    IP.rec_list_nv[name_rec].is_stop = true;
                }

                if(id_sim == 0 && myrank == 0 && (name_rec == "1753590" || std::isnan(IP.rec_list_nv[name_rec].tau_opt))){
                    std::cout << "name_rec: " << name_rec << ", change_dep: " << IP.rec_list_nv[name_rec].change_dep
                              << ", change_lat: " << IP.rec_list_nv[name_rec].change_lat
                              << ", change_lon: " << IP.rec_list_nv[name_rec].change_lon
                              << ", change_ortime: " << IP.rec_list_nv[name_rec].change_ortime
                              << std::endl;

                    std::cout << "kernel k(dep,km): " << grad_dep_km
                          << ", kernel i(lat,km): " << grad_lat_km
                          << ", kernel j(lon,km): " << grad_lon_km
                          << ", kernel time ref: " << grad_ortime
                          << std::endl; 

                    std::cout << "update dep (km): " << update_dep
                          << ", update lat (km): " << update_lat * DEG2RAD * R_earth 
                          << ", update lon (km): " << update_lon * DEG2RAD * R_earth * cos(IP.rec_list_nv[name_rec].lat * DEG2RAD)
                          << ", update ortime " << update_ortime
                          << std::endl; 
                    std::cout << "norm_grad: " << norm_grad << std::endl;     
                }

            

                // check if the new receiver position is within the domain
                // if not then set the receiver position to the closest point on the domain

                CUSTOMREAL mergin_lon = 1.01 * grid.get_delta_lon();
                CUSTOMREAL mergin_lat = 1.01 * grid.get_delta_lat();
                CUSTOMREAL mergin_r   = 1.01 * grid.get_delta_r();

                if (IP.rec_list_nv[name_rec].lon < IP.get_min_lon()*RAD2DEG)
                    IP.rec_list_nv[name_rec].lon = IP.get_min_lon()*RAD2DEG + mergin_lon;
                if (IP.rec_list_nv[name_rec].lon > IP.get_max_lon()*RAD2DEG)
                    IP.rec_list_nv[name_rec].lon = IP.get_max_lon()*RAD2DEG - mergin_lon;
                if (IP.rec_list_nv[name_rec].lat < IP.get_min_lat()*RAD2DEG)
                    IP.rec_list_nv[name_rec].lat = IP.get_min_lat()*RAD2DEG + mergin_lat;
                if (IP.rec_list_nv[name_rec].lat > IP.get_max_lat()*RAD2DEG)
                    IP.rec_list_nv[name_rec].lat = IP.get_max_lat()*RAD2DEG - mergin_lat;
                if (IP.rec_list_nv[name_rec].dep < IP.get_min_dep())
                    IP.rec_list_nv[name_rec].dep = IP.get_min_dep() + mergin_r;
                if (IP.rec_list_nv[name_rec].dep > IP.get_max_dep())
                    IP.rec_list_nv[name_rec].dep = IP.get_max_dep() - mergin_r;
            }
        }
    }
}


// void Receiver::update_source_location(InputParams& IP, Grid& grid, std::vector<SrcRec>& unique_rec_list) {

//     if (subdom_main) {
//         // get list of receivers from input parameters
//         // std::vector<SrcRec>& receivers = IP.get_rec_points(id_sim_src);

//         for(auto iter = rec_list_nv.begin(); iter != rec_list_nv.end(); iter++){
//             iter->second.dep -= 
//         }

//         for (auto& rec: receivers) {
//             SrcRec& rec_unique = unique_rec_list[rec.id_unique_list];
//             rec.dep -= rec_unique.grad_chi_k * step_length_src_reloc;
//             rec.lat -= rec_unique.grad_chi_j * step_length_src_reloc;
//             rec.lon -= rec_unique.grad_chi_i * step_length_src_reloc;

//             // check if the new receiver position is within the domain
//             // if not then set the receiver position to the closest point on the domain

//             CUSTOMREAL mergin_lon = 1.01 * grid.get_delta_lon();
//             CUSTOMREAL mergin_lat = 1.01 * grid.get_delta_lat();
//             CUSTOMREAL mergin_r   = 1.01 * grid.get_delta_r();

//             if (rec.lon < IP.get_min_lon()*RAD2DEG)
//                 rec.lon = IP.get_min_lon()*RAD2DEG + mergin_lon;
//             if (rec.lon > IP.get_max_lon()*RAD2DEG)
//                 rec.lon = IP.get_max_lon()*RAD2DEG - mergin_lon;
//             if (rec.lat < IP.get_min_lat()*RAD2DEG)
//                 rec.lat = IP.get_min_lat()*RAD2DEG + mergin_lat;
//             if (rec.lat > IP.get_max_lat()*RAD2DEG)
//                 rec.lat = IP.get_max_lat()*RAD2DEG - mergin_lat;
//             if (rec.dep < IP.get_min_dep())
//                 rec.dep = IP.get_min_dep() + mergin_r;
//             if (rec.dep > IP.get_max_dep())
//                 rec.dep = IP.get_max_dep() - mergin_r;
//        }
//     }
// }