#include "source.h"

Source::Source(InputParams &IP, Grid &grid, bool& is_teleseismic, const std::string& name_sim_src) {

    if (is_teleseismic) return;

    if (subdom_main) {
        if(if_verbose) stdout_by_main("--- start source initialization ---");

        // copy some parameters
        delta_lon = grid.get_delta_lon();
        delta_lat = grid.get_delta_lat();
        delta_r   = grid.get_delta_r();

        // set source position
        src_lon = IP.get_src_lon(   name_sim_src); // in radian
        src_lat = IP.get_src_lat(   name_sim_src); // in radian
        src_r   = IP.get_src_radius(name_sim_src); // radious

        // descretize source position (LOCAL ID)
        i_src_loc = std::floor((src_lon - grid.get_lon_min_loc()) / grid.get_delta_lon());
        j_src_loc = std::floor((src_lat - grid.get_lat_min_loc()) / grid.get_delta_lat());
        k_src_loc = std::floor((src_r   - grid.get_r_min_loc())   / grid.get_delta_r())  ;

        // check if the source is in this subdomain (including the ghost nodes)
        if (grid.get_lon_min_loc() <= src_lon && src_lon <= grid.get_lon_max_loc()  && \
            grid.get_lat_min_loc() <= src_lat && src_lat <= grid.get_lat_max_loc()  && \
            grid.get_r_min_loc()   <= src_r   && src_r   <= grid.get_r_max_loc()   ) {
            is_in_subdomain = true;
        }

        // check the rank where the source is located
        src_flags = new bool[nprocs];
        allgather_bool_single(&is_in_subdomain, src_flags);
        for (int i = 0; i < nprocs; i++) {
            if (src_flags[i]) {
                src_rank = i;
                break;
            }
        }

        delete[] src_flags;


        if (is_in_subdomain){
            // discretized source position (Global)
            dis_src_lon = grid.p_loc_1d[i_src_loc];
            dis_src_lat = grid.t_loc_1d[j_src_loc];
            dis_src_r   = grid.r_loc_1d[k_src_loc];

            // position error
            error_lon = (src_lon - dis_src_lon);
            error_lat = (src_lat - dis_src_lat);
            error_r   = (src_r   - dis_src_r);

            // relative position error
            dis_src_err_lon = std::min(error_lon / grid.get_delta_lon(), _1_CR);
            dis_src_err_lat = std::min(error_lat / grid.get_delta_lat(), _1_CR);
            dis_src_err_r   = std::min(error_r   / grid.get_delta_r()  , _1_CR);

            // std::cout << "src_lon, dis_src_lon, dis_src_err_lon : " << src_lon*RAD2DEG << ", " << dis_src_lon*RAD2DEG << ", " << dis_src_err_lon << std::endl;
            // std::cout << "src_lat, dis_src_lat, dis_src_err_lat : " << src_lat*RAD2DEG << ", " << dis_src_lat*RAD2DEG << ", " << dis_src_err_lat << std::endl;
            // std::cout << "src_r, dis_src_r, dis_src_err_r   : " << src_r   << ", " << dis_src_r   << ", " << dis_src_err_r   << std::endl;

            if (if_verbose){
                std::cout << "src positions lon lat r                  :    " << src_lon     << " " << src_lat     << " " << src_r     << std::endl;
                std::cout << "src positions lon(deg) lat(deg) depth(km):    " << src_lon*RAD2DEG << " " << src_lat*RAD2DEG << " " << radius2depth(src_r) << std::endl;
                std::cout << "src discretized position id i j k        :    " << i_src_loc       << " " << j_src_loc       << " " << k_src_loc    << std::endl;
                std::cout << "src discretized position lon lat r       :    " << dis_src_lon << " " << dis_src_lat << " " << dis_src_r << std::endl;
                std::cout << "src position error lon lat r             :    " << error_lon   << " " << error_lat   << " " << error_r   << std::endl;
                std::cout << "src relative position error lon lat r    :    " << dis_src_err_lon << " " << dis_src_err_lat << " " << dis_src_err_r << std::endl;
                std::cout << "delta lon lat r                          :    " << delta_lon   << " " << delta_lat   << " " << delta_r   << std::endl;
            }
            if(if_verbose) std::cout << "source is in the subdomain of rank " << myrank << std::endl;
            if(if_verbose) std::cout << "local index of the source: " << i_src_loc << " " << j_src_loc << " " << k_src_loc << std::endl;
        }

        if(if_verbose) stdout_by_main("--- end source initialization ---");
    }
}


Source::~Source() {
}


CUSTOMREAL Source::get_fac_at_source(CUSTOMREAL *loc_array) {

    // calculate factor by the rank where the source is located and broadcast to all

    CUSTOMREAL fac = 0.0;

    if (is_in_subdomain) {
        fac = (_1_CR - dis_src_err_r)*(_1_CR - dis_src_err_lat)*(_1_CR - dis_src_err_lon) * get_fac_at_point(loc_array,0,0,0) \
            + (_1_CR - dis_src_err_r)*(_1_CR - dis_src_err_lat)*         dis_src_err_lon  * get_fac_at_point(loc_array,1,0,0) \
            + (_1_CR - dis_src_err_r)*         dis_src_err_lat *(_1_CR - dis_src_err_lon) * get_fac_at_point(loc_array,0,1,0) \
            + (_1_CR - dis_src_err_r)*         dis_src_err_lat *         dis_src_err_lon  * get_fac_at_point(loc_array,1,1,0) \
            +          dis_src_err_r *(_1_CR - dis_src_err_lat)*(_1_CR - dis_src_err_lon) * get_fac_at_point(loc_array,0,0,1) \
            +          dis_src_err_r *(_1_CR - dis_src_err_lat)*         dis_src_err_lon  * get_fac_at_point(loc_array,1,0,1) \
            +          dis_src_err_r *         dis_src_err_lat *(_1_CR - dis_src_err_lon) * get_fac_at_point(loc_array,0,1,1) \
            +          dis_src_err_r *         dis_src_err_lat *         dis_src_err_lon  * get_fac_at_point(loc_array,1,1,1);
    } else {
        // do nothing
    }

    broadcast_cr_single(fac, src_rank);
    // std::cout << "interp: " << dis_src_err_r << ' ' << dis_src_err_lat << ' ' << dis_src_err_lon << ' ' << get_fac_at_point(loc_array,0,0,0) << std::endl;
    return fac;


}


CUSTOMREAL Source::get_fac_at_point(CUSTOMREAL* loc_array, int iplus, int jplus, int kplus){
    return loc_array[I2V(i_src_loc+iplus, j_src_loc+jplus, k_src_loc+kplus)];
}