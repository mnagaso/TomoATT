#include "inv_grid.h"


//
// utilities for inversion grid setup
//

#define ON_R 0
#define ON_T 1
#define ON_P 2

CUSTOMREAL get_max_inv(InputParams& IP, const int id_axis) {
    switch (id_axis) {
        case ON_R:
            return depth2radius(IP.get_max_dep_inv());
        case ON_T:
            return IP.get_max_lat_inv();
        case ON_P:
            return IP.get_max_lon_inv();
        default:
            // should not reach here. stop
            std::cout << "unknown id_axis" << std::endl;
            exit(1);
    }
}

CUSTOMREAL get_min_inv(InputParams& IP, const int id_axis) {
    switch (id_axis) {
        case ON_R:
            return depth2radius(IP.get_min_dep_inv());
        case ON_T:
            return IP.get_min_lat_inv();
        case ON_P:
            return IP.get_min_lon_inv();
        default:
            // should not reach here. stop
            std::cout << "unknown id_axis" << std::endl;
            exit(1);
    }
}

CUSTOMREAL get_max_inv_ani(InputParams& IP, const int id_axis) {
    switch (id_axis) {
        case ON_R:
            return depth2radius(IP.get_max_dep_inv_ani());
        case ON_T:
            return IP.get_max_lat_inv_ani();
        case ON_P:
            return IP.get_max_lon_inv_ani();
        default:
            // should not reach here. stop
            std::cout << "unknown id_axis" << std::endl;
            exit(1);
    }
}

CUSTOMREAL get_min_inv_ani(InputParams& IP, const int id_axis) {
    switch (id_axis) {
        case ON_R:
            return depth2radius(IP.get_min_dep_inv_ani());
        case ON_T:
            return IP.get_min_lat_inv_ani();
        case ON_P:
            return IP.get_min_lon_inv_ani();
        default:
            // should not reach here. stop
            std::cout << "unknown id_axis" << std::endl;
            exit(1);
    }
}


InvGrid1dBase inv_1d_grid_selector(InputParams& IP, const int grid_type, const int id_axis, const int n_points) {
    switch (grid_type) {
        case INV_GRID_REGULAR:
            if (id_axis == ON_R)
                return InvGrid1dRegular(IP, n_points, get_min_inv(IP, id_axis), get_max_inv(IP, id_axis)); // regular grid for radius
            else
                return InvGrid1dRegular(IP, n_points, get_min_inv(IP, id_axis), get_max_inv(IP, id_axis), n_inv_K_loc); // regular grid for lat and lon
        case INV_GRID_FLEX:
            if (id_axis == ON_R)
                return InvGrid1dFlexible(IP, n_points, IP.get_dep_inv(), id_axis); // flexible grid for radius
            else if (id_axis == ON_T)
                return InvGrid1dFlexible(IP, n_points, IP.get_lat_inv(), id_axis, n_inv_K_loc); // flexible grid for lat
            else
                return InvGrid1dFlexible(IP, n_points, IP.get_lon_inv(), id_axis, n_inv_K_loc); // flexible grid for lon
        case INV_GRID_TRAPE:
            if (id_axis == ON_R) { // trapezoidal grid cannot be defined on r axis
                std::cout << "trapezoidal grid cannot be defined on r axis" << std::endl;
                exit(1);
            } else if (id_axis == ON_T) {
                CUSTOMREAL domain_size = (IP.get_min_lat() + IP.get_max_lat()) / _2_CR; // RAD
                return InvGrid1dTrapezoidal(IP, n_points, domain_size, IP.get_lat_spacing_inv(), id_axis, n_inv_K_loc); // trapezoidal grid for lat
            } else if (id_axis == ON_P) {
                CUSTOMREAL domain_size = (IP.get_min_lon() + IP.get_max_lon()) / _2_CR; // RAD
                return InvGrid1dTrapezoidal(IP, n_points, domain_size, IP.get_lon_spacing_inv(), id_axis, n_inv_K_loc); // trapezoidal grid for lon
            } else {
                std::cout << "unknown id_axis" << std::endl;
                exit(1);
            }

        default:
            // should not reach here. stop
            std::cout << "unknown type of inversion grid" << std::endl;
            exit(1);
    }
}


InvGrid1dBase inv_1d_grid_selector_ani(InputParams& IP, const int grid_type, const int id_axis, const int n_points) {
     switch (grid_type) {
        case INV_GRID_REGULAR:
            if (id_axis == ON_R)
                return InvGrid1dRegular(IP, n_points, get_min_inv_ani(IP, id_axis), get_max_inv_ani(IP, id_axis)); // regular grid for radius
            else
                return InvGrid1dRegular(IP, n_points, get_min_inv_ani(IP, id_axis), get_max_inv_ani(IP, id_axis), n_inv_K_loc_ani); // regular grid for lat and lon
        case INV_GRID_FLEX:
            if (id_axis == ON_R)
                return InvGrid1dFlexible(IP, n_points, IP.get_dep_inv_ani(), id_axis); // flexible grid for radius
            else if (id_axis == ON_T)
                return InvGrid1dFlexible(IP, n_points, IP.get_lat_inv_ani(), id_axis, n_inv_K_loc_ani); // flexible grid for lat
            else
                return InvGrid1dFlexible(IP, n_points, IP.get_lon_inv_ani(), id_axis, n_inv_K_loc_ani); // flexible grid for lon
        case INV_GRID_TRAPE:
            if (id_axis == ON_R) { // trapezoidal grid cannot be defined on r axis
                std::cout << "trapezoidal grid cannot be defined on r axis" << std::endl;
                exit(1);
            } else if (id_axis == ON_T) {
                CUSTOMREAL domain_size = (IP.get_min_lat() + IP.get_max_lat()) / _2_CR; // RAD
                return InvGrid1dTrapezoidal(IP, n_points, domain_size, IP.get_lat_spacing_inv_ani(), id_axis, n_inv_K_loc_ani); // trapezoidal grid for lat
            } else if (id_axis == ON_P) {
                CUSTOMREAL domain_size = (IP.get_min_lon() + IP.get_max_lon()) / _2_CR; // RAD
                return InvGrid1dTrapezoidal(IP, n_points, domain_size, IP.get_lon_spacing_inv_ani(), id_axis, n_inv_K_loc_ani); // trapezoidal grid for lon
            } else {
                std::cout << "unknown id_axis" << std::endl;
                exit(1);
            }

        default:
            // should not reach here. stop
            std::cout << "unknown type of inversion grid" << std::endl;
            exit(1);
    }
}


//
// Inversion grid class definition
//


InvGrid::InvGrid(InputParams& IP) {

    // read the parameters and store them in this class
    get_inv_grid_params(IP);

    // create 1D inversion grid for each dimension.
    // t and p grids are not 1D but 2D grids for a trapezoidal definition of the grid.

    // r (radius) grid
    r = inv_1d_grid_selector(IP, IP.get_type_invgrid_dep(), ON_R, n_inv_K_loc);
    // t (latitude) grid
    t = inv_1d_grid_selector(IP, IP.get_type_invgrid_lat(), ON_T, n_inv_J_loc);
    // p (longitude) grid
    p = inv_1d_grid_selector(IP, IP.get_type_invgrid_lon(), ON_P, n_inv_I_loc);

    // inversion grid for anisotropic parameters
    if (IP.get_invgrid_ani()){
        // r (radius) grid
        r_ani = inv_1d_grid_selector_ani(IP, IP.get_type_invgrid_dep_ani(), ON_R, n_inv_K_loc_ani);
        // t (latitude) grid
        t_ani = inv_1d_grid_selector_ani(IP, IP.get_type_invgrid_lat_ani(), ON_T, n_inv_J_loc_ani);
        // p (longitude) grid
        p_ani = inv_1d_grid_selector_ani(IP, IP.get_type_invgrid_lon_ani(), ON_P, n_inv_I_loc_ani);
    } else {
        // use the same inversion grid with slowessw
        r_ani = r;
        t_ani = t;
        p_ani = p;
    }
}


InvGrid::~InvGrid() {
}


void InvGrid::get_inv_grid_params(InputParams& IP) {

    n_inv_grids = IP.get_n_inversion_grid();

    // inversion grid for velocity
    if(IP.get_type_invgrid_dep() == 0){
        n_inv_K_loc = IP.get_n_inv_r();
    } else if (IP.get_type_invgrid_dep() == 1) {
        n_inv_K_loc = IP.get_n_inv_r_flex();
    } else {
        std::cout << "unknown type of inversion grid" << std::endl;
        exit(1);
    }

    if(IP.get_type_invgrid_lat() == 0){
        n_inv_J_loc = IP.get_n_inv_t();
    } else if (IP.get_type_invgrid_lat() == 1) {
        n_inv_J_loc = IP.get_n_inv_t_flex();
    } else if (IP.get_type_invgrid_lat() == 2) {
        n_inv_J_loc = IP.get_n_inv_t_trape();
    } else {
        std::cout << "unknown type of inversion grid" << std::endl;
        exit(1);
    }

    if(IP.get_type_invgrid_lon() == 0){
        n_inv_I_loc = IP.get_n_inv_p();
    } else if (IP.get_type_invgrid_lon() == 1) {
        n_inv_I_loc = IP.get_n_inv_p_flex();
    } else if (IP.get_type_invgrid_lon() == 2) {
        n_inv_I_loc = IP.get_n_inv_p_trape();
    } else {
        std::cout << "unknown type of inversion grid" << std::endl;
        exit(1);
    }

    // inversion grid for anisotropy (requires flex inversion grid setup)
    if(IP.get_invgrid_ani()){
        if (IP.get_type_invgrid_dep_ani() == 0)
            n_inv_K_loc_ani =IP.get_n_inv_r_ani();
        else if (IP.get_type_invgrid_dep_ani() == 1)
            n_inv_K_loc_ani = IP.get_n_inv_r_flex_ani();
        else {
            std::cout << "unknown type of inversion grid" << std::endl;
            exit(1);
        }
        if (IP.get_type_invgrid_lat_ani() == 0)
            n_inv_J_loc_ani = IP.get_n_inv_t_ani();
        else if (IP.get_type_invgrid_lat_ani() == 1)
            n_inv_J_loc_ani = IP.get_n_inv_t_flex_ani();
        else if (IP.get_type_invgrid_lat_ani() == 2)
            n_inv_J_loc_ani = IP.get_n_inv_t_trape_ani();
        else {
            std::cout << "unknown type of inversion grid" << std::endl;
            exit(1);
        }
        if (IP.get_type_invgrid_lon_ani() == 0)
            n_inv_I_loc_ani = IP.get_n_inv_p_ani();
        else if (IP.get_type_invgrid_lon_ani() == 1)
            n_inv_I_loc_ani = IP.get_n_inv_p_flex_ani();
        else if (IP.get_type_invgrid_lon_ani() == 2)
            n_inv_I_loc_ani = IP.get_n_inv_p_trape_ani();
        else {
            std::cout << "unknown type of inversion grid" << std::endl;
            exit(1);
        }
    } else {
        // use the same inversion grid with slowessw
        n_inv_I_loc_ani = n_inv_I_loc;
        n_inv_J_loc_ani = n_inv_J_loc;
        n_inv_K_loc_ani = n_inv_K_loc;
    }

}


void InvGrid::write_inversion_grid_to_file(){
    std::ofstream ofs;

    std::string inversion_grid_file_out = output_dir + "/inversion_grid.txt";
    ofs.open(inversion_grid_file_out);

    if(subdom_main && id_subdomain == 0){       // main processor of subdomain && the first id of subdoumains
        // inversion grid of velocity
        for(int l = 0; l < n_inv_grids; l++){
            ofs << l << " " << n_inv_K_loc << " " << n_inv_J_loc << " " << n_inv_I_loc << std::endl;    // number of ivnersion grid

            for(int k = 0; k < n_inv_K_loc; k++){
                ofs << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ')
                    << radius2depth(r.arr[I2V_INV_GRIDS_1DK(k,l)]) << " ";
                ofs << std::endl;
                for(int j =0; j < n_inv_J_loc; j++)
                    ofs << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ')
                        << t.arr[I2V_INV_GRIDS_2DJ(j,k,l)]*RAD2DEG << " ";
                ofs << std::endl;
                for(int i =0; i < n_inv_I_loc; i++)
                    ofs << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ')
                        << p.arr[I2V_INV_GRIDS_2DI(i,k,l)]*RAD2DEG << " ";
                ofs << std::endl;
            }
        }
        // inversion grid of anisotropy
        for(int l = 0; l < n_inv_grids; l++){
            ofs << l << " " << n_inv_K_loc_ani << " " << n_inv_J_loc_ani << " " << n_inv_I_loc_ani << std::endl;    // number of ivnersion grid
            for(int k = 0; k < n_inv_K_loc_ani; k++){
                ofs << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ')
                    << radius2depth(r_ani.arr[I2V_INV_ANI_GRIDS_1DK(k,l)]) << " ";
                ofs << std::endl;
                for(int j =0; j < n_inv_J_loc_ani; j++)
                    ofs << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ')
                        << t_ani.arr[I2V_INV_ANI_GRIDS_1DJ(j,k,l)]*RAD2DEG << " ";
                ofs << std::endl;
                for(int i =0; i < n_inv_I_loc_ani; i++)
                    ofs << std::fixed << std::setprecision(4) << std::setw(9) << std::right << std::setfill(' ')
                        << p_ani.arr[I2V_INV_ANI_GRIDS_1DI(i,k,l)]*RAD2DEG << " ";
                ofs << std::endl;
            }
        }
    }
}


// regular grid positions for radius
InvGrid1dRegular::InvGrid1dRegular(InputParams& IP, const int n_inv_loc,
                                  const CUSTOMREAL min_inv, const CUSTOMREAL max_inv) {
    // set n
    n = n_inv_loc*n_inv_grids;

    // allocate memory for the grid
    arr = new CUSTOMREAL[n];

    // inversion grid is defined for all processes which covers entire domain
    dinv = (max_inv - min_inv) / (n_inv_loc-2);
    // shift of each set of inversion grid
    dinv_l = dinv/n_inv_grids;

    for (int l = 0; l < n_inv_grids; l++) {
        for (int i = 0; i < n_inv_loc; i++)
            arr[I2V_INV_GRIDS_1D_GENERIC(i, l, n_inv_loc)] = min_inv + i*dinv - l*dinv_l;
            //arr[I2V_INV_GRIDS_1D_GENERIC(n_inv_loc-1-i, l, n_inv_loc)] = min_inv + i*dinv - l*dinv_l;
    }
}


// regular grid positions for lat and lon
InvGrid1dRegular::InvGrid1dRegular(InputParams& IP, const int n_inv_loc,
                                  const CUSTOMREAL min_inv, const CUSTOMREAL max_inv,
                                  const int n_inv_loc_k) {

    // set n
    n = n_inv_loc*n_inv_loc_k*n_inv_grids;

    // allocate memory for the grid
    arr = new CUSTOMREAL[n];

    // inversion grid is defined for all processes which covers entire domain
    dinv = (max_inv - min_inv) / (n_inv_loc-2);
    // shift of each set of inversion grid
    dinv_l = dinv/n_inv_grids;

    for (int l = 0; l < n_inv_grids; l++) {
        for (int k = 0; k < n_inv_loc_k; k++) {
            for (int i = 0; i < n_inv_loc; i++)
                arr[I2V_INV_GRIDS_2D_GENERIC(i,k,l,n_inv_loc,n_inv_loc_k)] = min_inv + i*dinv - l*dinv_l;
        }
    }

}


// flexible grid positions for radius
InvGrid1dFlexible::InvGrid1dFlexible(InputParams& IP, const int n_inv_loc,
                                  const CUSTOMREAL* input_pos, const int axis_id) {

    // set n
    n = n_inv_loc*n_inv_grids;

    // allocate memory for the grid
    arr = new CUSTOMREAL[n];

    // shift of each set of inversion grid
    dinv_l = dinv/n_inv_grids;

    // create a copy array of input_pos
    CUSTOMREAL* arr_input_pos = new CUSTOMREAL[n_inv_loc];

    // change the unit depending on the axis
    if (axis_id == ON_R){
        for (int i = 0; i < n_inv_loc; i++)
            arr_input_pos[i] = depth2radius(input_pos[i]);
    } else {
        for (int i = 0; i < n_inv_loc; i++)
            arr_input_pos[i] = input_pos[i]*DEG2RAD;
    }

    for (int i = 0; i < n_inv_loc; i++){
        if (i < n_inv_loc-1)
            dinv_l = (arr_input_pos[i+1] - arr_input_pos[i]) / n_inv_grids;
        else
            dinv_l = (arr_input_pos[n_inv_loc-1] - arr_input_pos[n_inv_loc-2]) / n_inv_grids;

        for (int l = 0; l < n_inv_grids; l++)
            arr[I2V_INV_GRIDS_1D_GENERIC(i, l, n_inv_loc)] = arr_input_pos[i] - l*dinv_l;

    }

    delete[] arr_input_pos;
}


// flexible grid positions for lat and lon
InvGrid1dFlexible::InvGrid1dFlexible(InputParams& IP, const int n_inv_loc,
                                  const CUSTOMREAL* input_pos, const int axis_id,
                                  const int n_inv_loc_k) {

    // set n
    n = n_inv_loc*n_inv_loc_k*n_inv_grids;

    // allocate memory for the grid
    arr = new CUSTOMREAL[n];

    // shift of each set of inversion grid
    dinv_l = dinv/n_inv_grids;

    // create a copy array of input_pos
    CUSTOMREAL* arr_input_pos = new CUSTOMREAL[n_inv_loc];

    // change the unit depending on the axis
    if (axis_id == ON_R){
        for (int i = 0; i < n_inv_loc; i++)
            arr_input_pos[i] = depth2radius(input_pos[i]);
    } else {
        for (int i = 0; i < n_inv_loc; i++)
            arr_input_pos[i] = input_pos[i]*DEG2RAD;
    }

    for (int i = 0; i < n_inv_loc; i++){
        if (i < n_inv_loc-1)
            dinv_l = (arr_input_pos[i+1] - arr_input_pos[i]) / n_inv_grids;
        else
            dinv_l = (arr_input_pos[n_inv_loc-1] - arr_input_pos[n_inv_loc-2]) / n_inv_grids;

        for (int l = 0; l < n_inv_grids; l++) {
            for (int k = 0; k < n_inv_loc_k; k++)
                arr[I2V_INV_GRIDS_2D_GENERIC(i,k,l,n_inv_loc,n_inv_loc_k)] = arr_input_pos[i] - l*dinv_l;
        }
    }

    delete[] arr_input_pos;
}


// trapezoidal grid positions for lat and lon
InvGrid1dTrapezoidal::InvGrid1dTrapezoidal(InputParams& IP, const int n_inv_loc,
                                  const CUSTOMREAL domain_size,
                                  const CUSTOMREAL* input_spacing, const int axis_id,
                                  const int n_inv_loc_k) {

    // set n
    n = n_inv_loc*n_inv_loc_k*n_inv_grids;

    // allocate memory for the grid
    arr = new CUSTOMREAL[n];

    // shift of each set of inversion grid
    dinv_l = dinv/n_inv_grids;

    // create a copy array of input_pos
    CUSTOMREAL* arr_input_spacing = new CUSTOMREAL[n_inv_loc_k];

    // change the unit depending on the axis
    for (int i = 0; i < n_inv_loc_k; i++)
        arr_input_spacing[i] = input_spacing[i]*DEG2RAD;

    for (int l = 0; l < n_inv_grids; l++) {
        for (int i = 0; i < n_inv_loc; i++)
            for(int k = 0; k < n_inv_loc_k; k++){
                dinv   = arr_input_spacing[k]; // already in RAD
                dinv_l = dinv/n_inv_grids;     // already in RAD
                arr[I2V_INV_GRIDS_2D_GENERIC(i,k,l,n_inv_loc,n_inv_loc_k)] = domain_size + (i-n_inv_loc/_2_CR)*dinv - l*dinv_l;
            }
    }

    delete[] arr_input_spacing;
}




