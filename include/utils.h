#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream>
#include <sys/stat.h>
#include <iomanip>
#include <sstream>

#include "config.h"


inline void create_output_dir(std::string dir_path){
    // create output directory
    if (mkdir(dir_path.c_str(), 0777) == -1){
        std::cout << "Warning : directory " << dir_path << " can not be created. Maybe already exists (no problem in this case)." << std::endl;
    }
}


inline bool is_file_exist(const char* fileName)
{
    return static_cast<bool>(std::ifstream(fileName));
}


inline void stdout_by_main(char const* str){
    if (sim_rank == 0 && inter_sub_rank == 0 && sub_rank == 0)
        std::cout << str << std::endl;
}


inline void parse_options(int argc, char* argv[]){
    bool input_file_found = false;

    for (int i = 1; i < argc; i++){
        if(strcmp(argv[i], "-v") == 0)
            if_verbose = true;
        else if (strcmp(argv[i],"-i") == 0){
            input_file = argv[i+1];
            input_file_found = true;
        }
    }

    // error if input_file is  not found
    if(!input_file_found){
        stdout_by_main("usage: mpirun -np 4 ./TOMOATT -i input_params.yaml");
        std::cout << argc << std::endl;
        exit(EXIT_FAILURE);
    }
}


template <typename T>
inline int check_data_type(T const& data){
    if (std::is_same<T,bool>::value){
        return 0;
    } else if (std::is_same<T, int>::value){
        return 1;
    } else if (std::is_same<T, float>::value){
        return 2;
    } else if (std::is_same<T, double>::value){
        return 3;
    } else {
        std::cout << "Error: custom real type is not float or double" << std::endl;
        std::cout << "Please check the definition of CUSTOMREAL in io.h" << std::endl;
        std::cout << "Exiting..." << std::endl;
        exit(1);
    }
}


template <typename T>
inline T my_square(T const& a){
    return a * a;
}


// defined function is more than 2 times slower than inline function
//#define my_square(a) (a * a)


inline void RLonLat2xyz(CUSTOMREAL lon, CUSTOMREAL lat, CUSTOMREAL R, CUSTOMREAL& x, CUSTOMREAL& y, CUSTOMREAL& z){
    /*
    lon : longitude in radian
    lat : latitude in radian
    R : radius

    x : x coordinate in m
    y : y coordinate in m
    z : z coordinate in m
    */
    //x = R*cos(lat*DEG2RAD)*cos(lon*DEG2RAD);
    //y = R*cos(lat*DEG2RAD)*sin(lon*DEG2RAD);
    //z = R*sin(lat*DEG2RAD);
    x = R*cos(lat)*cos(lon);
    y = R*cos(lat)*sin(lon);
    z = R*sin(lat);

}


// calculate epicentral distance in radian from lon lat in radian
inline void Epicentral_distance_sphere(CUSTOMREAL lat1, CUSTOMREAL lon1, \
                                       CUSTOMREAL lat2, CUSTOMREAL lon2, \
                                       CUSTOMREAL& dist) {
    if (isZero(my_square((lat1-lat2)) \
     &&      + my_square((lon1-lon2)))){
        dist = _0_CR;
    } else {
        dist = std::abs(acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon2-lon1)));
    }
}


//calculate azimuth
inline void Azimuth_sphere(CUSTOMREAL lat1, CUSTOMREAL lon1, \
                             CUSTOMREAL lat2, CUSTOMREAL lon2, \
                             CUSTOMREAL& azi) {

    if (isZero(my_square((lat1-lat2)) \
        &&   + my_square((lon1-lon2)))){
        azi = _0_CR;
    } else {
        azi = atan2(sin(lon2-lon1)*cos(lat2), \
                    cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(lon2-lon1));
        if (azi < _0_CR)
            azi += _2_CR * PI;
    }
}


inline void WGS84ToCartesian(CUSTOMREAL& lon, CUSTOMREAL& lat, CUSTOMREAL R, \
                             CUSTOMREAL& x, CUSTOMREAL& y, CUSTOMREAL& z, \
                             CUSTOMREAL& lon_center, CUSTOMREAL& lat_center){

    // equatorial radius WGS84 major axis
    const static CUSTOMREAL equRadius = R_earth; // in m in this function
    const static CUSTOMREAL flattening = 1.0 / 298.257222101;

    const static CUSTOMREAL sqrEccentricity = flattening * (2.0 - flattening);

    const CUSTOMREAL lat_rad = lat - lat_center; // input should be already in radian
    const CUSTOMREAL lon_rad = lon - lon_center;
    const CUSTOMREAL alt = R - equRadius; // altitude (height above sea level)

    const CUSTOMREAL sinLat = sin(lat_rad);
    const CUSTOMREAL cosLat = cos(lat_rad);
    const CUSTOMREAL sinLon = sin(lon_rad);
    const CUSTOMREAL cosLon = cos(lon_rad);

    // Normalized radius
    const CUSTOMREAL normRadius = equRadius / sqrt(1.0 - sqrEccentricity * sinLat * sinLat);

    x = (normRadius + alt) * cosLat * cosLon;
    y = (normRadius + alt) * cosLat * sinLon;
    z = (normRadius * (1.0 - sqrEccentricity) + alt) * sinLat;

}


template <typename T>
inline T dot_product(T const* const a, T const* const b, int const& n){
    CUSTOMREAL result = _0_CR;
    for (int i = 0; i < n; i++){
        result += a[i] * b[i];
    }
    return result;
}

template <typename T>
inline T find_max(T const* const a, int const& n){
    T max = a[0];
    for (int i = 1; i < n; i++){
        if (a[i] > max){
            max = a[i];
        }
    }
    return max;
}

template <typename T>
inline T find_absmax(T const* const a, int const& n){
    T max = fabs(a[0]);
    for (int i = 1; i < n; i++){
        if (fabs(a[i]) > max){
            max = fabs(a[i]);
        }
    }
    return max;
}

template <typename T>
inline T calc_l2norm(T const* const a, int const& n){
    T result = _0_CR;
    for (int i = 0; i < n; i++){
        result += a[i] * a[i];
    }
    return result;
}


inline std::string int2string_zero_fill(int i) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(4) << i;
    return ss.str();
}



#endif // UTILS_H