#include <iostream>
#include <vector>

#include "utils.h"
#include "src_rec.h"
#include "mpi_funcs.h"



/*
    caluclate geographical weight for sources and receivers.
    At first, we calculate inversed weight for each source from summation of epicentral distances to all sources.
    Then, we calculate inversed weight for each receiver from summation of epicentral distances to all receivers.

    ref: Youyi Ruan(2019) https://academic.oup.com/gji/article/219/2/1225/5542709
*/

void calc_dist_min_max(std::vector<SrcRec>& src_or_rec, CUSTOMREAL& d_min, CUSTOMREAL& d_max){
    // calculate the min and max epicentral distance of all src_or_rec

    // check the min and max epicentral distance of all src_or_rec
    CUSTOMREAL tmp_min = 1.0e+10;
    CUSTOMREAL tmp_max = 0.0;

    int n_elm = src_or_rec.size();

    for (int i = 0; i < n_elm-1; i++) {
        for (int j = i+1; j < n_elm; j++) {
            CUSTOMREAL d_ij = 0.0;
            Epicentral_distance_sphere(src_or_rec[i].lat*DEG2RAD, \
                                       src_or_rec[i].lon*DEG2RAD, \
                                       src_or_rec[j].lat*DEG2RAD, \
                                       src_or_rec[j].lon*DEG2RAD, \
                                       d_ij);

            if (d_ij < tmp_min) tmp_min = d_ij;
            if (d_ij > tmp_max) tmp_max = d_ij;
        }
    }

    // set result
    d_min = tmp_min;
    d_max = tmp_max;

    // print the min and max epicentral distance
    //std::cout << "min epicentral distance = " << d_min << std::endl;
    //std::cout << "max epicentral distance = " << d_max << std::endl;

}

void init_weight(std::vector<SrcRec>& src_or_rec){
    // initialize all the source and receiver weights to be zero
    int n_elm = src_or_rec.size();
    for (int i = 0; i < n_elm; i++) {
        src_or_rec[i].weight = 0.0;
    }
}


CUSTOMREAL _calc_weight(std::vector<SrcRec>& src_or_rec, CUSTOMREAL& d_zero){
    int n_elm = src_or_rec.size();
    for (int i = 0; i < n_elm-1; i++) {
        for (int j = i+1; j < n_elm; j++) {
            CUSTOMREAL d_ij = 0.0;
            Epicentral_distance_sphere(src_or_rec[i].lat*DEG2RAD, \
                                       src_or_rec[i].lon*DEG2RAD, \
                                       src_or_rec[j].lat*DEG2RAD, \
                                       src_or_rec[j].lon*DEG2RAD, \
                                       d_ij);

            d_ij *= RAD2DEG;

            CUSTOMREAL w_inv_tmp = std::exp(-(d_ij/d_zero)*(d_ij/d_zero));

            src_or_rec[i].weight += w_inv_tmp;
            src_or_rec[j].weight += w_inv_tmp;
        }
    }

    // here the weight is inversed weight, so we need to invert it
    for (int i = 0; i < n_elm; i++) {
        src_or_rec[i].weight = 1.0/src_or_rec[i].weight;
    }

    // calculate condition number (min/max of weight)
    CUSTOMREAL w_min = 1.0e+10;
    CUSTOMREAL w_max = 0.0;
    for (int i = 0; i < n_elm; i++) {
        if (src_or_rec[i].weight < w_min) w_min = src_or_rec[i].weight;
        if (src_or_rec[i].weight > w_max) w_max = src_or_rec[i].weight;
    }

    std::cout << "weight min = " << w_min << std::endl;
    std::cout << "weight max = " << w_max << std::endl;
    std::cout << "condition number = " << w_max/w_min << std::endl;

    return w_max/w_min;

}


int find_good_ncond(std::vector<CUSTOMREAL>& condition_numbers){
    // select the best condition number

    int n_try = condition_numbers.size();

    CUSTOMREAL ncond_min = 1.0e+10;
    CUSTOMREAL ncond_max = 0.0;

    for (int i = 0; i < n_try; i++) {
        // ignore abnormal value
        if (std::isnan(condition_numbers[i])) continue;
        if (std::isinf(condition_numbers[i])) continue;

        if (condition_numbers[i] < ncond_min) {
            ncond_min = condition_numbers[i];
        }
        if (condition_numbers[i] > ncond_max) {
            ncond_max = condition_numbers[i];
        }
    }

    // find the first one-third point between min and max condition number
    CUSTOMREAL ncond_1_3 = ncond_min + (ncond_max - ncond_min) / 3.0;
    int ncond_1_3_idx = 0;

    for (int i = 0; i < n_try; i++) {
        // ignore abnormal value
        if (std::isnan(condition_numbers[i])) continue;
        if (std::isinf(condition_numbers[i])) continue;

        if (condition_numbers[i] > ncond_1_3) {
            ncond_1_3_idx = i;
            break;
        }
    }

    return ncond_1_3_idx;
}


void normalize_weight(std::vector<SrcRec>& src_or_rec){

    // scale weight values to be the total sum becomes 1.0*number of elements
    int n_elm = src_or_rec.size();
    CUSTOMREAL w_sum = 0.0;
    for (int i = 0; i < n_elm; i++) {
        w_sum += src_or_rec[i].weight;
    }

    for (int i = 0; i < n_elm; i++) {
        src_or_rec[i].weight = (src_or_rec[i].weight / w_sum) * n_elm;
    }

//    // scale the weight values to be within a range of 0.1 ~ 10.0
//    CUSTOMREAL vmax = 10.0;
//    CUSTOMREAL vmin = 0.1;
//
//    int n_elm = src_or_rec.size();
//
//    CUSTOMREAL w_min = 1.0e+10;
//    CUSTOMREAL w_max = 0.0;
//
//    for (int i = 0; i < n_elm; i++) {
//        if (src_or_rec[i].weight < w_min) w_min = src_or_rec[i].weight;
//        if (src_or_rec[i].weight > w_max) w_max = src_or_rec[i].weight;
//    }
//
//    CUSTOMREAL w_scale = vmax / (w_max - w_min);
//
//    for (int i = 0; i < n_elm; i++) {
//        src_or_rec[i].weight = (src_or_rec[i].weight - w_min) * w_scale + vmin;
//    }

}


void calc_weight(std::vector<SrcRec>& src_or_rec){
    // calculate the weights of sources or receivers
//    int n_try = 30; // number of tries

    // get the min and max epicentral distance of all src_or_rec
//    CUSTOMREAL d_min, d_max;
//    calc_dist_min_max(src_or_rec, d_min, d_max);

    // convert rad to deg
    //d_min *= RAD2DEG;
    //d_max *= RAD2DEG;

//    d_max = 5.0; // in deg
//    d_min = 0.01; // shouldn't be 0.0
//
//    std::cout << "ref distance min = " << d_min << std::endl;
//    std::cout << "ref distance max = " << d_max << std::endl;
//
//    std::vector<CUSTOMREAL> condition_numbers; // store the condition numbers
//    std::vector<CUSTOMREAL> test_d_zero; // store the test d_zero
//
//    // prepare the test d_zero (d_min ~ d_max by n_try times)
//    CUSTOMREAL d_zero_step = (d_max - d_min) / (n_try - 1);
//    CUSTOMREAL d_zero_tmp = d_min; // void d_zero = 0.0
//    for (int i = 0; i < n_try; i++) {
//        test_d_zero.push_back(d_zero_tmp);
//        d_zero_tmp += d_zero_step;
//    }
//
    // calculate the weights iteratively
//    int c = 0;
//    for (auto& d_zero_try : test_d_zero) {
//        // initialize all the source or receiver weights to be zero
//        init_weight(src_or_rec);
//        CUSTOMREAL ncond = _calc_weight(src_or_rec, d_zero_try);
//
//        // store the condition number
//        condition_numbers.push_back(ncond);
//
//        // output d_zero and condition number
//        std::cout << "i-iter/total: " << c << "/" << n_try << ", d_zero: " << d_zero_try << ", condition number: " << ncond << std::endl;
//
//        c++;
//    }
//
//    // find the one-third point between min and max condition number
//    CUSTOMREAL d_zero_fin = test_d_zero[find_good_ncond(condition_numbers)];

    // use largest distance as d_zero_fin
    CUSTOMREAL d_zero_fin = ref_value;

    // calculate the weights with the final d_zero
    init_weight(src_or_rec);
    _calc_weight(src_or_rec, d_zero_fin);

    // normalize the weights
    normalize_weight(src_or_rec);

}

void calculate_src_rec_weight(std::vector<SrcRec> &src_points, std::vector<std::vector<SrcRec>> &rec_points)
{

    int n_src = src_points.size();

    // make a unique list of receivers
    std::vector<SrcRec> rec_points_unique;
    for (int i = 0; i < n_src; i++) {
        int n_rec = rec_points[i].size();
        for (int j = 0; j < n_rec; j++) {
            bool is_unique = true;
            for (auto& rec: rec_points_unique) {
                // check if the receiver is the same as the new source
                if (rec_points[i][j].name_rec.compare(rec.name_rec)==0) {
                    is_unique = false;
                    break;
                }
            }
            if (is_unique) {
                rec_points_unique.push_back(rec_points[i][j]);
            }
        }
    }

    // calculate source weights
    std::cout << "calculating source weights..." << std::endl;
    calc_weight(src_points);

    // calculate receiver weights
    std::cout << "calculating receiver weights..." << std::endl;
    calc_weight(rec_points_unique);

    // assign the receiver weights to the corresponding receivers
    for (int i = 0; i < n_src; i++) {
        int n_rec = rec_points[i].size();
        for (int j = 0; j < n_rec; j++) {

            // search for the receiver in the unique list
            for (auto& rec: rec_points_unique) {
                // check if the receiver is the same as the new source
                if (rec_points[i][j].name_rec.compare(rec.name_rec)==0) {
                    rec_points[i][j].weight = rec.weight;
                    break;
                }
            }
        }
    }

    // end
}


void copy_arrival_times(std::vector<std::vector<SrcRec>>& rec_points){
    // copy arrival times read from input to output
    int n_src = rec_points.size();
    for (int i = 0; i < n_src; i++) {
        int n_rec = rec_points[i].size();
        for (int j = 0; j < n_rec; j++) {
            rec_points[i][j].arr_time = rec_points[i][j].arr_time_ori;
        }
    }
}


// function for calculating the source and receiver weight
int main(int argc, char *argv[])
{
    // parse options
    parse_options_srcrec_weight(argc, argv);

    // init mpi
    initialize_mpi();

    stdout_by_main("------------------------------------------------------");
    stdout_by_main("start Src Rec weight calculation.");
    stdout_by_main("------------------------------------------------------");

    std::vector<SrcRec> src_points;
    std::vector<std::vector<SrcRec>> rec_points;
    std::string output_file = "src_rec_with_weight.txt";

    // read src_rec file
    parse_src_rec_file(input_file, src_points, rec_points);

    // calculate source and receiver weight
    calculate_src_rec_weight(src_points, rec_points);

    // copy arrival times read from input to output
    copy_arrival_times(rec_points);

    // write src_rec file with calculated weights
    writeout_src_rec_file(output_file, src_points, rec_points);

    stdout_by_main("------------------------------------------------------");
    stdout_by_main("end Src Rec weight calculation.");
    stdout_by_main(("output file: " + output_file).c_str());
    stdout_by_main("------------------------------------------------------");
    return 0;

}