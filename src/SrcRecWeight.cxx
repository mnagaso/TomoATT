#include <iostream>
#include <vector>

#include "utils.h"
#include "src_rec.h"
#include "mpi_funcs.h"


// copy calculated weights from swapped sources and receivers to the original sources and receivers
void reverse_rec_weight(std::vector<SrcRec>& src_points,      std::vector<std::vector<SrcRec>>& rec_points, \
                        std::vector<SrcRec>& src_points_back, std::vector<std::vector<SrcRec>>& rec_points_back) {

    // loop swapped sources
    for (long unsigned int i_src = 0; i_src < src_points.size(); i_src++){
        // swapped only the regional events && really need swap
        CUSTOMREAL weight     = src_points[i_src].weight;

        // copy calculated weight to the original rec_points_back
        for (auto& i_src_orig: src_points[i_src].id_srcs_ori){
            for (long unsigned int i_rec = 0; i_rec < rec_points_back[i_src_orig].size(); i_rec++){
                rec_points_back[i_src_orig][i_rec].weight = weight;
            }
        }
    }
}


/*
    caluclate geographical weight for sources and receivers.
    At first, we calculate inversed weight for each source from summation of epicentral distances to all sources.
    Then, we calculate inversed weight for each receiver from summation of epicentral distances to all receivers.
*/
void calculate_src_rec_weight(std::vector<SrcRec> &src_points, std::vector<std::vector<SrcRec>> &rec_points)
{

    CUSTOMREAL d_zero = 1.0; // initial value of reference sdistance parameter
    int n_src = src_points.size();

    // initialize all the source and receiver weights to be zero
    for (int i = 0; i < n_src; i++) {
        src_points[i].weight = 0.0;
        int n_rec = rec_points[i].size();
        for (int j = 0; j < n_rec; j++) {
            rec_points[i][j].weight = 0.0;
        }
    }

    // calculate source weights
    for (int i = 0; i < n_src-1; i++) {
        for (int j = i+1; j < n_src; j++) {
            CUSTOMREAL d_ij = 0.0;
            Epicentral_distance_sphere(src_points[i].lat*DEG2RAD, \
                                       src_points[i].lon*DEG2RAD, \
                                       src_points[j].lat*DEG2RAD, \
                                       src_points[j].lon*DEG2RAD, \
                                       d_ij);

            CUSTOMREAL w_inv_tmp = std::exp(-(d_ij/d_zero)*(d_ij/d_zero));

            src_points[i].weight += w_inv_tmp;
            src_points[j].weight += w_inv_tmp;
        }
    }

    // here the weight is inversed weight, so we need to invert it
    for (int i = 0; i < n_src; i++) {
        src_points[i].weight = 1.0/src_points[i].weight;
    }

    // before calculating receiver weights, we need to make a unique list of receivers
    // so do src-rec swap and do weight calculation for swapped sources

    // swap src-rec
    std::vector<SrcRec> src_points_back = src_points;
    std::vector<std::vector<SrcRec>> rec_points_back = rec_points;

    do_swap_src_rec(src_points, rec_points, src_points_back, rec_points_back);

    // calculate swapped source weights
    for (int i = 0; i < n_src-1; i++) {
        for (int j = i+1; j < n_src; j++) {
            CUSTOMREAL d_ij = 0.0;
            Epicentral_distance_sphere(src_points[i].lat*DEG2RAD, \
                                       src_points[i].lon*DEG2RAD, \
                                       src_points[j].lat*DEG2RAD, \
                                       src_points[j].lon*DEG2RAD, \
                                       d_ij);

            CUSTOMREAL w_inv_tmp = std::exp(-(d_ij/d_zero)*(d_ij/d_zero));

            src_points[i].weight += w_inv_tmp;
            src_points[j].weight += w_inv_tmp;
        }
    }

    // here the weight is inversed weight, so we need to invert it
    for (int i = 0; i < n_src; i++) {
        src_points[i].weight = 1.0/src_points[i].weight;
    }


    // reverse swapped src-rec
    reverse_rec_weight(src_points, rec_points, \
                       src_points_back, rec_points_back);

    //
    src_points = src_points_back;
    rec_points = rec_points_back;

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

    // write src_rec file with calculated weights
    writeout_src_rec_file(output_file, src_points, rec_points);

    stdout_by_main("------------------------------------------------------");
    stdout_by_main("end Src Rec weight calculation.");
    stdout_by_main(("output file: " + output_file).c_str());
    stdout_by_main("------------------------------------------------------");
    return 0;

}