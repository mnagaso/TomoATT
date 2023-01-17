#include <iostream>
#include <vector>

#include "utils.h"
#include "src_rec.h"
#include "mpi_funcs.h"

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

    // calculate receiver weight
    for (int isrc = 0; isrc < n_src; isrc++){
        int n_rec = rec_points[isrc].size();
        for (int irec = 0; irec < n_rec-1; irec++){
            for (int jrec = irec+1; jrec < n_rec; jrec++){
                CUSTOMREAL d_ij = 0.0;
                Epicentral_distance_sphere(rec_points[isrc][irec].lat*DEG2RAD, \
                                           rec_points[isrc][irec].lon*DEG2RAD, \
                                           rec_points[isrc][jrec].lat*DEG2RAD, \
                                           rec_points[isrc][jrec].lon*DEG2RAD, \
                                           d_ij);

                CUSTOMREAL w_inv_tmp = std::exp(-(d_ij/d_zero)*(d_ij/d_zero));

                rec_points[isrc][irec].weight += w_inv_tmp;
                rec_points[isrc][jrec].weight += w_inv_tmp;
            }
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

    // write src_rec file with calculated weights
    writeout_src_rec_file(output_file, src_points, rec_points);

    stdout_by_main("------------------------------------------------------");
    stdout_by_main("end Src Rec weight calculation.");
    stdout_by_main(("output file: " + output_file).c_str());
    stdout_by_main("------------------------------------------------------");
    return 0;

}