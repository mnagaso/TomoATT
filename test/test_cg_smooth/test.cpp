// unit test code for ../include/smooth.h

#include <iostream>
#include "smooth_test.h"
#include "../../include/config.h"


void out_file(CUSTOMREAL* arr, std::string fname){
    // output initial and final grid in file
    std::ofstream out_file;
    out_file.open(fname);
    for (int k=0; k<loc_K; k++){
        for (int j=0; j<loc_J; j++){
            for (int i=0; i<loc_I; i++){
                out_file << arr[i+j*loc_I+k*loc_I*loc_J] << "   ";
            }
            out_file << std::endl;
        }
        //out_file << std::endl;
    }
    out_file.close();
}

bool test_cg_smooth_0() {
    // create a 3D grid for test
    loc_I = 50;
    loc_J = 50;
    loc_K = 50;
    smooth_lp = 0.1;
    smooth_lr = 0.1;
    smooth_lt = 0.1;
    // checkerboard pattern size
    int I = 10;
    int J = 10;
    int K = 10;

    CUSTOMREAL* arr_in = new CUSTOMREAL[loc_I*loc_J*loc_K];
    CUSTOMREAL* arr_out = new CUSTOMREAL[loc_I*loc_J*loc_K];



    for (int k=0; k<loc_K; k++){
        for (int j=0; j<loc_J; j++){
            for (int i=0; i<loc_I; i++){
                // set a checkerboard pattern in 3D
                if ((i/I+j/J+k/K)%2 == 0) {
                    arr_in[i+j*loc_I+k*loc_I*loc_J] = 1.0;
                } else {
                    arr_in[i+j*loc_I+k*loc_I*loc_J] = 0.0;
                }
            }
        }
    }

    // run smoothing
    CG_smooth(arr_in, arr_out, smooth_lr, smooth_lt, smooth_lp);

    // output initial and final grid in file
    out_file(arr_in, "test_cg_smooth_0_in");
    out_file(arr_out, "test_cg_smooth_0_out");

    // delete arrays
    delete[] arr_in;
    delete[] arr_out;

    // return true if test passed
    return true;
}

int main(int argc, char **argv) {
    if (test_cg_smooth_0()) {
        std::cout << "test_cg_smooth_0 passed" << std::endl;
    } else {
        std::cout << "test_cg_smooth_0 failed" << std::endl;
    }


    return 0;

}
