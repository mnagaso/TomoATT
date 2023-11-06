#include "eikonal_solver_2d.h"

void prepare_teleseismic_boundary_conditions(InputParams& IP, Grid& grid, IO_utils& io) {

    // teleseismic source is not supported in the public version
    //warning_teleseismic_use();
}

//
// here PlanGrid class is omitted for public version v2
//



void load_2d_traveltime(InputParams& IP, Source& src, Grid& grid, IO_utils& io) {

    warning_teleseismic_use();

}


