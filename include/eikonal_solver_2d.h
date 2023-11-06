#ifndef EIKONAL_SOLVER_2D_H
#define EIKONAL_SOLVER_2D_H

#include <fstream>
#include <string>
#include <iostream>

#include "utils.h"
#include "input_params.h"
#include "grid.h"
#include "1d_models.h"
#include "mpi_funcs.h"
#include "io.h"

//
// here PlanGrid class is omitted for public version v2
//

void prepare_teleseismic_boundary_conditions(InputParams&, Grid&, IO_utils&);
//void run_2d_solver(InputParams&, Source&, IO_utils&);
//void interp2d(PlainGrid&, CUSTOMREAL, CUSTOMREAL, CUSTOMREAL&);
void load_2d_traveltime(InputParams&, Source&, Grid&, IO_utils&);
//std::string get_2d_tt_filename(const std::string&, Source&);

#endif // EIKONAL_SOLVER_2D_H