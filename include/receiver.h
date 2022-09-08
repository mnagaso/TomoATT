#ifndef RECEIVER_H
#define RECEIVER_H

#include <vector>
#include "input_params.h"
#include "grid.h"

class Receiver {
public:
    Receiver();
    ~Receiver();

    void calculate_arrival_time(InputParams&, Grid&);
    // regional source
    CUSTOMREAL calculate_adjoint_source(InputParams&);
    // teleseismic source
    CUSTOMREAL calculate_adjoint_source_teleseismic(InputParams&);
    // Gradient of traveltime
    void calculate_T_gradient(InputParams&, Grid&);
    // initialize variables for source relocation
    void init_vars_src_reloc(InputParams&);
    // approximated optimal origin time
    void calculate_optimal_origin_time(InputParams&);
    // divide optimal origin time by summed weight
    void divide_optimal_origin_time_by_summed_weight(InputParams&);
    // Gradient of objective function
    void calculate_grad_obj_src_reloc(InputParams&);
    // update source location
    void update_source_location(InputParams&);

private:
    void interpolate_travel_time(Grid&, SrcRec&);
    void calculate_T_gradient_one_rec(Grid&, SrcRec&);
};

#endif // RECEIVER_H