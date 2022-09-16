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

    // adjoint source
    CUSTOMREAL calculate_adjoint_source(InputParams&);
    // teleseismic source
    CUSTOMREAL calculate_adjoint_source_teleseismic(InputParams&);
    // Gradient of traveltime
    void calculate_T_gradient(InputParams&, Grid&);
    // initialize variables for source relocation
    void init_vars_src_reloc(InputParams&, std::vector<SrcRec>&);
    // approximated optimal origin time
    void calculate_optimal_origin_time(InputParams&, std::vector<SrcRec>&);
    // divide optimal origin time by summed weight
    void divide_optimal_origin_time_by_summed_weight(InputParams&, std::vector<SrcRec>&);
    // Gradient of objective function
    void calculate_grad_obj_src_reloc(InputParams&, std::vector<SrcRec>&);
    // update source location
    void update_source_location(InputParams&, Grid&, std::vector<SrcRec>&);

private:
    void interpolate_travel_time(Grid&, SrcRec&);
    void calculate_T_gradient_one_rec(Grid&, SrcRec&);
    void interpolate_differential_travel_time(Grid&, SrcRec&);
};

#endif // RECEIVER_H