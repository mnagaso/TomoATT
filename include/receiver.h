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
    void interpolate_and_store_arrival_times_at_rec_position(InputParams&, Grid&, const std::string&);       // only for common receiver differential traveltime

    // adjoint source
<<<<<<< HEAD
    void calculate_adjoint_source(InputParams&, const std::string&);
    // objective function and residual
    std::vector<CUSTOMREAL> calculate_obj_and_residual(InputParams&);
    // teleseismic source
    // CUSTOMREAL calculate_adjoint_source_teleseismic(InputParams&, const std::string&);
=======
    std::vector<CUSTOMREAL> calculate_adjoint_source(InputParams&, const std::string&);
>>>>>>> eac89ca3a2cefe512853ef92a09e2f258e485fbb
    // Gradient of traveltime
    void calculate_T_gradient(InputParams&, Grid&, const std::string&);
    // initialize variables for source relocation
    void init_vars_src_reloc(InputParams&);
    // divide optimal origin time by summed weight
    void divide_optimal_origin_time_by_summed_weight(InputParams&);
    // Gradient of objective function
    void calculate_grad_obj_src_reloc(InputParams&, const std::string&);
    // objective function
    void calculate_obj_reloc(InputParams&, int);
    // update source location
    void update_source_location(InputParams&, Grid&);

private:
    CUSTOMREAL interpolate_travel_time(Grid&, InputParams&, std::string, std::string);
    void calculate_T_gradient_one_rec(Grid&, SrcRecInfo&, CUSTOMREAL*);
};

#endif // RECEIVER_H