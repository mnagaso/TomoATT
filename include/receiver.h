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

private:
    void interpolate_travel_time(Grid&, SrcRec&);

    void interpolate_differential_travel_time(Grid&, SrcRec&);
};

#endif // RECEIVER_H