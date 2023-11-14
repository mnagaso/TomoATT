#ifndef INV_GRID_H
#define INV_GRID_H

#include <algorithm>

#include "config.h"
#include "input_params.h"


// Base class for 1D inversion grid ()
class InvGrid1dBase {
public:
    InvGrid1dBase(){};
    virtual ~InvGrid1dBase(){
        delete[] arr;
    };

    // copy constructor
    InvGrid1dBase(const InvGrid1dBase& other){
        n = other.n;
        dinv = other.dinv;
        dinv_l = other.dinv_l;
        arr = new CUSTOMREAL[n];
        std::copy(other.arr, other.arr+n, arr);
    }

    // assignment operator
    InvGrid1dBase& operator=(const InvGrid1dBase& other){
        if (this != &other){
            n = other.n;
            dinv = other.dinv;
            dinv_l = other.dinv_l;
            arr = new CUSTOMREAL[n];
            std::copy(other.arr, other.arr+n, arr);
        }
        return *this;
    }

    CUSTOMREAL* arr = nullptr; // 1d or 2d array storing the grid coordinates
    int n; // number of grid points
    CUSTOMREAL dinv; // grid spacing
    CUSTOMREAL dinv_l; // amount of shift for each inversion grid
};

// Derived class for 1D inversion grid
class InvGrid1d : public InvGrid1dBase {
public:
    InvGrid1d(InputParams&, const int, const CUSTOMREAL*); // for r grid
    InvGrid1d(InputParams&, const int, const CUSTOMREAL*, const CUSTOMREAL*, const int, const CUSTOMREAL*); // function overload for t and p grids
    ~InvGrid1d() override {};
};


// Base class for 3D inversion grid
class InvGrid {
public:
    InvGrid(InputParams&);
    ~InvGrid();

    void write_inversion_grid_to_file();

    InvGrid1dBase r;
    InvGrid1dBase t;
    InvGrid1dBase p;
    InvGrid1dBase r_ani;
    InvGrid1dBase t_ani;
    InvGrid1dBase p_ani;

private:
    void get_inv_grid_params(InputParams&);
};



#endif // INV_GRID_H
