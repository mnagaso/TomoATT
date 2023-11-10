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
    InvGrid1dBase(const InvGrid1dBase& rhs) : arr(new CUSTOMREAL[rhs.n]), n(rhs.n) {
        std::copy(rhs.arr, rhs.arr + rhs.n, arr);
    };

    // assignment operator
    InvGrid1dBase& operator=(const InvGrid1dBase& rhs){
        if (this != &rhs){
            CUSTOMREAL* newArr = new CUSTOMREAL[rhs.n];
            std::copy(rhs.arr, rhs.arr + rhs.n, newArr);
            delete[] arr;
            arr = newArr;
            n = rhs.n;
        }

        return *this;
    };

    CUSTOMREAL* arr; // 1d or 2d array storing the grid coordinates
    int n; // number of grid points
    CUSTOMREAL dinv; // grid spacing
    CUSTOMREAL dinv_l; // amount of shift for each inversion grid
};


// Derived class for 1D inversion grid (regular grid)
class InvGrid1dRegular : public InvGrid1dBase {
public:
    InvGrid1dRegular(InputParams&, const int, const CUSTOMREAL, const CUSTOMREAL); // for r grid
    InvGrid1dRegular(InputParams&, const int, const CUSTOMREAL, const CUSTOMREAL,
                                   const int); // function overload for t and p grids
    ~InvGrid1dRegular() override {};
};


// Derived class for 1D inversion grid (flexible grid)
class InvGrid1dFlexible : public InvGrid1dBase {
public:
    InvGrid1dFlexible(InputParams&, const int, const CUSTOMREAL*, const int); // for r grid
    InvGrid1dFlexible(InputParams&, const int, const CUSTOMREAL*, const int,
                                   const int); // function overload for t and p grids
    ~InvGrid1dFlexible() override {};
};


// Derived class for 1D inversion grid (trapezoidal grid)
class InvGrid1dTrapezoidal : public InvGrid1dBase {
public:
    InvGrid1dTrapezoidal(InputParams&, const int, const CUSTOMREAL, const CUSTOMREAL*, const int, const int); // function overload for t and p grids
    ~InvGrid1dTrapezoidal() override {};
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
