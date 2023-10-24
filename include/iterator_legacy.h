#ifndef ITERATOR_legacy_H
#define ITERATOR_legacy_H

#include "iterator.h"


class Iterator_legacy : public Iterator {
public:
    Iterator_legacy(InputParams&, Grid&, Source&, IO_utils&, const std::string&, bool, bool, bool);
protected:
    void do_sweep_adj(int, Grid&, InputParams&) override ; // do sweeping for adjoint routine
    virtual void do_sweep(int, Grid&, InputParams&) {}; // do sweeping
    virtual void do_sweep_attenuation(int, Grid&, InputParams&) {}; // do sweeping for attenuation
};


class Iterator_legacy_tele : public Iterator {
public:
    Iterator_legacy_tele(InputParams& , Grid&, Source&, IO_utils&, const std::string&, bool, bool, bool);
protected:
    void do_sweep_adj(int, Grid&, InputParams&) override ; // do sweeping for adjoint routine
    virtual void do_sweep(int, Grid&, InputParams&) {}; // do sweeping
    virtual void do_sweep_attenuation(int, Grid&, InputParams&) {}; // do sweeping for attenuation
};


class Iterator_legacy_1st_order : public Iterator_legacy {

public:
    Iterator_legacy_1st_order(InputParams&, Grid&, Source&, IO_utils&, const std::string&, bool, bool, bool);
private:
    void do_sweep(int, Grid&, InputParams&) override ; // do sweeping
    void do_sweep_attenuation(int, Grid&, InputParams&) override; // do sweeping for attenuation

};

class Iterator_legacy_3rd_order : public Iterator_legacy {

public:
    Iterator_legacy_3rd_order(InputParams&, Grid&, Source&, IO_utils&, const std::string&, bool, bool, bool);
    //void run_iteration_forward(InputParams&, Grid&, IO_utils&, bool&);
private:
    void do_sweep(int, Grid&, InputParams&) override ; // do sweeping
    void do_sweep_attenuation(int, Grid&, InputParams&) override; // do sweeping for attenuation

};

class Iterator_legacy_1st_order_upwind : public Iterator_legacy {

public:
    Iterator_legacy_1st_order_upwind(InputParams&, Grid&, Source&, IO_utils&, const std::string&, bool, bool, bool);
private:
    void do_sweep(int, Grid&, InputParams&) override ; // do sweeping
    void do_sweep_attenuation(int, Grid&, InputParams&) override; // do sweeping for attenuation

};

class Iterator_legacy_1st_order_tele : public Iterator_legacy_tele {

public:
    Iterator_legacy_1st_order_tele(InputParams&, Grid&, Source&, IO_utils&, const std::string&, bool, bool, bool);
    //void run_iteration_forward(InputParams&, Grid&, IO_utils&, bool&);
private:
    void do_sweep(int, Grid&, InputParams&) override ; // do sweeping
    void do_sweep_attenuation(int, Grid&, InputParams&) override; // do sweeping for attenuation

};

class Iterator_legacy_3rd_order_tele : public Iterator_legacy_tele {

public:
    Iterator_legacy_3rd_order_tele(InputParams&, Grid&, Source&, IO_utils&, const std::string&, bool, bool, bool);
    //void run_iteration_forward(InputParams&, Grid&, IO_utils&, bool&);
private:
    void do_sweep(int, Grid&, InputParams&) override ; // do sweeping
    void do_sweep_attenuation(int, Grid&, InputParams&) override; // do sweeping for attenuation

};

class Iterator_legacy_1st_order_upwind_tele : public Iterator_legacy_tele {

public:
    Iterator_legacy_1st_order_upwind_tele(InputParams&, Grid&, Source&, IO_utils&, const std::string&, bool, bool, bool);
    //void run_iteration_forward(InputParams&, Grid&, IO_utils&, bool&);
private:
    void do_sweep(int, Grid&, InputParams&) override ; // do sweeping
    void do_sweep_attenuation(int, Grid&, InputParams&) override; // do sweeping for attenuation

};

#endif // ITERATOR_LEGACY_H