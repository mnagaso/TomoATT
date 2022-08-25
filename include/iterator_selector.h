#ifndef ITERATOR_SELECTOR_H
#define ITERATOR_SELECTOR_H

#include <memory>
#include "iterator.h"
#include "iterator_legacy.h"
#include "iterator_level.h"


std::unique_ptr<Iterator> select_iterator(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic) {
    // initialize iterator object
    std::unique_ptr<Iterator> It;

    if (!is_teleseismic){
        if (IP.get_sweep_type() == 0) {
            if (IP.get_stencil_order() == 1)
                It = std::make_unique<Iterator_legacy_1st_order>(IP, grid, src, io, first_init, is_teleseismic);
            else if (IP.get_stencil_order() == 3)
                It = std::make_unique<Iterator_legacy_3rd_order>(IP, grid, src, io, first_init, is_teleseismic);
            else{
                std::cout << "ERROR: Stencil order not supported" << std::endl;
                exit(1);
            }
        } else if (IP.get_sweep_type() == 2){
            if (IP.get_stencil_order() == 1)
                It = std::make_unique<Iterator_level_1st_order>(IP, grid, src, io, first_init, is_teleseismic);
            else if (IP.get_stencil_order() == 3)
                It = std::make_unique<Iterator_level_3rd_order>(IP, grid, src, io, first_init, is_teleseismic);
            else{
                std::cout << "ERROR: Stencil order not supported" << std::endl;
                exit(1);
            }
        }
    } else { // teleseismic event
        if (IP.get_sweep_type() == 0) {
            if (IP.get_stencil_order() == 1)
                It = std::make_unique<Iterator_legacy_1st_order_tele>(IP, grid, src, io, first_init, is_teleseismic);
            else if (IP.get_stencil_order() == 3)
                It = std::make_unique<Iterator_legacy_3rd_order_tele>(IP, grid, src, io, first_init, is_teleseismic);
            else{
                std::cout << "ERROR: Stencil order not supported" << std::endl;
                exit(1);
            }
        } else if (IP.get_sweep_type() == 2){
            if (IP.get_stencil_order() == 1)
                It = std::make_unique<Iterator_level_1st_order_tele>(IP, grid, src, io, first_init, is_teleseismic);
            else if (IP.get_stencil_order() == 3)
                It = std::make_unique<Iterator_level_3rd_order_tele>(IP, grid, src, io, first_init, is_teleseismic);
            else{
                std::cout << "ERROR: Stencil order not supported" << std::endl;
                exit(1);
            }
        }
    }

    return It;
}

#endif

