#ifndef ITERATOR_SELECTOR_H
#define ITERATOR_SELECTOR_H

#include <memory>
#include "iterator.h"
#include "iterator_legacy.h"
#include "iterator_level.h"


void select_iterator(InputParams& IP, Grid& grid, Source& src, IO_utils& io, bool first_init, bool is_teleseismic, std::unique_ptr<Iterator>& It, bool is_second_run) {
    // initialize iterator object

    if (!is_teleseismic){
        if (IP.get_sweep_type() == SWEEP_TYPE_LEGACY) {
            if (IP.get_stencil_order() == 1){
                if (IP.get_stencil_type() == 1)
                    It = std::make_unique<Iterator_legacy_1st_order_upwind>(IP, grid, src, io, first_init, is_teleseismic, is_second_run);
                else
                    It = std::make_unique<Iterator_legacy_1st_order>(IP, grid, src, io, first_init, is_teleseismic, is_second_run);
            } else if (IP.get_stencil_order() == 3){
                It = std::make_unique<Iterator_legacy_3rd_order>(IP, grid, src, io, first_init, is_teleseismic, is_second_run);
                if (IP.get_stencil_type() == 1)
                    std::cout << "WARNING: Upwind Stencil type not supported, using non upwind scheme (LF)" << std::endl;
            } else{
                std::cout << "ERROR: Stencil order not supported" << std::endl;
                exit(1);
            }
        } else if (IP.get_sweep_type() == SWEEP_TYPE_LEVEL){
            if (IP.get_stencil_order() == 1){
                It = std::make_unique<Iterator_level_1st_order>(IP, grid, src, io, first_init, is_teleseismic, is_second_run);
                if (IP.get_stencil_type() == 1)
                    std::cout << "WARNING: Upwind Stencil type not supported, using non upwind scheme (LF)" << std::endl;
            } else if (IP.get_stencil_order() == 3){
                It = std::make_unique<Iterator_level_3rd_order>(IP, grid, src, io, first_init, is_teleseismic, is_second_run);
                if (IP.get_stencil_type() == 1)
                    std::cout << "WARNING: Upwind Stencil type not supported, using non upwind scheme (LF)" << std::endl;
            } else{
                std::cout << "ERROR: Stencil order not supported" << std::endl;
                exit(1);
            }
        } else {
            std::cout << "ERROR: Sweep type not supported" << std::endl;
            exit(1);
        }
    } else { // teleseismic event
        if (IP.get_sweep_type() == SWEEP_TYPE_LEGACY) {
            if (IP.get_stencil_order() == 1){
                if (IP.get_stencil_type() == 1)
                    It = std::make_unique<Iterator_legacy_1st_order_upwind_tele>(IP, grid, src, io, first_init, is_teleseismic, is_second_run);
                else
                    It = std::make_unique<Iterator_legacy_1st_order_tele>(IP, grid, src, io, first_init, is_teleseismic, is_second_run);
            } else if (IP.get_stencil_order() == 3){
                It = std::make_unique<Iterator_legacy_3rd_order_tele>(IP, grid, src, io, first_init, is_teleseismic, is_second_run);
                if (IP.get_stencil_type() == 1)
                    std::cout << "WARNING: Upwind Stencil type not supported, using non upwind scheme (LF)" << std::endl;
            } else{
                std::cout << "ERROR: Stencil order not supported" << std::endl;
                exit(1);
            }
        } else if (IP.get_sweep_type() == SWEEP_TYPE_LEVEL){
            if (IP.get_stencil_order() == 1){
                It = std::make_unique<Iterator_level_1st_order_tele>(IP, grid, src, io, first_init, is_teleseismic, is_second_run);
                if (IP.get_stencil_type() == 1)
                    std::cout << "WARNING: Upwind Stencil type not supported, using non upwind scheme (LF)" << std::endl;
            } else if (IP.get_stencil_order() == 3){
                It = std::make_unique<Iterator_level_3rd_order_tele>(IP, grid, src, io, first_init, is_teleseismic, is_second_run);
                if (IP.get_stencil_type() == 1)
                    std::cout << "WARNING: Upwind Stencil type not supported, using non upwind scheme (LF)" << std::endl;
            } else{
                std::cout << "ERROR: Stencil order not supported" << std::endl;
                exit(1);
            }
        } else {
            std::cout << "ERROR: Sweep type not supported" << std::endl;
            exit(1);
        }
    }

}

#endif

