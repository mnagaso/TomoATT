#ifndef MAIN_ROUTINES_CALLING_H
#define MAIN_ROUTINES_CALLING_H

#include <iostream>
#include <memory>
#include "mpi_funcs.h"
#include "config.h"
#include "utils.h"
#include "input_params.h"
#include "grid.h"
#include "io.h"
#include "main_routines_inversion_mode.h"
#include "model_optimization_routines.h"
#include "main_routines_earthquake_relocation.h"
#include "iterator_selector.h"
#include "iterator.h"
#include "iterator_legacy.h"
#include "iterator_level.h"
#include "source.h"
#include "receiver.h"
#include "kernel.h"
#include "model_update.h"
#include "lbfgs.h"


// prepare header line of objective_funciton.txt
inline void prepare_header_line(InputParams &IP, std::ofstream &out_main) {
    // prepare output for iteration status
    if(myrank == 0 && id_sim ==0){
        out_main.open(output_dir + "/objective_function.txt");
        if (optim_method == GRADIENT_DESCENT || optim_method == HALVE_STEPPING_MODE){
            bool have_abs    = ( IP.data_type.find("abs")    != IP.data_type.end() );
            bool have_cs_dif = ( IP.data_type.find("cs_dif") != IP.data_type.end() );
            bool have_cr_dif = ( IP.data_type.find("cr_dif") != IP.data_type.end() );
            bool have_tele   = ( IP.data_type.find("tele")   != IP.data_type.end() );


            out_main << std::setw(8) << std::right << "# iter,";
            if (optim_method == HALVE_STEPPING_MODE)
                out_main << std::setw(8) << std::right << "subiter,";
            std::string tmp = "obj(";
            tmp.append(std::to_string(IP.N_data));
            tmp.append("),");
            out_main << std::setw(20) << tmp;
            if (have_abs){
                tmp = "obj_abs(";
                tmp.append(std::to_string(IP.N_abs_local_data));
                tmp.append("),");
                out_main << std::setw(20) << tmp;
            }
            if (have_cs_dif){
                tmp = "obj_cs_dif(";
                tmp.append(std::to_string(IP.N_cs_dif_local_data));
                tmp.append("),");
                out_main << std::setw(20) << tmp;
            }
            if (have_cr_dif){
                tmp = "obj_cr_dif(";
                tmp.append(std::to_string(IP.N_cr_dif_local_data));
                tmp.append("),");
                out_main << std::setw(20) << tmp;
            }
            if (have_tele){
                tmp = "obj_tele(";
                tmp.append(std::to_string(IP.N_teleseismic_data));
                tmp.append("),");
                out_main << std::setw(20) << tmp;
            }
            out_main << std::setw(20) << "misfit,";
            if (have_abs){
                out_main << std::setw(20) << "misfit_abs,";
            }
            if (have_cs_dif){
                out_main << std::setw(20) << "misfit_cs_dif,";
            }
            if (have_cr_dif){
                out_main << std::setw(20) << "misfit_cr_dif,";
            }
            if (have_tele){
                out_main << std::setw(20) << "misfit_tele,";
            }
            out_main << std::setw(20) << "step_size," << std::endl;
        } else if (optim_method == LBFGS_MODE)
            out_main << std::setw(6)  << "it,"        << std::setw(6)  << "subit,"  << std::setw(16) << "step_size," << std::setw(16) << "qpt," << std::setw(16) << "v_obj_new," \
                     << std::setw(16) << "v_obj_reg," << std::setw(16) << "q_new,"  << std::setw(16) << "q_k,"       << std::setw(16) << "td,"  << std::setw(16) << "tg," \
                     << std::setw(16) << "c1*q_k,"    << std::setw(16) << "c2*q_k," << std::setw(6)  << "step ok"    << std::endl;

    }
}


//
// run forward-only or inversion mode
//
inline void run_forward_only_or_inversion(InputParams &IP, Grid &grid, IO_utils &io) {

    // for check if the current source is the first source
    bool first_src = true;

    if(myrank == 0)
        std::cout << "id_sim: " << id_sim << ", size of src_map: " << IP.src_map.size() << std::endl;

    // prepare objective_function file
    std::ofstream out_main; // close() is not mandatory
    prepare_header_line(IP, out_main);

    if (subdom_main && id_sim==0 && IP.get_is_output_model_dat()==1) {
        io.write_concerning_parameters(grid, 0, IP);
    }

    // output station correction file (only for teleseismic differential data)
    IP.write_station_correction_file(0);

    synchronize_all_world();

    /////////////////////
    // loop for inversion
    /////////////////////

    bool line_search_mode = false; // if true, run_simulation_one_step skips adjoint simulation and only calculates objective function value

    // objective function for all src
    CUSTOMREAL v_obj = 0.0, old_v_obj = 0.0;
    std::vector<CUSTOMREAL> v_obj_misfit(10, 0.0);

    for (int i_inv = 0; i_inv < IP.get_max_iter_inv(); i_inv++) {

        if(myrank == 0 && id_sim ==0){
            std::cout << "iteration " << i_inv << " starting ... " << std::endl;
        }

        old_v_obj = v_obj;

        // prepare inverstion iteration group in xdmf file
        io.prepare_grid_inv_xdmf(i_inv);

        ///////////////////////////////////////////////////////
        // run (forward and adjoint) simulation for each source
        ///////////////////////////////////////////////////////

        // run forward and adjoint simulation and calculate current objective function value and sensitivity kernel for all sources
        line_search_mode = false;
        // skip for the mode with sub-iteration
        if (i_inv > 0 && optim_method != GRADIENT_DESCENT) {
        } else {
            v_obj_misfit = run_simulation_one_step(IP, grid, io, i_inv, first_src, line_search_mode);
            v_obj = v_obj_misfit[0];
        }

        // wait for all processes to finish
        synchronize_all_world();

        // check if v_obj is nan
        if (std::isnan(v_obj)) {
            if (myrank == 0)
                std::cout << "v_obj is nan, stop inversion" << std::endl;

            // stop inversion
            break;
        }

        // output src rec file with the result arrival times
        if (i_inv == IP.get_max_iter_inv()-1)
            IP.write_src_rec_file(i_inv);

        ///////////////
        // model update
        ///////////////
        if (IP.get_run_mode() == DO_INVERSION) {
            if (optim_method == GRADIENT_DESCENT)
                model_optimize(IP, grid, io, i_inv, v_obj, old_v_obj, v_obj_misfit, first_src, out_main);
            else if (optim_method == LBFGS_MODE)
                model_optimize_lbfgs(IP, grid, io, i_inv, v_obj, first_src, out_main);
            else if (optim_method == HALVE_STEPPING_MODE)
                model_optimize_halve_stepping(IP, grid, io, i_inv, v_obj, first_src, out_main);
        }

        // output station correction file (only for teleseismic differential data)
        IP.write_station_correction_file(i_inv + 1);

        // output updated model
        if (subdom_main && id_sim==0) {
            //io.change_xdmf_obj(0); // change xmf file for next src
            io.change_group_name_for_model();

            // write out model info
            if (IP.get_is_output_in_process() || i_inv >= IP.get_max_iter_inv() - 2){
                io.write_vel(grid, i_inv+1);
                io.write_xi( grid, i_inv+1);
                io.write_eta(grid, i_inv+1);
            }
            //io.write_zeta(grid, i_inv); // TODO

            if (IP.get_is_verbose_output()){
                io.write_a(grid,   i_inv+1);
                io.write_b(grid,   i_inv+1);
                io.write_c(grid,   i_inv+1);
                io.write_f(grid,   i_inv+1);
                io.write_fun(grid, i_inv+1);
            }

            // output model_parameters_inv_0000.dat
            if (IP.get_is_output_model_dat() \
            && (IP.get_is_output_in_process() || i_inv >= IP.get_max_iter_inv() - 2))
                io.write_concerning_parameters(grid, i_inv + 1, IP);

        }

        // writeout temporary xdmf file
        io.update_xdmf_file();

        // wait for all processes to finish
        synchronize_all_world();

    } // end loop inverse

    // close xdmf file
    io.finalize_data_output_file();

}


// run earthquake relocation mode
inline void run_earthquake_relocation(InputParams& IP, Grid& grid, IO_utils& io) {

    Receiver recs;
    int nrec_total = IP.rec_map_all.size();
    // broadcast
    broadcast_i_single_inter_and_intra_sim(nrec_total, 0);

    // calculate traveltime for each receiver (swapped from source) and write in output file
    calculate_traveltime_for_all_src_rec(IP, grid, io);

    //std::cout << "cpk-sp, id_sim: " << id_sim << ", myrank: " << myrank << ", Nrec: " << IP.rec_map.size() << std::endl;
    //for(auto iter = IP.rec_map.begin(); iter != IP.rec_map.end(); iter++){
    //    std::cout << "cpk-sp, id_sim: " << id_sim << ", myrank: " << myrank << ", name: " << iter->first << std::endl;
    //}

    // prepare output for iteration status
    std::ofstream out_main;
    if(myrank == 0 && id_sim ==0){
        out_main.open(output_dir + "/objective_function_reloc.txt");
        out_main << std::setw(8) << std::right << "# iter,";
        out_main << std::setw(16) << std::right << "N_reloc,";
        out_main << std::setw(16) << std::right << "N_located,";
        out_main << std::setw(16) << std::right << "obj_weighted,";
        // out_main << std::setw(12) << std::right << "obj_noweight,";
        out_main << std::endl;
    }

    // objective function and its gradient
    CUSTOMREAL v_obj      = 999999999.0;
    CUSTOMREAL v_obj_old  = 0.0;
    CUSTOMREAL v_obj_grad = 0.0;
    int        i_iter     = 0;

    // iterate
    while (true) {

        v_obj_old  = v_obj;
        v_obj      = 0.0;
        v_obj_grad = 0.0;

        // determine which earthquake should be located
        IP.name_for_reloc.clear();
        for(auto iter = IP.rec_map.begin(); iter != IP.rec_map.end(); iter++){
            if (!iter->second.is_stop)
                IP.name_for_reloc.push_back(iter->first);
        }

        // calculate gradient of objective function at sources
        calculate_gradient_objective_function(IP, grid, io, i_iter);

        // update source location
        recs.update_source_location(IP, grid);

        synchronize_all_world();

        // calculate sum of objective function and gradient
        for (auto iter = IP.rec_map.begin(); iter != IP.rec_map.end(); iter++) {
            v_obj      += iter->second.vobj_src_reloc;
            v_obj_grad += iter->second.vobj_grad_norm_src_reloc;
        }

        // reduce
        //allreduce_cr_sim_single_inplace(v_obj);
        //allreduce_cr_sim_single_inplace(v_obj_grad);

        // check convergence
        int count_loc = 0;
        bool finished = false;

        if (subdom_main && id_subdomain==0) {
            if (IP.name_for_reloc.size() == 0){
                std::cout << "DEBUG: finished relocation because all receivers have been located." << std::endl;
                finished = true;
            }

            if (i_iter >= N_ITER_MAX_SRC_RELOC){
                std::cout << "DEBUG: finished relocation because iteration number is larger than " << N_ITER_MAX_SRC_RELOC << std::endl;
                finished = true;
            }
            allreduce_bool_single_inplace_sim(finished); //LAND
        }

        synchronize_all_world();

        // check if all processes have finished
        broadcast_bool_inter_and_intra_sim(finished, 0);

        synchronize_all_world();

        // output location information
        if(id_sim == 0 && myrank == 0){
            // number of receiver which have been completed
            int n_relocated = IP.rec_map.size() - IP.name_for_reloc.size();

            // write objective function
            std::cout << "iteration: " << i_iter << ", objective function: "              << v_obj
                                                 << ", objective function old: "          << v_obj_old
                                                 << ", norm grad of relocating: "         << v_obj_grad // BUG: strangely high when simultaneous run
                                                 << ", average norm grad of relocating: " << v_obj_grad/IP.name_for_reloc.size()
                                                 << ", v_obj/n_src: "                     << v_obj/nrec_total
                                                 << ", diff_v/v_obj_old "                 << std::abs(v_obj-v_obj_old)/v_obj_old << std::endl;
            std::cout << "Earthquakes require location: " << n_relocated << " / " << nrec_total << " completed." << std::endl;

            // the last 10 sources under location
            if (nrec_total - count_loc < 10){
                std::cout << "Source under location. names: ";
                for (int i = 0; i < (int)IP.name_for_reloc.size(); i++){
                    std::cout << IP.name_for_reloc[i] << ", ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        // write objective functions
        if (myrank==0 && id_sim==0) {
            out_main << std::setw(8) << std::right << i_iter - 1;
            out_main << std::setw(16) << std::right << (int)IP.rec_map.size();
            out_main << std::setw(16) << std::right << count_loc;
            out_main << std::setw(16) << std::right << v_obj;
            out_main << std::endl;
        }

        if (finished)
            break;

        // new iteration
        i_iter++;


    }


    // modify the receiver's location
    IP.modift_swapped_source_location();
    // write out new src_rec_file
    IP.write_src_rec_file(0);
    // close xdmf file
    io.finalize_data_output_file();


}


#endif // MAIN_ROUTINES_CALLING_H