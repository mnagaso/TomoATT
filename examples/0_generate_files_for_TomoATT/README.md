# Generating model and data files for examples 

This is an example to generating model file and data file for the examples in the folder. Also users can mofidy the script for generating the own model for TomoATT.

1. Run all cells of `1_build_initial_model.ipynb` or python script `1_build_initial_model.py` for creating the initial model:
    - model with a mesh of 61 * 61 * 61 grid nodes (models/model_init_N61_61_61.h5)
    - the region is [-10 km,50 km] * [30, 32] * [30, 32] in depth, latitude, and longitude

The velocity of this model is:
    - 6.0 km/s,             if          dep < 0 km;
    - 6.0 + dep/20 km/s,    if 0  km <  dep < 40 km;
    - 8.0 km/s,             if 40 km <  dep;

2. Run all cells of `2_build_checkerboard_model.ipynb` or python script `2_build_checkerboard_model.py` for creating the checkerboard model:
    - model with a mesh of 61 * 61 * 61 grid nodes (models/model_ckb_N61_61_61.h5)

<!-- The checkerboard model is built by assigning perturbations to velocity and azimuthal anisotropy. You can run `` to show the horizontal sections in:
    - img/ckb_model_ani.jpg
    - img/ckb_model_vel.jpg -->

3. Run all cells of `3_generate_src_rec_file.ipynb` or python script `3_generate_src_rec_file.py` for creating the data file `src_rec_files/src_rec_config.dat`, which including:
    - absoulte traveltime data;
    - common-source differential traveltime data;
    - common-receiver differential traveltime data;
    
4. You can run `src_rec_visulization_pygmt_is_needed.ipynb` or python script to check the distribution of earthquakes (dot) and stations (triangle) in img/src_rec.jpg

You can run `bash run_this_example.sh` to proceed steps 1-3



