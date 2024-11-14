#!/bin/bash

# Some script to plot figures for the output file of TomoATT

python prepare_files.py             # download the files for plotting
tar -xf files_for_plotting.tar.gz   # extract the files


# Test 1: plot velocity perturbation and azimuthal anisotropy fields to generate
# img/1_dep_vel.png     2D velocity perturbation at 20 km depth
# img/1_dep_ani.png     2D azimuthal anisotropy at 20 km depth
# img/1_sec_vel.png     2D velocity perturbation along vertical section
python 1_plot_model.py

# Test 2: plot traveltime and adjoint fields to generate
# img/2_dep_time.png     2D traveltime field at 20 km depth
# img/2_sec_time.png     2D traveltime field along vertical section
# img/2_dep_adjoint.png  2D adjoint field at XX depth
# img/2_sec_adjoint.png  2D adjoint field along vertical section
python 2_plot_time_field.py

# Test 3: plot kernels to generate
# img/3_dep_Ks_inv_0007.png             Ks:                 original kernel w.r.t slowness at 20 km depth
# img/3_dep_Ks_density_inv_0007.png     Kden:               kernel density w.r.t slowness at 20 km depth
# img/3_dep_Ks_over_Kden_inv_0007.png   Ks/Kden^{\zeta}:    normalized kernel w.r.t slowness at 20 km depth
# img/3_dep_Ks_update_inv_0007.png                          smoothed normalized kernel w.r.t slowness at 20 km depth
# img/3_sec_Ks_inv_0007.png             Ks:                 original kernel w.r.t slowness along vertical section
# img/3_sec_Ks_density_inv_0007.png     Kden:               kernel density w.r.t slowness along vertical section
# img/3_sec_Ks_over_Kden_inv_0007.png   Ks/Kden^{\zeta}:    normalized kernel w.r.t slowness along vertical section
# img/3_sec_Ks_update_inv_0007.png                          smoothed normalized kernel w.r.t slowness along vertical section
# and the same for kernels w.r.t xi and eta (azimuthal anisotropy)
python 3_plot_kernel.py

# Test 4: plot earthquakes and stations to generate 
# img/4_earthquakes_and_stations.png     the location of earthquakes and stations
python 4_plot_earthquake_station.py

# Test 5: plot objective function reduction to generate
# img/5_objective_function_reduction.png     the reduction of objective function value
python 5_plot_objective_function.py

# Test 6: plot traveltime residuals to generate
# img/6_data_residual.png     the traveltime residuals
python 6_plot_data_residual.py

# Test 7: plot inversion grid to generate
# img/7_inversion_grid.png     the inversion grid
python 7_plot_inversion_grid.py


