# TomoATT User Manual version 1.0

## Introduction

TomoATT is a library which implements an eikonal equation solver based Adjoint-state Travel-Time Tomography for a very large-scale computation, following a published article [Ping Tong (2021)](https://doi.org/10.1029/2021JB021818) and [Jing Chen (2022)](add_here_when_published).
The slowness field and anisotropy fields are computed in spherical coordinate system.

Thanks to the efficiency of an eikonal equation solver, the computation of the travel-time is very fast and requires less amount of computational resources.
As an input data for TomoATT is travel times at seismic stations, we can easily prepare a great amount of input data for the computation.

This library is developped to be used for modeling a very-large domain. For this purpose, 3-layer parallelization is applied, which are:
- layer 1: simulutaneous run parallelization (travel times for multiple seismic sources may be calculated simultaneously)
- layer 2: subdomain decomposition (If the number of computational nodes requires too large memory, we can separate the domain into subdomains and run each subdomain in a separate compute node)
- layer 3: sweeping parallelization (in each subdomain, sweeping layers are also parallelized)

The details of the parallelization method applied in this library are described in the paper [Miles Detrixhe and Frédéric Gibou (2016)](https://doi.org/10.1016/j.jcp.2016.06.023).

Regional events (sources within the global domain) and teleseismic events (sources outside the global domain) may be used for inversion.

---

## Input files

TomoATT requires 3 input files : 
1. input parameter file (setup file for simulation parameters)
2. source receiver file (source and receiver definitions and observation arrival times)
3. initial model file (3d initial model)


### 1. input parameter file

All the necessary parameters for setuping a calculation are described in input parameter file in [yaml format](https://en.wikipedia.org/wiki/YAML).

Below is an example of input parameter file for making a forward simulation.

``` yaml
version : 2

domain :
  min_max_dep : [5740.6370,5790.6370] # depth in km
  min_max_lat : [45.0,55.0]           # latitude in degree
  min_max_lon : [35.0,45.0]           # longitude in degree
  n_rtp : [40,40,40]                  # number of nodes

source :
  src_rec_file : 'src_rec_test.dat' # source receiver file

model :
  init_model_path : './test_model_init.h5' # path to initial model file

inversion :
  run_mode : 0 # 0 for forward simulation only, 1 for inversion

parallel :
  n_sims : 1         # number of simultaneous run
  ndiv_rtp : [2,2,1] # number of subdomains
  nproc_sub : 1      # number of subprocess used for each subdomain

calculation :
  convergence_tolerance : 1e-6
  max_iterations : 100
  stencil_order : 3 # 1 or 3
  sweep_type : 2   # 0: legacy, 1: cuthill-mckee. 2: cuthill-mckee with shm parallelization
```

There are categories and sub-categories with setup parameters.
Below is the explanation of each category and sub-category.
If a category tag or sub-category tag is missing from the input parameter file, the default value will be used.

#### domain :
The domain category is for setting a global domain of the simulation.
- `min_max_dep` : `[dep1, dep2]` minimum and maximum depth of the domain in kilo meter
- `min_max_lat` : `[lat1, lat2]` minimum and maximum latitude in degree
- `min_max_lon` : `[lon1, lon2]` minimum and maximum longitude in degree
- `n_rtp` : `[nr, nt, np]` number of computation nodes for r (depth),t (latitude) and p (longitude) direction

#### source :
Source category is used for setting about source and receiver definitions.
- `src_rec_file` : for specifying a path to a source receiver file (details of this file is in the following section). The calculated travel time at receivers will be output as a new source receiver file in OUTPUT_FILES directory.
- `swap_src_rec` : `0` or `1`. Set 1 for using receivers as sources, which reduces computation time when the number of sources are larger than the number of receivers. (usual when using a large dataset.)

#### model :
- `init_model_path` : File path for 3D initial model file. Details will be explained in the following section.
-  `model_1d_name` : the name of 1d model used for teleseismic tomography. `ak135` and `iasp91` are available. User defined model can be used by modifying the file `1d_models.h` in `include` directory, by changing `model_1d_prem` part and setting this `model_1d_name` parameter as `user_defined`. 
#### inversion :
- `run_mode` : 
  - `0` for running only a forward simulation. 
  - `1` for do inversion. 
  - `2` for only precalculation of 2d traveltime field for teleseismic sources. This is an optional step for teleseismic case. If 2d traveltime field is not precalculated, it will be calculated during mode `0` or `1`. This is useful for calculations on HPC as 2d eikonal solver has not been parallelized. Users can run it on a local machine beforehand to reduce CPU time on HPC.
  - `3` run earthquake relocation. For using this mode, `sources : swap_src_rec` need to be set to 1. Otherwise the program will exit with an error message.
- `n_inversion_grid` :  the number of inversion grid.
- `n_inv_dep_lat_lon` : the numbers of inversion grids for r, t and p direction.
- `min_max_dep_inv` : `[dep1, dep2]` minimum and maximum depth of inversion grid in kilo meter.
- `min_max_lat_inv` : `[lat1, lat2]` minimum and maximum latitude of inversion grid in degree.
- `min_max_lon_inv` : `[lon1, lon2]` minimum and maximum longitude of inversion grid in degree.

Currently, TomoATT provide two ways for defining the inversion grid. One is a regular grid (even intervals for all axis) type and another uses used defined intervals for each axis.

For setting the user defined intervals, please use the flags below. 
- `type_dep_inv` : `0` for regular grid (default), `1` for user defined intervals.
- `type_lat_inv` : `0` for regular grid (default), `1` for user defined intervals.
- `type_lon_inv` : `0` for regular grid (default), `1` for user defined intervals.

and the coordinates where the main inversion grids are placed,
- `dep_inv` : `[z1,z2,z3,...]` for user defined intervals.
- `lat_inv` : `[y1,y2,y3,...]` for user defined intervals.
- `lon_inv` : `[x1,x2,x3,...]` for user defined intervals.

other parameers for inversion seetting.
- `max_iterations_inv` :  The limit of iteration number for inversion.
- `step_size` :  Maximum step size ratio for updating model.
- `smooth_method` : `0` or `1`. 0 for multigrid parametrization ([Ping Tong 2019](https://doi.org/10.1093/gji/ggz151)). 1 for laplacian smoothing with CG method.
- `l_smooth_rtp` : `[lr, lt, lp]` smoothing coefficients for laplacian smoothing on r,t and p direction.
- `optim_method` : `0`, `1` or `2`. 0 for gradient descent, 1 for gradient descent with halve-stepping, 2 for L-BFGS (THIS MODE IS EXPERIMENTAL AND CURRENTLY NOT WORKING.).
- `regularization_weight`: regularization weight, USED ONLY L-BFGS mode.
- `max_sub_iterations` : maximum number of sub iteration. Used for optim_method 1 and 2.

#### parallel :
- `n_sims` : number of simulutaneous run
- `ndiv_rtp` : `[ndiv_r, ndiv_t, ndiv_p]`  number of domain decomposition on r, t and p direction
- `nproc_sub` : number of processes used for sweeping parallelization
- `use_gpu` :  `0` or `1`. Set 1 for using GPU. (Currently gpu mode is used only for a forward simulation. n_sims, ndiv_rtp and nproc_sub need to be 1.)

The total number of mpi processes (i.e. mpirun -n NUMBER) must be n_sims\*ndiv_r\*ndiv_t\*ndiv_p\*nproc_sub, otherwise the code will stop instantly.

#### calculation :
- `convergence_tolerance` : convergence criterion for forward and adjoint run.
- `max_iterations` : number of maximum iteration for forward and adjoint run
- `stencil_order` :  `1` or `3`. The order of stencil for sweeping.
- `sweep_type` :  `0`or `1`. 0 is for sweeping in legacy order (threefold loops on r,t and p), 1 for cuthill-mckee node ordering with sweeping parallelization.
- `output_file_format` : `0` or `1` for selecting input and output file format. `0` is for HDF5 format, `1` is for ASCII format.


### 2. source receiver file
Source receiver file is a file which defines source and receiver positions and arrival times.

Below is an example:
```
1   1992  1  1  2 43  56.900    1.8000     98.9000 137.00  2.80    8    305644 <- source line   　
1      1      PCBI       1.8900     98.9253   1000.0000  P   10.40  18.000     <- receiver 1
1      2      MRPI       1.6125     99.3172   1100.0000  P   50.84  19.400     <- receiver 2
1      3      HUTI       2.3153     98.9711   1600.0000  P   57.84  19.200
     ....

```

```
Source line : id_src year month day hour min sec lat lon dep_km magnitude num_recs id_event (weight)
Receiver line : id_src id_rec name_rec lat lon elevation_m phase epicentral_distance_km arrival_time_sec (weight)
```

`num_recs` (number of receivers for this source) need to be the same with the number of receiver lines.
The last column of both source and receiver line is for put weight (on objective function). These are the optional column and set to be 1.0 if not written.
`name_rec` need to be different for each station.
`lon` and `lat` are in degree.

If the source position is out of the global domain (defined in input parameter file.), the code will flag this event as a teleseismic event and run the dedicated routine for teleseimic event. For teleseismic event, swap_src_rec will be ignored for this event (as the teleseismic case, a source is not a point but boundary surfaces).


### 3. initial model file

Initial model file is used for defining parameters of input mode.
Necessary parameters are `fun` (slowness), `eta`, `xi`, `zeta`, `fac_a`, `fac_b`, `fac_c`, `fac_f`.

#### inital model file in HDF5 format

In HDF5 I/O mode (`output_file_format`: 0 in input parameter file), all the necessary parameters should be saved in one single `.h5` file, with exactly the same name of dataset with parameter name written above.
The dimension of dataset should be the same with `n_rtp` in input parameter file.
Please refer the `examples/inversion_small/make_test_model.py` for details.


#### initial model file in ASCII format

In ASCII I/O mode (`output_file_format`: 1 in input parameter file), all the necessary parameters should be save in one single ASCII file.
The number of rows in the file need to be equivalent with the number of global nodes (i.e. n_rtp[0]\*n_rtp[1]\*n_rtp[2]).

The node order should be:
```python
    # write nodes in rtp
    for ir in range(n_rtp[0]): # number of nodes on r direction
        for it in range(n_rtp[1]): # number of nodes on t direction
            for ip in range(n_rtp[2]): # number of nodes on p direction
                # write out parameters
                # eta xi zeta fun fac_a fac_b fac_c fac_f

```

Complete example may be found `examples/inversion_small_ASCII/make_test_model.py`.

---

## Output files

Calculated travel times at the stations will be writen in `(source receiver file)_out.dat` on the column for travel times.

Volumetric result data files are saved in OUTPUT_FILES directory.

As the node order in the output file is not in the global domain but each subdomains, it is necessary to do a small post-processing for extracting slices.
`utils/tomoatt_data_retrieval.py` includes functions for this post processing.
Please refer the concrete example in `inversion_small/data_post_process.py` for HDF5 mode, and `inversion_small_ASCII/data_post_process.py` for ASCII mode.


### HDF5 I/O mode
In HDF5 mode, the code will carry out collective writing from all MPI processes into one single output file, which will try to maximize the I/O bandwidth for efficient I/O.

TomoATT produces output files like (for i-th event):
- out_data_grid.h5 : grid coordinate and connectivity data
- out_data_sim_i.h5 : field data
- out_data_sim_i.xmf : XDMF index data for visualizing 3D data. This may be open by Paraview.

Travel time field for i-th source may be visualized by reading `OUTPUT_FILES/out_data_sim_i.xmf`.
All the inversed parameters from all the sources and receivers are saved in `out_data_sim_0.xmf`.


Internal composition of .h5 data may be indicated by `h5ls -r` command.
The composition of out_data_grid.h5 is :
```
/                        Group
/Mesh                    Group
/Mesh/elem_conn          Dataset {21609, 9}  # node connectivity used for visualization
/Mesh/node_coords_p      Dataset {26010}     # longitude [degree]
/Mesh/node_coords_r      Dataset {26010}     # radious [km]
/Mesh/node_coords_t      Dataset {26010}     # latiude [degree]
/Mesh/node_coords_x      Dataset {26010}     # xyz coordinate in WGS84
/Mesh/node_coords_y      Dataset {26010}
/Mesh/node_coords_z      Dataset {26010}
/Mesh/procid             Dataset {26010}     # mpi processor id
```


out_data_sim_i.h5 is :
```
/                          Group
/Data                      Group
/Data/Keta_inv_000j        Dataset {26010}         # Kernel eta of j-th inversion iteration     
/Data/Keta_update_inv_000j Dataset {26010}         # Smoothed kernel eta of j-th inversion iteration
/Data/Ks_inv_000j          Dataset {26010}         # Kernel s(lowness)
/Data/Ks_update_inv_000j   Dataset {26010}         # Smoothed Kernel s
/Data/Kxi_inv_000j         Dataset {26010}         # Kernel xi
/Data/Kxi_update_inv_000j  Dataset {26010}         # Smoothed kernel xi
/Data/T0v_src_i_inv_000j   Dataset {26010}         # initial travel time field for i-th source
/Data/T_res_src_i_inv_000j Dataset {26010}         # calculated travel time field
/Data/adjoint_field_src_0_inv_000j Dataset {26010} # adjoint field
/Data/eta_inv_000j         Dataset {26010}         # eta
/Data/fac_b_inv_000j       Dataset {26010}         # factor b
/Data/fac_c_inv_000j       Dataset {26010}         # factor c
/Data/fac_f_inv_000j       Dataset {26010}         # factor f
/Data/fun_inv_000j         Dataset {26010}         # fun (slowness)
/Data/tau_src_0_inv_000j   Dataset {26010}         # tau (T0*tau = T_res)
/Data/xi_inv_000j          Dataset {26010}         # xi
```


### ASCII I/O mode
In ASCII mode, code will be do independent writing, (i.e. each MPI process do I/O process sequencially) in a single output file.

The files that TomoATT creates in OUPUT_FILES directory are :

```
adjoint_field_inv_000j_src_000i.dat  # adjoint field for i-th src and j-th inversion
eta_inv_000j.dat                     # eta
fac_b_inv_000j.dat                   # factor b
fac_c_inv_000j.dat                   # factor c
fac_f_inv_000j.dat                   # factor f
fun_inv_000j.dat                     # fun (slowness)
Keta_inv_000j.dat                    # Kernel eta
Keta_update_inv_000j.dat             # Smoothed kernel eta
Ks_inv_000j.dat                      # Kernel s(lownless)
Ks_update_inv_000j.dat               # Smoothed Kernel s
Kxi_inv_000j.dat                     # Kernel xi
Kxi_update_inv_000j.dat              # Smoothed kernel xi
out_grid_conn.dat                    # node connectivity
out_grid_ptr.dat                     # grid coordinate lon (degree), lat (degree), radius (km)
out_grid_xyz.dat                     # grid coordinate in WGS84
T0v_inv_000j_src_000i.dat            # inital travel time field
tau_inv_000j_src_000i.dat            # tau (T0*tau = T_res)
T_res_inv_000j_src_000i.dat          # calculated travel time field
xi_inv_000j.dat                      # xi 
```

