# TomoATT User Manual version 1.0

## Introduction

TomoATT is a library which implements an eikonal equation solver based Adjoint-state Travel-Time Tomography for a very large-scale computation, following a published article [Ping Tong (2021)](https://doi.org/10.1029/2021JB021818) and [Jing Chen (2022)](add_here_when_published).
The slowness field and anisotropy field are computed in spherical coordinate system.

Thanks to the efficiency of an eikonal equation solver, the computation of the travel-time is very fast and requires less amount of computational resources.
As an input data for TomoATT is travel times at seismic stations, we can easily prepare a great amount of data for the computation.

This library aims to be used for modeling a very-large domain. For this purpose, 3-layer parallelization is applied, which are:
- layer 1: simulutaneous run parallelization (travel times for multiple seismic sources may be calculated simultaneously)
- layer 2: subdomain decomposition (If the number of computational nodes requires too large memory, we can separate the domain into subdomains and run each subdomain in a separate compute node)
- layer 3: sweeping parallelization (in each subdomain, sweeping layers are also parallelized)

The details of the parallelization method are described in the paper [Miles Detrixhe and Frédéric Gibou (2016)](https://doi.org/10.1016/j.jcp.2016.06.023).

Regional events (sources within the global domain) and teleseismic events (sources outside the global domain) may be used for inversion.



## Input files

TomoATT requires 3 input files : 
- input parameter file (setup file for simulation parameters)
- source receiver file (source and receiver definitions and observation arrival times)
- initial model file (3d initial model)


### input parameter file

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
  do_inversion : 0 # 0 for forward simulation only, 1 for inversion

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

There are 6 category tags and some parameters under each category.
Belows are all the parameters which can be set.

#### domain
The domain category is for setting a global domain.
- `min_max_dep` :  minimum and maximum depth of the domain in kilo meter
- `min_max_lat` :  minimum and maximum latitude in degree
- `min_max_lon` :  minimum and maximum longitude in degree
- `n_rtp` : number of computation nodes for r (depth),t (latitude) and p (longitude) direction

#### source
- `src_rec_file` : for specifying a path to source receiver file (details of this file is in the following section). The calculated travel time at receivers will output as a new source receiver file. If src_rec_file : `src_rec.dat`, TomoATT will output a file `src_rec_out.dat`.
- `swap_src_rec` : `0` or `1`. Set 1 for using receivers as sources, which is useful when the number of sources are larger than the number of receivers. (quite usual when using large dataset.)

#### model
- `init_model_path` : File path for initial model file. Details will be explained in the following section.

#### inversion
- `do_inversion` : `0` for running only a forward simulation. `1` for do inversion.
- `n_inversion_grid` :  the number of inversion grid.
- `n_inv_rtp` : the numbers of inversion grids for r, t and p direction/
- `min_max_dep_inv` :  minimum and maximum depth of inversion grid in kilo meter.
- `min_max_lat_inv` :  minimum and maximum latitude of inversion grid in degree.
- `min_max_lon_inv` :  minimum and maximum longitude of inversion grid in degree.
- `max_iterations_inv` :  The limit of iteration number for inversion.
- `step_size` :  Maximum step size ratio for updating model.
- `smooth_method` : `0` or `1`. 0 for multigrid parametrization ([Ping Tong 2019](https://doi.org/10.1093/gji/ggz151)). 1 for laplacian smoothing with CG method.
- `l_smooth_rtp` :  smoothing coefficients for laplacian smoothing on r,t and p direction.
- `optim_method` : `0`, `1` or `2`. 0 for gradient descent, 1 for gradient descent with halve-stepping, 2 for L-BFGS (THIS MODE IS EXPERIMENTAL AND CURRENTLY NOT WORKING.).
- `regularization_weight`: regularization weight, USED ONLY L-BFGS mode.
- `max_sub_iterations` : maximum number of sub iteration. Used for optim_method 1 and 2.

#### parallel
- `n_sims` :  number of simulutaneous run
- `ndiv_rtp` :  number of domain decomposition on r, t and p direction
- `nproc_sub` :  number of processes used for sweeping parallelization
- `use_gpu` :  `0` or `1`. Set 1 for using GPU. (Currently gpu mode is used only for a forward simulation. n_sims, ndiv_rtp and nproc_sub need to be 1.)

The total number of mpi processes (i.e. mpirun -n NUMBER) must be n_sims\*ndiv_r\*ndiv_t\*ndiv_p\*nproc_sub, otherwise the code will stop instantly.

#### calculation
- `convergence_tolerance` :  convergence criterion for forward and adjoint run.
- `max_iterations` :  number of maximum iteration for forward and adjoint run
- `stencil_order` :  `1` or `3`. The order of stencil for sweeping.
- `sweep_type` :  `0`, `1` or `2`. 0 is for sweeping in legacy order (threefold loops on r,t and p), 1 for cuthill-mckee ordering without parallelization (for debuging purpose), 2 for cuthill-mckee node ordering with sweeping parallelization.
- `output_file_format` : `0` or `1` for selecting input and output file format. `0` is for HDF5 format, `1` is for ASCII format.


### source receiver file
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


### initial model file

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


## Output files

Calculated travel times at the stations will be writen in `(source receiver file)_out.dat` on the column for travel times.

Volumetric result data files are saved in OUTPUT_FILES directory.
As the node order in the output file is not in the global domain but each subdomains, it is necessary to do a small post-processing for extracting slices.
`utils/tomoatt_data_retrieval.py` includes functions for this post processing.
Please refer the concrete example in `inversion_small/data_post_process.py` for HDF5 mode, and `inversion_small_ASCII/data_post_process.py` for ASCII mode.

In HDF5 mode, the code will carry out collective writing from all MPI processes into one single output file, which will try to maximize the I/O bandwidth for efficient I/O.

In ASCII mode, code will be do independent writing, (i.e. each MPI process do I/O process sequencially) in a single output file.

Additianlly with HDF5 I/O mode, TomoATT creates xdmf index file for easily visualizing the result files.
Travel time field for i-th source may be visualized by reading `OUTPUT_FILES/out_data_sim_(i).xmf`.
All the inversed parameters from all the sources and receivers are saved in `out_data_sim_0.xmf`.
