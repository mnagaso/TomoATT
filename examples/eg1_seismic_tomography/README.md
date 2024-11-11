# seismic tomography

This is a toy model to invert traveltimes for Vp and anisotropy (Paper XXX, Figure 8c.)

Reference:
citation

Python modules are required to initiate the inversion and to plot final results:
- h5py
- PyTomoAT
- Pygmt
- gmt

Run this example:

1. Run bash script `bash run_this_example.sh` to execute the test.

2. After inversion, run `plot_output.py` to plot the results.

The target model:

![](img/ckb_model.jpg)

The inversion result:

![](img/OUTPUT_FILES_inv_abs_0040_ani.jpg)


