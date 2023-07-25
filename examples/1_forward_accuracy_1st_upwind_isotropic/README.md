# Forward modeling test 

This is an example for checking the accuracy of the 1st order eikonal solver with upwind scheme in isotropic media (Chen et al., 2023, GJI)

![](img/Chen_et_al_2023_GJI.png)

1. Run all cells of `make_test_model.ipynb` or python script `make_test_mode.py` for creating necessary input files: 
    - source, receiver file (src_rec_true.dat)
    - model with different mesh (model_N%d_%d_%d.h5)


2. then run TOMOATT forward with `input_params/input_params_N%d.yml` for calculating the numerical arrival times in different mesh
Following command will run the forward simulation with one mesh for example
``` bash
mpirun --oversubscribe -n 2 ../../build/bin/TOMOATT -i input_params/input_params_N41.yml
```
The computed src_rec_file is saved as `OUTPUT_FILES_N41/src_rec_file_forward.dat`.

3. compare the true traveltime with the numerical solution by running all cells of `compare_src_rec.ipynb` or python script `compare_src_rec.py`

You will obtain the numerical error in different mesh and the convergence order of the eikonal solver.

You can run `bash run_this_example.sh` to proceed Step 1-3. Result is correct is you have the following output:

The numerical error in a mesh of 41^3 is:  0.113 s (mean),  0.085 s (std),  0.382 (max)
The numerical error in a mesh of 61^3 is:  0.075 s (mean),  0.059 s (std),  0.268 (max), order of accuracy is  1.01
The numerical error in a mesh of 81^3 is:  0.056 s (mean),  0.046 s (std),  0.206 (max), order of accuracy is  0.99
The numerical error in a mesh of 121^3 is:  0.038 s (mean),  0.032 s (std),  0.143 (max), order of accuracy is  0.99
The numerical error in a mesh of 161^3 is:  0.028 s (mean),  0.025 s (std),  0.110 (max), order of accuracy is  1.00

(PS: order of accuracy is close to 1, indicating the eikonal solver is correct and has a firts order of accuracy)

