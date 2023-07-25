# Forward modeling test 

This is an example for checking whether the code can read and output src_rec.dat file during the inversion.


1. Run all cells of `make_test_model.ipynb` or python script `make_test_mode.py` for creating necessary input files: 
    - source, receiver file (src_rec_obs.dat)
    - model with a mesh (model_N61_61_61.h5)


2. then run TOMOATT forward with `input_params/input_params_XXX.yml` for reading and output src_rec.dat files:
``` bash
mpirun --oversubscribe -n 2 ../../build/bin/TOMOATT -i input_params/input_params_no_swap_parallel.yml
```
The computed src_rec_file is saved as `OUTPUT_FILES_no_swap_parallel/src_rec_file_step0000.dat`.

3. compare the input src_rec.data with the output src_rec.data by running all cells of `compare_src_rec.ipynb` or python script `compare_src_rec.py`

You can run `bash run_this_example.sh` to proceed Step 1-3. Result is correct is you have the following output:

For no_swap_serial,   src_rec.dat read and output are CORRECT! (match coeffient: 0.9986460342912408)
For no_swap_parallel, src_rec.dat read and output are CORRECT! (match coeffient: 0.9986460342912408)
For swap_serial,      src_rec.dat read and output are CORRECT! (match coeffient: 0.9994013789951578)
For swap_parallel,    src_rec.dat read and output are CORRECT! (match coeffient: 0.9994013789951578)


