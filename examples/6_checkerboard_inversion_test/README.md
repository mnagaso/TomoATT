# checkerboard inversion test 

This is an example of checkerboard test to invert for Vp and anisotropy

1. this example use the model files, src_rec_file and input_params.yml in `0_generate_files_for_TomoATT`
    - `0_generate_files_for_TomoATT/2_models/model_init_N61_61_61.h5`
    - `0_generate_files_for_TomoATT/2_models/model_ckb_N61_61_61.h5`
    - `0_generate_files_for_TomoATT/1_src_rec_files/src_rec_config.dat`
    - `0_generate_files_for_TomoATT/3_input_params/0_input_params_forward_simulation.yaml`
    - `0_generate_files_for_TomoATT/3_input_params/1_input_params_inversion.yaml`

You can check the distribution of earthquakes (star) and stations (triangle) in `0_generate_files_for_TomoATT/img/src_rec.jpg`

![](../0_generate_files_for_TomoATT/img/src_rec.jpg)

2. run all cells of "1_generate_input_params.ipynb" to generate the necessary input_params files:
    - `input_params/input_params_signal.yaml`
    - `input_params/input_params_inv_abs.yaml`

3. run TOMOATT forward with `input_params/input_params_signal.yaml` to compute traveltime data in checkerboard model
``` bash
mpirun --oversubscribe -n 8 ../../build/bin/TOMOATT -i input_params/input_params_signal.yaml
```
4. run TOMOATT forward with `input_params/input_params_inv_abs.yaml` to invert for Vp and anisotropy using absolute traveltime data from the initial model
``` bash
mpirun --oversubscribe -n 8 ../../build/bin/TOMOATT -i input_params/input_params_inv_abs.yaml
```
5. finally, you can run all cells of `3_plot_ckb_model.ipynb` to plot the checkerboard model 
``` bash
python  plot_ckb_model.py
```

The checkerboard model:

![](img/ckb_model_vel.jpg)

![](img/ckb_model_ani.jpg)

run all cells of `4_plot_inversion_result.ipynb` to plot the inversion result 

![](img/OUTPUT_FILES_inv_abs_0040_vel.jpg)

![](img/OUTPUT_FILES_inv_abs_0040_ani.jpg)

Users can run `bash run_this_example.sh` to proceed above Steps.


