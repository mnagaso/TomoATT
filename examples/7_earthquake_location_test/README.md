# Forward modeling test 

This is an example to test earthquake location in the checkerboard model using absolute traveltime and common receiver differential traveltime data.

1. Run all cells of `make_test_model.ipynb` or python script `make_test_mode.py` for creating necessary input files: 
    - source, receiver file (src_rec_config.dat)
    - model with a mesh (models/model_init_N61_61_61.h5 and models/model_ckb_N61_61_61.h5)

You can check the distribution of earthquakes (star) and stations (triangle) in img/src_rec.jpg

![](img/src_rec.jpg)

2. then run TOMOATT forward with `input_params/input_params_signal.yml` to compute traveltime data in checkerboard model
``` bash
mpirun --oversubscribe -n 2 ../../build/bin/TOMOATT -i input_params/input_params_signal.yml
```
3. then run all cells of `generate_syn_obs_src_rec_data.ipynb` or python script `generate_syn_obs_src_rec_data.py` for deviate source location and origin time.
``` bash
python  generate_syn_obs_src_rec_data.py
```

4. then run TOMOATT forward with `input_params/input_params_reloc_XX.yml` to do relocation using different types of data
``` bash
mpirun --oversubscribe -n 2 ../../build/bin/TOMOATT -i input_params/input_params_reloc_abs.yml
```
``` bash
mpirun --oversubscribe -n 2 ../../build/bin/TOMOATT -i input_params/input_params_reloc_cr.yml
```
``` bash
mpirun --oversubscribe -n 2 ../../build/bin/TOMOATT -i input_params/input_params_reloc_abs_cr.yml
```

3. finally, you can run `compare_location_result.ipynb` or `compare_location_result.py` to show the location result and evaluate the location error. 
``` bash
python  compare_location_result.py
```

You can get the location result as below in img

location using absolute traveltime data:
![](img/reloc_abs.jpg)

location using common receiver traveltime data:
![](img/reloc_cr.jpg)

location using two types of data:
![](img/reloc_abs_cr.jpg)


You can run `bash run_this_example.sh` to proceed all steps



