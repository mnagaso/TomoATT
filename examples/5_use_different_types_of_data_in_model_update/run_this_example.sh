#!/bin/bash

python make_test_model.py

mkdir OUTPUT_FILES

echo "run for src rec no swapping cases"
# run for no swap abs_cr_cs
mpirun -n 2 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_no_swap_abs_cr_cs.yml

# run for no swap abs_cr
mpirun -n 2 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_no_swap_abs_cr.yml

# run for no swap abs_cs
mpirun -n 2 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_no_swap_abs_cs.yml

# run for no swap abs
mpirun -n 2 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_no_swap_abs.yml

# run for no swap cr_cs
mpirun -n 2 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_no_swap_cr_cs.yml

# run for no swap cr
mpirun -n 2 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_no_swap_cr.yml

# run for no swap cs
mpirun -n 2 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_no_swap_cs.yml

# run for no swap no data
mpirun -n 2 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_no_swap_no_data.yml


echo "run for src rec swapping cases"
# run for no swap abs_cr_cs
mpirun -n 2 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_swap_abs_cr_cs.yml

# run for no swap abs_cr
mpirun -n 2 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_swap_abs_cr.yml

# run for no swap abs_cs
mpirun -n 2 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_swap_abs_cs.yml

# run for no swap abs
mpirun -n 2 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_swap_abs.yml

# run for no swap cr_cs
mpirun -n 2 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_swap_cr_cs.yml

# run for no swap cr
mpirun -n 2 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_swap_cr.yml

# run for no swap cs
mpirun -n 2 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_swap_cs.yml

# run for no swap no data
mpirun -n 2 --allow-run-as-root ../../build/bin/TOMOATT -i input_params/input_params_swap_no_data.yml
