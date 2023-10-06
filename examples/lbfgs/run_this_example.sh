#!/bin/bash

python make_test_model.py

# run for preparing true travel times
mpirun -n 2 ../../build/bin/TOMOATT -i ./input_params/input_params_pre.yml

NPROC=4
# run inversions for l-curve analysis
mpirun -n $NPROC ../../build/bin/TOMOATT -i ./input_params/input_params_l_0.00001.yml
mv OUTPUT_FILES/objective_function.txt OUTPUT_FILES/objective_function_l_0.00001.txt

mpirun -n $NPROC ../../build/bin/TOMOATT -i ./input_params/input_params_l_0.0001.yml
mv OUTPUT_FILES/objective_function.txt OUTPUT_FILES/objective_function_l_0.0001.txt

mpirun -n $NPROC ../../build/bin/TOMOATT -i ./input_params/input_params_l_0.001.yml
mv OUTPUT_FILES/objective_function.txt OUTPUT_FILES/objective_function_l_0.001.txt

mpirun -n $NPROC ../../build/bin/TOMOATT -i ./input_params/input_params_l_0.01.yml
mv OUTPUT_FILES/objective_function.txt OUTPUT_FILES/objective_function_l_0.01.txt

mpirun -n $NPROC ../../build/bin/TOMOATT -i ./input_params/input_params_l_0.1.yml
mv OUTPUT_FILES/objective_function.txt OUTPUT_FILES/objective_function_l_0.1.txt

mpirun -n $NPROC ../../build/bin/TOMOATT -i ./input_params/input_params_l_1.0.yml
mv OUTPUT_FILES/objective_function.txt OUTPUT_FILES/objective_function_l_1.0.txt


# read the objective function values and plot
python plot_l_curve.py
