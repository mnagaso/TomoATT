#!/bin/bash

# python make_test_model.py
# run make_test_model.ipynb


mkdir OUTPUT_FILES

# compute synthetic traveltime in the ckb model
# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_signal.yaml

# do ckb inversion using abs data
mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_inv_abs.yaml

# # do ckb inversion using abs common receiver double difference 
# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_inv_cr.yml

# # do ckb inversion using abs common source double difference
# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_inv_cs.yml

# # do ckb inversion using abs common source double difference
# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_inv_abs_cs.yml

# # do ckb inversion using abs common source double difference
# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_inv_abs_cr.yml

# # do ckb inversion using abs common source double difference
# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_inv_cs_cr.yml

# do ckb inversion using abs + cr + cs
# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_inv_abs_cs_cr.yml


# run model_visualization.ipynb to plot ckb model and inv model
