mkdir OUTPUT_FILES

# update model parameters and relocation
mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params_ckb_inv_error.yml

# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params_ckb_inv_good_1.yml

# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params_ckb_inv_good_2.yml
