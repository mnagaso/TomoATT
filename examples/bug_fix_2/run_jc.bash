mkdir OUTPUT_FILES

# update model parameters and relocation
mpirun -n 4 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params_jc.yml
