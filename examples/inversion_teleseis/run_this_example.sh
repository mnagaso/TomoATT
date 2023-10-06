# make test model source receiver file
python make_test_model.py

TOMOATT=../../build/bin/TOMOATT

# pre run
mpirun -n 4 $TOMOATT -i input_params_pre.yml

# run iversio
mpirun --oversubscribe -n 8 $TOMOATT -i input_params.yml