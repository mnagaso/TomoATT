
# generate the model data and src_rec_file
python make_test_model.py

# run forward simulation
mpirun --oversubscribe -np 8 ../../build/bin/TOMOATT -i input_params_pre.yml

# modify the source position
python modify_source_location.py

# run forward simulation in earthquake relocation mode
mpirun --oversubscribe -np 8 ../../build/bin/TOMOATT -i input_params.yml

# result can be indicated in check_relocate_sourece.py/ipynb
python check_relocated_sourece.py