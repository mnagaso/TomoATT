# Common receiver (swapped common source double difference) example

1. run make_test_model.ipynb to generate a test model
2. calculate the true travel times by running   
`mpirun --oversubscribe -n 8 ../../build/bin/TOMOATT -i input_params_pre.yml`
3. generate one receiver pair by running gen_rec_pairs.ipynb
4. run the inversion by running  
`mpirun --oversubscribe -n 8 ../../build/bin/TOMOATT -i input_params.yml`