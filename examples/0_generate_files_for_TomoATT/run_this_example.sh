#!/bin/bash

python 1_build_initial_model.py
# python 2_build_checkerboard_model.py
sh 2_build_checkerboard_model.sh
python 3_generate_src_rec_file.py