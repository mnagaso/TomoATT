# target directory to be run
TARGET_DIR=$1

# option to run only cpp code
CPP_ONLY=$2

if [ -z $TARGET_DIR ]; then
    echo "Usage: $0 <target_dir>"
    exit 1
fi

if [ ! -d $TARGET_DIR ]; then
    echo "Target directory $TARGET_DIR does not exist"
    exit 1
fi

# cd to fortran_code directory
cd $TARGET_DIR/fortran_code

# skip here if  not cpp only
if [ -z $CPP_ONLY ]; then
   # compile ega5_sphe_3d_kernel.f90
    mpif90 -O3 ega5_sphe_3d_kernel.f90
    # mkdir ega5/output directory
    mkdir ega5/output
    # mkdir ega5/output/tele_traveltime_field directory
    mkdir ega5/output/tele_traveltime_field
    # run ega5_sphe_3d_kernel.f90
    ./a.out
fi

# cd to one directory above
cd ../

# run make_test_model.py
python make_test_model.py
# run synthetic test with nproc 8
mpirun --oversubscribe -np 8 ../../build/TOMOATT -i ./input_params_pre.yml
# run inverse
mpirun --oversubscribe -np 8 ../../build/TOMOATT -i ./input_params.yml


