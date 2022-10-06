# target directory to be cleaned
TARGET_DIR=$1

if [ -z $TARGET_DIR ]; then
    echo "Usage: $0 <target_dir>"
    exit 1
fi

if [ ! -d $TARGET_DIR ]; then
    echo "Target directory $TARGET_DIR does not exist" exit 1
fi

# .h5 files
find $TARGET_DIR -name "*.h5" -exec rm {} \;
# erase test_model_*
find $TARGET_DIR -name "test_model_*" -exec rm {} \;
# erase OUTPUT_FILES directory
find $TARGET_DIR -name "OUTPUT_FILES" -exec rm -rf {} \;
# erase src_rec_test*
find $TARGET_DIR -name "src_rec_test*" -exec rm {} \;
# erase objective_function.txt
find $TARGET_DIR -name "objective_function.txt" -exec rm {} \;
# erase files in fortran_code/ega5/output/
find $TARGET_DIR/fortran_code/ega5/output -name "*" -exec rm {} \;
# erase files in fortran_code/a.out
find $TARGET_DIR/fortran_code/a.out -name "*" -exec rm {} \;
# erase files in fortran_code/ega5/output/tele_traveltime_field
find $TARGET_DIR/fortran_code/ega5/output/tele_traveltime_field -name "*" -exec rm {} \;
# erase cuda_device_info.txt
find $TARGET_DIR -name "cuda_device_info.txt" -exec rm {} \;
# erase error_message_*.txt
find $TARGET_DIR -name "error_message_*.txt" -exec rm {} \;
# erase time*.txt
find $TARGET_DIR -name "time*.txt" -exec rm {} \;
