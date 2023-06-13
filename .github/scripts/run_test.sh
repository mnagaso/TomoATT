#!/bin/bash
#
# runs a test example case
#

# getting updated environment
if [ -f $HOME/.tmprc ]; then
    source $HOME/.tmprc
fi

WORKDIR=`pwd`
dir=${TESTDIR}

# print info
echo "Running test case in a work directory: $WORKDIR"
echo `date`
echo
echo "******************************************************************"
echo
echo "test directory: $dir"
echo
echo "******************************************************************"
echo


# test script for checking the result
my_test(){
    echo "Checking the result..."
    ln -s ../../utils/compare_models.py
    python compare_models.py -t ./test_model_true.h5 -r ./OUTPUT_FILES/final_model.h5
    if [[ $? -ne 0 ]]; then
        echo "Test failed!"
        exit 1
    fi
    echo "Test passed!"
}

# test example
cd $dir

if [ "$TESTDIR" == "examples/inversion_small/" ]; then
    echo "Running inversion_small example..."
fi

# run the test
./run_this_example.sh

# check the exit code
if [[ $? -ne 0 ]]; then
    echo "Test failed!"
    exit 1
fi

# smulation finished, check the result
echo
echo "simulation finished: `pwd`"
echo `date`
echo

# check the result
my_test

# clean up
rm -rf OUTPUT_FILES* src_rec_test.dat time.txt *h5

echo
echo "test finished successfully!"
echo `date`
echo
