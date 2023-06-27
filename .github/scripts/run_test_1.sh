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


# test example
cd $dir

if [ "$TESTDIR" == "examples/1_forward_accuracy_1st_upwind_isotropic/" ]; then
    echo "Running 1_forward_accuracy_1st_upwind_isotropic example..."
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

# clean up
rm -rf OUTPUT_FILES* src_rec_true.dat time.txt *h5 models

echo
echo "test finished successfully!"
echo `date`
echo
