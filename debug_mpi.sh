#!/bin/sh

# general mpi debugging

# instantly start the debugger
#mpirun --oversubscribe -np $1 xterm -e gdb  -ex run --args ../../build/bin/TOMOATT -v -i $2
# for break point insertion
mpirun --oversubscribe -np $1 xterm -e gdb  --args ../../build/bin/TOMOATT -v -i $2
# tui mode
#mpirun --oversubscribe -np $1 xterm -e gdb --tui --args ../../build/bin/TOMOATT -i $2


# valgrind memory leak check
#mpirun --oversubscribe -np $1 valgrind --log-file="log_val" --leak-check=yes --track-origins=yes ../../build/bin/TOMOATT -i $2

# cuda debug
#mpirun --oversubscribe -np $1 xterm -e cuda-gdb  -ex run --args ../../build/bin/TOMOATT $2
# nvprof
#nvprof mpirun --oversubscribe -np 1 ../../build/bin/TOMOATT $2
