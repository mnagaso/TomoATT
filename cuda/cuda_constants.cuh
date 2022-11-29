#ifndef CUDA_CONSTANTS_H
#define CUDA_CONSTANTS_H

#include <cuda_runtime.h>
#include <cuda.h>

#define CUSTOMREAL double // need here for cuda kernels
#define MPI_CR     MPI_DOUBLE
//#define CUSTOMREAL float // need here for cuda kernels
//#define MPI_CR     MPI_FLOAT

#define MPI_DUMMY_TAG_CUDA 9999

// maximum grid dimension in one direction of GPU
//#define MAXIMUM_GRID_DIM 65535

#define CUDA_MAX_BLOCK_SIZE 1024
#define CUDA_MAX_GRID_SIZE 65535
#define CUDA_MAX_THREADS_PER_BLOCK 1024

//#define CUDA_SWEEPING_BLOCK_SIZE 16 // s
//#define CUDA_SWEEPING_BLOCK_SIZE 32 // 15.254 s
//#define CUDA_SWEEPING_BLOCK_SIZE 64 // 15.281 s
//#define CUDA_SWEEPING_BLOCK_SIZE 128 // 15.378 s
//#define CUDA_SWEEPING_BLOCK_SIZE 256 //　s
#define CUDA_SWEEPING_BLOCK_SIZE 512 //
//#define CUDA_SWEEPING_BLOCK_SIZE 1024 //


#define CUDA_L1_BLOCK_SIZE 128
//#define CUDA_L1_BLOCK_SIZE 256

#define CUDA_MAX_NUM_STREAMS 32

#endif // CUDA_CONSTANTS_H