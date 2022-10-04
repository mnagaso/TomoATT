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

//#define CUDA_SWEEPING_BLOCK_SIZE 16 // 60.177s
//#define CUDA_SWEEPING_BLOCK_SIZE 32 // 33.103s
//#define CUDA_SWEEPING_BLOCK_SIZE 64 // 33.414s
#define CUDA_SWEEPING_BLOCK_SIZE 128 // 33.960s
//#define CUDA_SWEEPING_BLOCK_SIZE 256 //　35.937
//#define CUDA_SWEEPING_BLOCK_SIZE 512 //
//#define CUDA_SWEEPING_BLOCK_SIZE 1024 //


#define CUDA_L1_BLOCK_SIZE 128

// store device properties
//inline cudaDeviceProp deviceProp;

// do not use level-unified-kernel (currently set as true for strange block synchronization behavior)
#define FORCE_UNUSE_UNIFIED_KERNEL true

#endif // CUDA_CONSTANTS_H