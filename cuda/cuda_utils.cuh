#ifndef CUDA_UTILS_H
#define CUDA_UTILS_H

#include <mpi.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
#include <vector>

#include "cuda_constants.cuh"

// function to convert kernel,grid to ijk
#define I2V_cuda(i,j,k,II,JJ) ((k)*(JJ)*(II)+(j)*(II)+(i))


// allocate memory on device
cudaError_t allocate_memory_on_device_i(void** d_ptr, size_t size);
cudaError_t allocate_memory_on_device_cv(void** d_ptr, size_t size);
cudaError_t allocate_memory_on_device_bl(void** d_ptr, size_t size);
cudaError_t allocate_memory_on_device_cv_pinned(void** d_ptr, size_t size);


// deallocate memory on device
cudaError_t deallocate_memory_on_device_i(int *&d_ptr);
cudaError_t deallocate_memory_on_device_cv(CUSTOMREAL *&d_ptr);
cudaError_t deallocate_memory_on_device_bl(bool *&d_ptr);


// copy memory from host to device
cudaError_t copy_host_to_device_i(int *d_ptr, int *h_ptr, const size_t size);
cudaError_t copy_host_to_device_cv(CUSTOMREAL *d_ptr, CUSTOMREAL *h_ptr, const size_t size);
cudaError_t copy_host_to_device_bl(bool *d_ptr, bool *h_ptr, const size_t size);


// copy memory from device to host
cudaError_t copy_device_to_host_i(int *h_ptr, int *d_ptr, size_t size);
cudaError_t copy_device_to_host_cv(CUSTOMREAL *h_ptr, CUSTOMREAL *d_ptr, size_t size);

// allocate and copy to device
void* allocate_and_copy_host_to_device_i(int* h_ptr, size_t size, int num);
void* allocate_and_copy_host_to_device_cv(CUSTOMREAL* h_ptr, size_t size, int num);

// allocate, flatten and copy to device
void flatten_arr_i(int* h_ptr_flattened, std::vector<int*> &h_v, int size_total, int* size_each);
void flatten_arr_cv(CUSTOMREAL* h_ptr_flattened, std::vector<CUSTOMREAL*> &h_v, int size_total, int* size_each);

void* allocate_and_copy_host_to_device_flattened_i(std::vector<int*>&vh, int size_total, int* size_each, int num);
void* allocate_and_copy_host_to_device_flattened_cv(std::vector<CUSTOMREAL*>& vh, int size_total, int* size_each, int num);


// mpi send recv
static inline void cuda_send_cr(CUSTOMREAL* buf, int count, int dest, MPI_Comm inter_sub_comm){
    MPI_Send(buf, count, MPI_CR, dest, MPI_DUMMY_TAG_CUDA, inter_sub_comm);
}

static inline void cuda_recv_cr(CUSTOMREAL* buf, int count, int source, MPI_Comm inter_sub_comm){
    MPI_Status stat;
    MPI_Recv(buf, count, MPI_CR, source, MPI_DUMMY_TAG_CUDA, inter_sub_comm, &stat);
}

static inline void cuda_synchronize_all_sub(MPI_Comm& sub_comm){
    MPI_Barrier(sub_comm);
}

inline void cuda_wait_req(MPI_Request& req){
    MPI_Status status;
    MPI_Wait(&req, &status);
}



static inline void get_block_xy(int num_blocks, int* num_blocks_x, int* num_blocks_y) {
    // at first, the num_blocks_x is set with equal value of num_blocks, and num_blocks_y is set with value 1
    // when the num_blocks_x exceeds the block size limit of 65535, the num_blocks_x is divided by 2 and num_blocks_y is increased by 1
    *num_blocks_x = num_blocks;
    *num_blocks_y = 1;

    while (*num_blocks_x > CUDA_MAX_GRID_SIZE) {
        *num_blocks_x = (int) ceil(*num_blocks_x * 0.5f);
        *num_blocks_y = *num_blocks_y * 2;;
    }

}


static inline void get_thread_block_for_3d_loop(int nx, int ny, int nz, dim3* threads, dim3* blocks)  {
    threads->x = 8; threads->y = 8; threads->z = 8; // use 512 threads in total
    blocks->x = (nx + threads->x - 1)/threads->x;
    blocks->y = (ny + threads->y - 1)/threads->y;
    blocks->z = (nz + threads->z - 1)/threads->z;
}


static inline void get_thread_block_for_ibound(int nx, int ny, int nz, dim3* threads, dim3* blocks)  {
    threads->x = nx; threads->y = 8; threads->z = 8;
    blocks->x = (nx + threads->x - 1)/threads->x;
    blocks->y = (ny + threads->y - 1)/threads->y;
    blocks->z = (nz + threads->z - 1)/threads->z;
}


static inline void get_thread_block_for_jbound(int nx, int ny, int nz, dim3* threads, dim3* blocks)  {
    threads->x = 8; threads->y = ny; threads->z = 8;
    blocks->x = (nx + threads->x - 1)/threads->x;
    blocks->y = (ny + threads->y - 1)/threads->y;
    blocks->z = (nz + threads->z - 1)/threads->z;
}


static inline void get_thread_block_for_kbound(int nx, int ny, int nz, dim3* threads, dim3* blocks)  {
    threads->x = 8; threads->y = 8; threads->z = nz;
    blocks->x = (nx + threads->x - 1)/threads->x;
    blocks->y = (ny + threads->y - 1)/threads->y;
    blocks->z = (nz + threads->z - 1)/threads->z;
}


inline void cuda_isend_cr(CUSTOMREAL* buf, int count, int dest, MPI_Comm& comm, MPI_Request& request){
    //MPI_Request request = MPI_REQUEST_NULL;
    //std::cout << "sending from : " << inter_sub_rank << ", to : " << dest <<", size : " << count << std::endl;
    int DUMMY_TAG = 9999;
    MPI_Isend(buf, count, MPI_CR, dest, DUMMY_TAG, comm, &request);
}

inline void cuda_irecv_cr(CUSTOMREAL* buf, int count, int source, MPI_Comm& comm, MPI_Request& request){
    //MPI_Request request = MPI_REQUEST_NULL;
    //std::cout << "receiving by : " << inter_sub_rank << ", from : " << source << ", size : " << count << std::endl;
    int DUMMY_TAG = 9999;
    MPI_Irecv(buf, count, MPI_CR, source, DUMMY_TAG, comm, &request);
}


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

inline void print_memory_usage(){
    size_t free_byte ;
    size_t total_byte ;
    cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;
    if ( cudaSuccess != cuda_status ){
        printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
        exit(1);
    }

    double free_db = (double)free_byte ;
    double total_db = (double)total_byte ;
    double used_db = total_db - free_db ;

    printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
            used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);
}

inline void print_CUDA_error_if_any(cudaError_t err, int num) {
  if (cudaSuccess != err)
  {
    printf("\nCUDA error !!!!! <%s> !!!!! \nat CUDA call error code: # %d\n",cudaGetErrorString(err),num);
    fflush(stdout);

    // outputs error file
    FILE* fp;
    int myrank;
    char filename[BUFSIZ];
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    sprintf(filename,"error_message_%06d.txt",myrank);
    fp = fopen(filename,"a+");
    if (fp != NULL) {
      fprintf(fp,"\nCUDA error !!!!! <%s> !!!!! \nat CUDA call error code: # %d\n",cudaGetErrorString(err),num);
      fclose(fp);
    }

    // check memory usage
    size_t free_byte ;
    size_t total_byte ;
    cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;

    if ( cudaSuccess != cuda_status ){
      printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
      fflush(stdout);
      exit(1);
    }

    // print usage
    double free_db = (double)free_byte ;
    double total_db = (double)total_byte ;
    double used_db = total_db - free_db ;
    printf("GPU memory usage: used = %f, free = %f MB, total = %f MB", used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);


    // stops program
    //MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
}






#endif // CUDA_UTILS_H