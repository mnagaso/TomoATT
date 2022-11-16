#include "cuda_utils.cuh"



// allocate memory on device
cudaError_t allocate_memory_on_device_i(void** d_ptr, size_t size)
{
    return cudaMalloc((void**) d_ptr, size * sizeof(int));
}

cudaError_t allocate_memory_on_device_cv(void** d_ptr, size_t size)
{
    return cudaMalloc((void**) d_ptr, size * sizeof(CUSTOMREAL));
}

cudaError_t allocate_memory_on_device_bl(void** d_ptr, size_t size)
{
    return cudaMalloc((void**) d_ptr, size * sizeof(bool));
}


// device-host shared memory (pinned memory) (maybe unnecessary for CUDA-aware MPI)
cudaError_t allocate_memory_on_device_cv_pinned(void** d_ptr, size_t size)
{
    return cudaMallocHost((void**) d_ptr, size * sizeof(CUSTOMREAL));
}


// deallocate memory on device
cudaError_t deallocate_memory_on_device_i(int*& d_ptr)
{
    return cudaFree(d_ptr);
}

cudaError_t deallocate_memory_on_device_cv(CUSTOMREAL*& d_ptr)
{
    return cudaFree(d_ptr);
}

cudaError_t deallocate_memory_on_device_bl(bool*& d_ptr)
{
    return cudaFree(d_ptr);
}


// copy memory between host and device
cudaError_t copy_host_to_device_i(int* d_ptr, int* h_ptr, const size_t size)
{
    return cudaMemcpy(d_ptr, h_ptr, size * sizeof(int), cudaMemcpyHostToDevice);
}

cudaError_t copy_host_to_device_cv(CUSTOMREAL* d_ptr, CUSTOMREAL* h_ptr, const size_t size)
{
    return cudaMemcpy(d_ptr, h_ptr, size * sizeof(CUSTOMREAL), cudaMemcpyHostToDevice);
}

cudaError_t copy_host_to_device_bl(bool* d_ptr, bool* h_ptr, const size_t size)
{
    return cudaMemcpy(d_ptr, h_ptr, size * sizeof(bool), cudaMemcpyHostToDevice);
}

// copy memory from device to host
cudaError_t copy_device_to_host_i(int* h_ptr, int* d_ptr, size_t size)
{
    return cudaMemcpy(h_ptr, d_ptr, size * sizeof(int), cudaMemcpyDeviceToHost);
}
cudaError_t copy_device_to_host_cv(CUSTOMREAL* h_ptr, CUSTOMREAL* d_ptr, size_t size)
{
    return cudaMemcpy(h_ptr, d_ptr, size * sizeof(CUSTOMREAL), cudaMemcpyDeviceToHost);
}


// allocate and copy to device
void* allocate_and_copy_host_to_device_i(int* h_ptr, size_t size, int num)
{
    void* d_ptr;

    print_CUDA_error_if_any(allocate_memory_on_device_i(&d_ptr, size), num);
    print_CUDA_error_if_any(copy_host_to_device_i((int*)d_ptr, h_ptr, size),num);

    return d_ptr;
}

void* allocate_and_copy_host_to_device_cv(CUSTOMREAL* h_ptr, size_t size, int num)
{
    void* d_ptr;
    print_CUDA_error_if_any(allocate_memory_on_device_cv(&d_ptr, size),num);
    print_CUDA_error_if_any(copy_host_to_device_cv((CUSTOMREAL*)d_ptr, h_ptr, size), num);

    return d_ptr;
}

// allocate, flatten and copy from host to device
void flatten_arr_i(int* h_ptr_flattened, std::vector<int*>&h_v, int size_total, int* size_each)
{
    // flatten
    int counter = 0;
    int n_v = h_v.size();

    for (int i = 0; i < n_v; i++) { // levels
        for (int j = 0; j < size_each[i]; j++) {
            h_ptr_flattened[counter] = h_v.at(i)[j];
            counter++;
        }
    }
}

void flatten_arr_cv(CUSTOMREAL* h_ptr_flattened, std::vector<CUSTOMREAL*> &h_v, int size_total, int* size_each)
{
    // flatten
    int counter = 0;
    int n_v = h_v.size();

    for (int i = 0; i < n_v; i++) { // levels
        for (int j = 0; j < size_each[i]; j++) {
            h_ptr_flattened[counter] = h_v.at(i)[j];
            counter++;
        }
    }
}

void* allocate_and_copy_host_to_device_flattened_i(std::vector<int*>& vh, int size_total, int* size_each, int num){
    // flatten
    int* h_ptr_flattened = new int[size_total];
    flatten_arr_i(h_ptr_flattened, vh, size_total, size_each);

    // allocate and copy
    void* d_ptr = allocate_and_copy_host_to_device_i(h_ptr_flattened, size_total, num);

    // free
    delete[] h_ptr_flattened;

    return d_ptr;
}

void* allocate_and_copy_host_to_device_flattened_cv(std::vector<CUSTOMREAL*>& vh, int size_total, int* size_each, int num){
    // flatten
    CUSTOMREAL* h_ptr_flattened = new CUSTOMREAL[size_total];
    flatten_arr_cv(h_ptr_flattened, vh, size_total, size_each);

    // allocate and copy
    void* d_ptr = allocate_and_copy_host_to_device_cv(h_ptr_flattened, size_total, num);

    // free
    delete[] h_ptr_flattened;

    return d_ptr;
}

