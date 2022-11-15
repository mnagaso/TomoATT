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


// copy memory from host to device
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
cudaError_t allocate_and_copy_host_to_device_i(int* d_ptr, int* h_ptr, size_t size)
{
    cudaError_t err0 = allocate_memory_on_device_i((void**)&d_ptr, size);
    cudaError_t err1 = copy_host_to_device_i(d_ptr, h_ptr, size);

    return err1;
}

cudaError_t allocate_and_copy_host_to_device_cv(CUSTOMREAL* d_ptr, CUSTOMREAL* h_ptr, size_t size)
{
    cudaError_t err0 = allocate_memory_on_device_cv((void**) &d_ptr, size);
    cudaError_t err1 = copy_host_to_device_cv(d_ptr, h_ptr, size);

    return err1;
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

cudaError_t allocate_and_copy_host_to_device_flattened_i(int* d_ptr, std::vector<int*>& vh, int size_total, int* size_each){
    // flatten
    int* h_ptr_flattened = new int[size_total];
    flatten_arr_i(h_ptr_flattened, vh, size_total, size_each);

    // allocate and copy
    cudaError_t err = allocate_and_copy_host_to_device_i(d_ptr, h_ptr_flattened, size_total);

    // free
    delete[] h_ptr_flattened;

    return err;
}

cudaError_t allocate_and_copy_host_to_device_flattened_cv(CUSTOMREAL* d_ptr, std::vector<CUSTOMREAL*>& vh, int size_total, int* size_each){
    // flatten
    CUSTOMREAL* h_ptr_flattened = new CUSTOMREAL[size_total];
    flatten_arr_cv(h_ptr_flattened, vh, size_total, size_each);

    // allocate and copy
    cudaError_t err = allocate_and_copy_host_to_device_cv(d_ptr, h_ptr_flattened, size_total);

    // free
    delete[] h_ptr_flattened;

    return err;
}

