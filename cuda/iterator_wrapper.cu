#include "iterator_wrapper.cuh"

__device__ const CUSTOMREAL PLUS = 1.0f;
__device__ const CUSTOMREAL MINUS = -1.0f;
__device__ const CUSTOMREAL v_eps = 1e-12;

__device__ const CUSTOMREAL _0_5_CR   = 0.5f;
__device__ const CUSTOMREAL _1_CR     = 1.0f;
__device__ const CUSTOMREAL _2_CR     = 2.0f;
__device__ const CUSTOMREAL _3_CR     = 3.0f;
__device__ const CUSTOMREAL _4_CR     = 4.0f;

__device__ CUSTOMREAL my_square_cu(CUSTOMREAL const& x) {
    return x*x;
}

__device__ CUSTOMREAL calc_stencil_1st(CUSTOMREAL const& a, CUSTOMREAL const& b, CUSTOMREAL const& Dinv){
    return Dinv*(a-b);
}

__device__ CUSTOMREAL calc_stencil_3rd(CUSTOMREAL const& a, CUSTOMREAL const& b, CUSTOMREAL const& c, CUSTOMREAL const& d, CUSTOMREAL const& Dinv_half, CUSTOMREAL const& sign){
    CUSTOMREAL tmp1 = v_eps + my_square_cu(a-_2_CR*b+c);
    CUSTOMREAL tmp2 = v_eps + my_square_cu(d-_2_CR*a+b);
    CUSTOMREAL ww   = _1_CR/(_1_CR+_2_CR*my_square_cu(tmp1/tmp2));
    return sign*((_1_CR-ww)* (b-d)*Dinv_half + ww*(-_3_CR*a+_4_CR*b-c)*Dinv_half);
}

__device__ CUSTOMREAL cuda_calc_LF_Hamiltonian( \
                                            CUSTOMREAL const& fac_a_, \
                                            CUSTOMREAL const& fac_b_, \
                                            CUSTOMREAL const& fac_c_, \
                                            CUSTOMREAL const& fac_f_, \
                                            CUSTOMREAL const& T0r_, \
                                            CUSTOMREAL const& T0t_, \
                                            CUSTOMREAL const& T0p_, \
                                            CUSTOMREAL const& T0v_, \
                                            CUSTOMREAL& tau_, \
                                            CUSTOMREAL const& pp1, CUSTOMREAL& pp2, \
                                            CUSTOMREAL const& pt1, CUSTOMREAL& pt2, \
                                            CUSTOMREAL const& pr1, CUSTOMREAL& pr2 \
                                            ) {
    // LF Hamiltonian for T = T0 * tau
    return sqrt(
              fac_a_ * my_square_cu(T0r_ * tau_ + T0v_ * (pr1+pr2)/_2_CR) \
    +         fac_b_ * my_square_cu(T0t_ * tau_ + T0v_ * (pt1+pt2)/_2_CR) \
    +         fac_c_ * my_square_cu(T0p_ * tau_ + T0v_ * (pp1+pp2)/_2_CR) \
    -   _2_CR*fac_f_ * (T0t_ * tau_ + T0v_ * (pt1+pt2)/_2_CR) \
                     * (T0p_ * tau_ + T0v_ * (pp1+pp2)/_2_CR) \
    );
}

__global__ void cuda_do_sweep_level_kernel_1st(\
    const int i__j__k__[],\
    const int ip1j__k__[],\
    const int im1j__k__[],\
    const int i__jp1k__[],\
    const int i__jm1k__[],\
    const int i__j__kp1[],\
    const int i__j__km1[],\
    const CUSTOMREAL fac_a[], \
    const CUSTOMREAL fac_b[], \
    const CUSTOMREAL fac_c[], \
    const CUSTOMREAL fac_f[], \
    const CUSTOMREAL T0v[], \
    const CUSTOMREAL T0r[], \
    const CUSTOMREAL T0t[], \
    const CUSTOMREAL T0p[], \
    const CUSTOMREAL fun[], \
    const CUSTOMREAL changed[], \
    CUSTOMREAL tau[], \
    const int loc_I, \
    const int loc_J, \
    const int loc_K, \
    const CUSTOMREAL dr, \
    const CUSTOMREAL dt, \
    const CUSTOMREAL dp, \
    const int n_nodes_this_level, \
    const int i_start \
){

    unsigned int i_node = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
    //unsigned int i_node = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;

    if (i_node >= n_nodes_this_level) return;

    i_node += i_start;

    //if (i_node >= loc_I*loc_J*loc_K) return;

    if (changed[i_node] != _1_CR) return;

    CUSTOMREAL sigr = _1_CR*sqrt(fac_a[i_node])*T0v[i_node];
    CUSTOMREAL sigt = _1_CR*sqrt(fac_b[i_node])*T0v[i_node];
    CUSTOMREAL sigp = _1_CR*sqrt(fac_c[i_node])*T0v[i_node];
    CUSTOMREAL coe  = _1_CR/((sigr/dr)+(sigt/dt)+(sigp/dp));

    CUSTOMREAL pp1 = calc_stencil_1st(tau[i__j__k__[i_node]],tau[im1j__k__[i_node]], _1_CR/dp);
    CUSTOMREAL pp2 = calc_stencil_1st(tau[ip1j__k__[i_node]],tau[i__j__k__[i_node]], _1_CR/dp);

    CUSTOMREAL pt1 = calc_stencil_1st(tau[i__j__k__[i_node]],tau[i__jm1k__[i_node]], _1_CR/dt);
    CUSTOMREAL pt2 = calc_stencil_1st(tau[i__jp1k__[i_node]],tau[i__j__k__[i_node]], _1_CR/dt);

    CUSTOMREAL pr1 = calc_stencil_1st(tau[i__j__k__[i_node]],tau[i__j__km1[i_node]], _1_CR/dr);
    CUSTOMREAL pr2 = calc_stencil_1st(tau[i__j__kp1[i_node]],tau[i__j__k__[i_node]], _1_CR/dr);

    // LF Hamiltonian
    CUSTOMREAL Htau = cuda_calc_LF_Hamiltonian(\
                                               fac_a[i_node], \
                                               fac_b[i_node], \
                                               fac_c[i_node], \
                                               fac_f[i_node], \
                                               T0r[i_node], \
                                               T0t[i_node], \
                                               T0p[i_node], \
                                               T0v[i_node], \
                                               tau[i__j__k__[i_node]], \
                                               pp1, pp2, pt1, pt2, pr1, pr2);

    tau[i__j__k__[i_node]] += coe*((fun[i_node] - Htau) \
                                  +(sigr*(pr2-pr1) + sigt*(pt2-pt1) + sigp*(pp2-pp1))/_2_CR);

}

__global__ void cuda_do_sweep_level_kernel_3rd(\
    const int i__j__k__[],\
    const int ip1j__k__[],\
    const int im1j__k__[],\
    const int i__jp1k__[],\
    const int i__jm1k__[],\
    const int i__j__kp1[],\
    const int i__j__km1[],\
    const int ip2j__k__[],\
    const int im2j__k__[],\
    const int i__jp2k__[],\
    const int i__jm2k__[],\
    const int i__j__kp2[],\
    const int i__j__km2[],\
    const CUSTOMREAL fac_a[], \
    const CUSTOMREAL fac_b[], \
    const CUSTOMREAL fac_c[], \
    const CUSTOMREAL fac_f[], \
    const CUSTOMREAL T0v[], \
    const CUSTOMREAL T0r[], \
    const CUSTOMREAL T0t[], \
    const CUSTOMREAL T0p[], \
    const CUSTOMREAL fun[], \
    const CUSTOMREAL changed[], \
    CUSTOMREAL tau[], \
    const int loc_I, \
    const int loc_J, \
    const int loc_K, \
    const CUSTOMREAL dr, \
    const CUSTOMREAL dt, \
    const CUSTOMREAL dp,  \
    const int n_nodes_this_level, \
    const int i_start \
){

    CUSTOMREAL pp1, pp2, pt1, pt2, pr1, pr2;

    unsigned int i_node = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
    //unsigned int i_node = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;

    if (i_node >= n_nodes_this_level) return;

    i_node += i_start;
    //if (i_node >= loc_I*loc_J*loc_K) return;


    if (changed[i_node] != _1_CR) return;

    int k = i__j__k__[i_node]/(loc_I*loc_J);
    int j = (i__j__k__[i_node] - k*loc_I*loc_J)/loc_I;
    int i = i__j__k__[i_node] - k*loc_I*loc_J - j*loc_I;


    CUSTOMREAL DRinv = _1_CR/dr;
    CUSTOMREAL DTinv = _1_CR/dt;
    CUSTOMREAL DPinv = _1_CR/dp;
    CUSTOMREAL DRinv_half = DRinv*_0_5_CR;
    CUSTOMREAL DTinv_half = DTinv*_0_5_CR;
    CUSTOMREAL DPinv_half = DPinv*_0_5_CR;

    CUSTOMREAL sigr = _1_CR*sqrt(fac_a[i_node])*T0v[i_node];
    CUSTOMREAL sigt = _1_CR*sqrt(fac_b[i_node])*T0v[i_node];
    CUSTOMREAL sigp = _1_CR*sqrt(fac_c[i_node])*T0v[i_node];
    CUSTOMREAL coe  = _1_CR/((sigr/dr)+(sigt/dt)+(sigp/dp));

    // direction p
    if (i == 1) {
        pp1 = calc_stencil_1st(tau[i__j__k__[i_node]],tau[im1j__k__[i_node]],DPinv);
        pp2 = calc_stencil_3rd(tau[i__j__k__[i_node]],tau[ip1j__k__[i_node]],tau[ip2j__k__[i_node]],tau[im1j__k__[i_node]],DPinv_half, PLUS);
    } else if (i == loc_I-2) {
        pp1 = calc_stencil_3rd(tau[i__j__k__[i_node]],tau[im1j__k__[i_node]],tau[im2j__k__[i_node]],tau[ip1j__k__[i_node]],DPinv_half, MINUS);
        pp2 = calc_stencil_1st(tau[ip1j__k__[i_node]],tau[i__j__k__[i_node]],DPinv);
    } else {
        pp1 = calc_stencil_3rd(tau[i__j__k__[i_node]],tau[im1j__k__[i_node]],tau[im2j__k__[i_node]],tau[ip1j__k__[i_node]],DPinv_half, MINUS);
        pp2 = calc_stencil_3rd(tau[i__j__k__[i_node]],tau[ip1j__k__[i_node]],tau[ip2j__k__[i_node]],tau[im1j__k__[i_node]],DPinv_half, PLUS);
    }

    // direction t
    if (j == 1) {
        pt1 = calc_stencil_1st(tau[i__j__k__[i_node]],tau[i__jm1k__[i_node]],DTinv);
        pt2 = calc_stencil_3rd(tau[i__j__k__[i_node]],tau[i__jp1k__[i_node]],tau[i__jp2k__[i_node]],tau[i__jm1k__[i_node]],DTinv_half, PLUS);
    } else if (j == loc_J-2) {
        pt1 = calc_stencil_3rd(tau[i__j__k__[i_node]],tau[i__jm1k__[i_node]],tau[i__jm2k__[i_node]],tau[i__jp1k__[i_node]],DTinv_half, MINUS);
        pt2 = calc_stencil_1st(tau[i__jp1k__[i_node]],tau[i__j__k__[i_node]],DTinv);
    } else {
        pt1 = calc_stencil_3rd(tau[i__j__k__[i_node]],tau[i__jm1k__[i_node]],tau[i__jm2k__[i_node]],tau[i__jp1k__[i_node]],DTinv_half, MINUS);
        pt2 = calc_stencil_3rd(tau[i__j__k__[i_node]],tau[i__jp1k__[i_node]],tau[i__jp2k__[i_node]],tau[i__jm1k__[i_node]],DTinv_half, PLUS);
    }

    // direction r
    if (k == 1) {
        pr1 = calc_stencil_1st(tau[i__j__k__[i_node]],tau[i__j__km1[i_node]],DRinv);
        pr2 = calc_stencil_3rd(tau[i__j__k__[i_node]],tau[i__j__kp1[i_node]],tau[i__j__kp2[i_node]],tau[i__j__km1[i_node]],DRinv_half, PLUS);
    } else if (k == loc_K-2) {
        pr1 = calc_stencil_3rd(tau[i__j__k__[i_node]],tau[i__j__km1[i_node]],tau[i__j__km2[i_node]],tau[i__j__kp1[i_node]],DRinv_half, MINUS);
        pr2 = calc_stencil_1st(tau[i__j__kp1[i_node]],tau[i__j__k__[i_node]],DRinv);
    } else {
        pr1 = calc_stencil_3rd(tau[i__j__k__[i_node]],tau[i__j__km1[i_node]],tau[i__j__km2[i_node]],tau[i__j__kp1[i_node]],DRinv_half, MINUS);
        pr2 = calc_stencil_3rd(tau[i__j__k__[i_node]],tau[i__j__kp1[i_node]],tau[i__j__kp2[i_node]],tau[i__j__km1[i_node]],DRinv_half, PLUS);
    }

    CUSTOMREAL Htau = cuda_calc_LF_Hamiltonian(\
                                               fac_a[i_node], \
                                               fac_b[i_node], \
                                               fac_c[i_node], \
                                               fac_f[i_node], \
                                               T0r[i_node], \
                                               T0t[i_node], \
                                               T0p[i_node], \
                                               T0v[i_node], \
                                               tau[i__j__k__[i_node]], \
                                               pp1, pp2, pt1, pt2, pr1, pr2);

    tau[i__j__k__[i_node]] += coe*((fun[i_node] - Htau) \
                                  +(sigr*(pr2-pr1) + sigt*(pt2-pt1) + sigp*(pp2-pp1))/_2_CR);


}


void initialize_sweep_params(Grid_on_device* grid_dv){

    // check the numBlockPerSm and set the block size accordingly
    int numBlocksPerSm = 0;
    int block_size = CUDA_SWEEPING_BLOCK_SIZE;

    int device;
    cudaGetDevice(&device);

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, device);
    if(grid_dv->if_3rd_order)
        cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm, cuda_do_sweep_level_kernel_3rd, CUDA_SWEEPING_BLOCK_SIZE, 0);
    else
        cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm, cuda_do_sweep_level_kernel_1st, CUDA_SWEEPING_BLOCK_SIZE, 0);

    int max_cooperative_blocks = deviceProp.multiProcessorCount*numBlocksPerSm;

    grid_dv->threads_sweep_host = dim3(block_size, 1, 1);
    grid_dv->grid_sweep_host = dim3(max_cooperative_blocks, 1, 1);

    // spawn streams
    //grid_dv->level_streams = (cudaStream_t*)malloc(CUDA_MAX_NUM_STREAMS*sizeof(cudaStream_t));
    //for (int i = 0; i < CUDA_MAX_NUM_STREAMS; i++) {
    grid_dv->level_streams = (cudaStream_t*)malloc(grid_dv->n_levels_host*sizeof(cudaStream_t));
    for (int i = 0; i < grid_dv->n_levels_host; i++) {
        //cudaStreamCreate(&(grid_dv->level_streams[i]));
        // add null
        //cudaStreamCreateWithFlags(&(grid_dv->level_streams[i]), cudaStreamNonBlocking);
        grid_dv->level_streams[i] = nullptr;

    }


}


void finalize_sweep_params(Grid_on_device* grid_on_dv){
    // destroy streams
    //for (int i = 0; i < CUDA_MAX_NUM_STREAMS; i++) {
    //for (int i = 0; i < grid_on_dv->n_levels_host; i++) {
    //    cudaStreamDestroy(grid_on_dv->level_streams[i]);
    //}

    free(grid_on_dv->level_streams);
}


void run_kernel(Grid_on_device* grid_dv, int const& iswp, int& i_node_offset, int const& i_level, \
                dim3& grid_each, dim3& threads_each, int& n_nodes_this_level){

        int id_stream = i_level;// % CUDA_MAX_NUM_STREAMS;

        if (grid_dv->if_3rd_order) {
           if (iswp == 0){
                void *kernelArgs[]{\
                    &(grid_dv->vv_i__j__k___0), \
                    &(grid_dv->vv_ip1j__k___0), \
                    &(grid_dv->vv_im1j__k___0), \
                    &(grid_dv->vv_i__jp1k___0), \
                    &(grid_dv->vv_i__jm1k___0), \
                    &(grid_dv->vv_i__j__kp1_0), \
                    &(grid_dv->vv_i__j__km1_0), \
                    &(grid_dv->vv_ip2j__k___0), \
                    &(grid_dv->vv_im2j__k___0), \
                    &(grid_dv->vv_i__jp2k___0), \
                    &(grid_dv->vv_i__jm2k___0), \
                    &(grid_dv->vv_i__j__kp2_0), \
                    &(grid_dv->vv_i__j__km2_0), \
                    &(grid_dv->vv_fac_a_0    ), \
                    &(grid_dv->vv_fac_b_0    ), \
                    &(grid_dv->vv_fac_c_0    ), \
                    &(grid_dv->vv_fac_f_0    ), \
                    &(grid_dv->vv_T0v_0      ), \
                    &(grid_dv->vv_T0r_0      ), \
                    &(grid_dv->vv_T0t_0      ), \
                    &(grid_dv->vv_T0p_0      ), \
                    &(grid_dv->vv_fun_0      ), \
                    &(grid_dv->vv_change_0   ), \
                    &(grid_dv->tau), \
                    &(grid_dv->loc_I_host), \
                    &(grid_dv->loc_J_host), \
                    &(grid_dv->loc_K_host), \
                    &(grid_dv->dr_host), \
                    &(grid_dv->dt_host), \
                    &(grid_dv->dp_host), \
                    &n_nodes_this_level, \
                    &i_node_offset \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);

            } else if (iswp == 1){
                void* kernelArgs[]{\
                    &(grid_dv->vv_i__j__k___1), \
                    &(grid_dv->vv_i__jp1k___1), \
                    &(grid_dv->vv_i__jm1k___1), \
                    &(grid_dv->vv_i__j__kp1_1), \
                    &(grid_dv->vv_i__j__km1_1), \
                    &(grid_dv->vv_ip1j__k___1), \
                    &(grid_dv->vv_im1j__k___1), \
                    &(grid_dv->vv_ip2j__k___1), \
                    &(grid_dv->vv_im2j__k___1), \
                    &(grid_dv->vv_i__jp2k___1), \
                    &(grid_dv->vv_i__jm2k___1), \
                    &(grid_dv->vv_i__j__kp2_1), \
                    &(grid_dv->vv_i__j__km2_1), \
                    &(grid_dv->vv_fac_a_1    ), \
                    &(grid_dv->vv_fac_b_1    ), \
                    &(grid_dv->vv_fac_c_1    ), \
                    &(grid_dv->vv_fac_f_1    ), \
                    &(grid_dv->vv_T0v_1      ), \
                    &(grid_dv->vv_T0r_1      ), \
                    &(grid_dv->vv_T0t_1      ), \
                    &(grid_dv->vv_T0p_1      ), \
                    &(grid_dv->vv_fun_1      ), \
                    &(grid_dv->vv_change_1   ), \
                    &(grid_dv->tau), \
                    &(grid_dv->loc_I_host), \
                    &(grid_dv->loc_J_host), \
                    &(grid_dv->loc_K_host), \
                    &(grid_dv->dr_host), \
                    &(grid_dv->dt_host), \
                    &(grid_dv->dp_host), \
                    &n_nodes_this_level, \
                    &i_node_offset \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);

            } else if (iswp == 2){
                void* kernelArgs[]{\
                    &(grid_dv->vv_i__j__k___2), \
                    &(grid_dv->vv_i__j__kp1_2), \
                    &(grid_dv->vv_i__j__km1_2), \
                    &(grid_dv->vv_ip1j__k___2), \
                    &(grid_dv->vv_im1j__k___2), \
                    &(grid_dv->vv_i__jp1k___2), \
                    &(grid_dv->vv_i__jm1k___2), \
                    &(grid_dv->vv_ip2j__k___2), \
                    &(grid_dv->vv_im2j__k___2), \
                    &(grid_dv->vv_i__jp2k___2), \
                    &(grid_dv->vv_i__jm2k___2), \
                    &(grid_dv->vv_i__j__kp2_2), \
                    &(grid_dv->vv_i__j__km2_2), \
                    &(grid_dv->vv_fac_a_2), \
                    &(grid_dv->vv_fac_b_2), \
                    &(grid_dv->vv_fac_c_2), \
                    &(grid_dv->vv_fac_f_2), \
                    &(grid_dv->vv_T0v_2), \
                    &(grid_dv->vv_T0r_2), \
                    &(grid_dv->vv_T0t_2), \
                    &(grid_dv->vv_T0p_2), \
                    &(grid_dv->vv_fun_2), \
                    &(grid_dv->vv_change_2), \
                    &(grid_dv->tau), \
                    &(grid_dv->loc_I_host), \
                    &(grid_dv->loc_J_host), \
                    &(grid_dv->loc_K_host), \
                    &(grid_dv->dr_host), \
                    &(grid_dv->dt_host), \
                    &(grid_dv->dp_host), \
                    &n_nodes_this_level, \
                    &i_node_offset \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);

            } else if (iswp == 3){
                void* kernelArgs[]{\
                    &(grid_dv->vv_i__j__k___3), \
                    &(grid_dv->vv_ip1j__k___3), \
                    &(grid_dv->vv_im1j__k___3), \
                    &(grid_dv->vv_i__jp1k___3), \
                    &(grid_dv->vv_i__jm1k___3), \
                    &(grid_dv->vv_i__j__kp1_3), \
                    &(grid_dv->vv_i__j__km1_3), \
                    &(grid_dv->vv_ip2j__k___3), \
                    &(grid_dv->vv_im2j__k___3), \
                    &(grid_dv->vv_i__jp2k___3), \
                    &(grid_dv->vv_i__jm2k___3), \
                    &(grid_dv->vv_i__j__kp2_3), \
                    &(grid_dv->vv_i__j__km2_3), \
                    &(grid_dv->vv_fac_a_3), \
                    &(grid_dv->vv_fac_b_3), \
                    &(grid_dv->vv_fac_c_3), \
                    &(grid_dv->vv_fac_f_3), \
                    &(grid_dv->vv_T0v_3), \
                    &(grid_dv->vv_T0r_3), \
                    &(grid_dv->vv_T0t_3), \
                    &(grid_dv->vv_T0p_3), \
                    &(grid_dv->vv_fun_3), \
                    &(grid_dv->vv_change_3), \
                    &(grid_dv->tau), \
                    &(grid_dv->loc_I_host), \
                    &(grid_dv->loc_J_host), \
                    &(grid_dv->loc_K_host), \
                    &(grid_dv->dr_host), \
                    &(grid_dv->dt_host), \
                    &(grid_dv->dp_host), \
                    &n_nodes_this_level, \
                    &i_node_offset \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);

            } else if (iswp == 4){
                void* kernelArgs[]{\
                    &(grid_dv->vv_i__j__k___4), \
                    &(grid_dv->vv_ip1j__k___4), \
                    &(grid_dv->vv_im1j__k___4), \
                    &(grid_dv->vv_i__jp1k___4), \
                    &(grid_dv->vv_i__jm1k___4), \
                    &(grid_dv->vv_i__j__kp1_4), \
                    &(grid_dv->vv_i__j__km1_4), \
                    &(grid_dv->vv_ip2j__k___4), \
                    &(grid_dv->vv_im2j__k___4), \
                    &(grid_dv->vv_i__jp2k___4), \
                    &(grid_dv->vv_i__jm2k___4), \
                    &(grid_dv->vv_i__j__kp2_4), \
                    &(grid_dv->vv_i__j__km2_4), \
                    &(grid_dv->vv_fac_a_4), \
                    &(grid_dv->vv_fac_b_4), \
                    &(grid_dv->vv_fac_c_4), \
                    &(grid_dv->vv_fac_f_4), \
                    &(grid_dv->vv_T0v_4), \
                    &(grid_dv->vv_T0r_4), \
                    &(grid_dv->vv_T0t_4), \
                    &(grid_dv->vv_T0p_4), \
                    &(grid_dv->vv_fun_4), \
                    &(grid_dv->vv_change_4), \
                    &(grid_dv->tau), \
                    &(grid_dv->loc_I_host), \
                    &(grid_dv->loc_J_host), \
                    &(grid_dv->loc_K_host), \
                    &(grid_dv->dr_host), \
                    &(grid_dv->dt_host), \
                    &(grid_dv->dp_host), \
                    &n_nodes_this_level, \
                    &i_node_offset \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);

            } else if (iswp == 5) {
                void* kernelArgs[]{\
                    &(grid_dv->vv_i__j__k___5), \
                    &(grid_dv->vv_ip1j__k___5), \
                    &(grid_dv->vv_im1j__k___5), \
                    &(grid_dv->vv_i__jp1k___5), \
                    &(grid_dv->vv_i__jm1k___5), \
                    &(grid_dv->vv_i__j__kp1_5), \
                    &(grid_dv->vv_i__j__km1_5), \
                    &(grid_dv->vv_ip2j__k___5), \
                    &(grid_dv->vv_im2j__k___5), \
                    &(grid_dv->vv_i__jp2k___5), \
                    &(grid_dv->vv_i__jm2k___5), \
                    &(grid_dv->vv_i__j__kp2_5), \
                    &(grid_dv->vv_i__j__km2_5), \
                    &(grid_dv->vv_fac_a_5), \
                    &(grid_dv->vv_fac_b_5), \
                    &(grid_dv->vv_fac_c_5), \
                    &(grid_dv->vv_fac_f_5), \
                    &(grid_dv->vv_T0v_5), \
                    &(grid_dv->vv_T0r_5), \
                    &(grid_dv->vv_T0t_5), \
                    &(grid_dv->vv_T0p_5), \
                    &(grid_dv->vv_fun_5), \
                    &(grid_dv->vv_change_5), \
                    &(grid_dv->tau), \
                    &(grid_dv->loc_I_host), \
                    &(grid_dv->loc_J_host), \
                    &(grid_dv->loc_K_host), \
                    &(grid_dv->dr_host), \
                    &(grid_dv->dt_host), \
                    &(grid_dv->dp_host), \
                    &n_nodes_this_level, \
                    &i_node_offset \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);

            } else if (iswp == 6) {
                void* kernelArgs[]{\
                    &(grid_dv->vv_i__j__k___6), \
                    &(grid_dv->vv_ip1j__k___6), \
                    &(grid_dv->vv_im1j__k___6), \
                    &(grid_dv->vv_i__jp1k___6), \
                    &(grid_dv->vv_i__jm1k___6), \
                    &(grid_dv->vv_i__j__kp1_6), \
                    &(grid_dv->vv_i__j__km1_6), \
                    &(grid_dv->vv_ip2j__k___6), \
                    &(grid_dv->vv_im2j__k___6), \
                    &(grid_dv->vv_i__jp2k___6), \
                    &(grid_dv->vv_i__jm2k___6), \
                    &(grid_dv->vv_i__j__kp2_6), \
                    &(grid_dv->vv_i__j__km2_6), \
                    &(grid_dv->vv_fac_a_6), \
                    &(grid_dv->vv_fac_b_6), \
                    &(grid_dv->vv_fac_c_6), \
                    &(grid_dv->vv_fac_f_6), \
                    &(grid_dv->vv_T0v_6), \
                    &(grid_dv->vv_T0r_6), \
                    &(grid_dv->vv_T0t_6), \
                    &(grid_dv->vv_T0p_6), \
                    &(grid_dv->vv_fun_6), \
                    &(grid_dv->vv_change_6), \
                    &(grid_dv->tau), \
                    &(grid_dv->loc_I_host), \
                    &(grid_dv->loc_J_host), \
                    &(grid_dv->loc_K_host), \
                    &(grid_dv->dr_host), \
                    &(grid_dv->dt_host), \
                    &(grid_dv->dp_host), \
                    &n_nodes_this_level, \
                    &i_node_offset \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);

            } else {
                void* kernelArgs[]{\
                    &(grid_dv->vv_i__j__k___7), \
                    &(grid_dv->vv_ip1j__k___7), \
                    &(grid_dv->vv_im1j__k___7), \
                    &(grid_dv->vv_i__jp1k___7), \
                    &(grid_dv->vv_i__jm1k___7), \
                    &(grid_dv->vv_i__j__kp1_7), \
                    &(grid_dv->vv_i__j__km1_7), \
                    &(grid_dv->vv_ip2j__k___7), \
                    &(grid_dv->vv_im2j__k___7), \
                    &(grid_dv->vv_i__jp2k___7), \
                    &(grid_dv->vv_i__jm2k___7), \
                    &(grid_dv->vv_i__j__kp2_7), \
                    &(grid_dv->vv_i__j__km2_7), \
                    &(grid_dv->vv_fac_a_7), \
                    &(grid_dv->vv_fac_b_7), \
                    &(grid_dv->vv_fac_c_7), \
                    &(grid_dv->vv_fac_f_7), \
                    &(grid_dv->vv_T0v_7), \
                    &(grid_dv->vv_T0r_7), \
                    &(grid_dv->vv_T0t_7), \
                    &(grid_dv->vv_T0p_7), \
                    &(grid_dv->vv_fun_7), \
                    &(grid_dv->vv_change_7), \
                    &(grid_dv->tau), \
                    &(grid_dv->loc_I_host), \
                    &(grid_dv->loc_J_host), \
                    &(grid_dv->loc_K_host), \
                    &(grid_dv->dr_host), \
                    &(grid_dv->dt_host), \
                    &(grid_dv->dp_host), \
                    &n_nodes_this_level, \
                    &i_node_offset \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);

            }
        } else { // 1st order
            if (iswp == 0){
                void* kernelArgs[]{\
                    &(grid_dv->vv_i__j__k___0), \
                    &(grid_dv->vv_ip1j__k___0), \
                    &(grid_dv->vv_im1j__k___0), \
                    &(grid_dv->vv_i__jp1k___0), \
                    &(grid_dv->vv_i__jm1k___0), \
                    &(grid_dv->vv_i__j__kp1_0), \
                    &(grid_dv->vv_i__j__km1_0), \
                    &(grid_dv->vv_fac_a_0), \
                    &(grid_dv->vv_fac_b_0), \
                    &(grid_dv->vv_fac_c_0), \
                    &(grid_dv->vv_fac_f_0), \
                    &(grid_dv->vv_T0v_0), \
                    &(grid_dv->vv_T0r_0), \
                    &(grid_dv->vv_T0t_0), \
                    &(grid_dv->vv_T0p_0), \
                    &(grid_dv->vv_fun_0), \
                    &(grid_dv->vv_change_0), \
                    &(grid_dv->tau), \
                    &(grid_dv->loc_I_host), \
                    &(grid_dv->loc_J_host), \
                    &(grid_dv->loc_K_host), \
                    &(grid_dv->dr_host), \
                    &(grid_dv->dt_host), \
                    &(grid_dv->dp_host), \
                    &n_nodes_this_level, \
                    &i_node_offset \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_1st, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30000);

            } else if (iswp == 1){
                void* kernelArgs[]{\
                    &(grid_dv->vv_i__j__k___1), \
                    &(grid_dv->vv_i__jp1k___1), \
                    &(grid_dv->vv_i__jm1k___1), \
                    &(grid_dv->vv_i__j__kp1_1), \
                    &(grid_dv->vv_i__j__km1_1), \
                    &(grid_dv->vv_ip1j__k___1), \
                    &(grid_dv->vv_im1j__k___1), \
                    &(grid_dv->vv_fac_a_1), \
                    &(grid_dv->vv_fac_b_1), \
                    &(grid_dv->vv_fac_c_1), \
                    &(grid_dv->vv_fac_f_1), \
                    &(grid_dv->vv_T0v_1), \
                    &(grid_dv->vv_T0r_1), \
                    &(grid_dv->vv_T0t_1), \
                    &(grid_dv->vv_T0p_1), \
                    &(grid_dv->vv_fun_1), \
                    &(grid_dv->vv_change_1), \
                    &(grid_dv->tau), \
                    &(grid_dv->loc_I_host), \
                    &(grid_dv->loc_J_host), \
                    &(grid_dv->loc_K_host), \
                    &(grid_dv->dr_host), \
                    &(grid_dv->dt_host), \
                    &(grid_dv->dp_host), \
                    &n_nodes_this_level, \
                    &i_node_offset \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_1st, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);

            } else if (iswp == 2){
                void* kernelArgs[]{\
                    &(grid_dv->vv_i__j__k___2), \
                    &(grid_dv->vv_i__j__kp1_2), \
                    &(grid_dv->vv_i__j__km1_2), \
                    &(grid_dv->vv_ip1j__k___2), \
                    &(grid_dv->vv_im1j__k___2), \
                    &(grid_dv->vv_i__jp1k___2), \
                    &(grid_dv->vv_i__jm1k___2), \
                    &(grid_dv->vv_fac_a_2), \
                    &(grid_dv->vv_fac_b_2), \
                    &(grid_dv->vv_fac_c_2), \
                    &(grid_dv->vv_fac_f_2), \
                    &(grid_dv->vv_T0v_2), \
                    &(grid_dv->vv_T0r_2), \
                    &(grid_dv->vv_T0t_2), \
                    &(grid_dv->vv_T0p_2), \
                    &(grid_dv->vv_fun_2), \
                    &(grid_dv->vv_change_2), \
                    &(grid_dv->tau), \
                    &(grid_dv->loc_I_host), \
                    &(grid_dv->loc_J_host), \
                    &(grid_dv->loc_K_host), \
                    &(grid_dv->dr_host), \
                    &(grid_dv->dt_host), \
                    &(grid_dv->dp_host), \
                    &n_nodes_this_level, \
                    &i_node_offset \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_1st, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30002);

            } else if (iswp == 3){
                void* kernelArgs[]{\
                    &(grid_dv->vv_i__j__k___3), \
                    &(grid_dv->vv_ip1j__k___3), \
                    &(grid_dv->vv_im1j__k___3), \
                    &(grid_dv->vv_i__jp1k___3), \
                    &(grid_dv->vv_i__jm1k___3), \
                    &(grid_dv->vv_i__j__kp1_3), \
                    &(grid_dv->vv_i__j__km1_3), \
                    &(grid_dv->vv_fac_a_3), \
                    &(grid_dv->vv_fac_b_3), \
                    &(grid_dv->vv_fac_c_3), \
                    &(grid_dv->vv_fac_f_3), \
                    &(grid_dv->vv_T0v_3), \
                    &(grid_dv->vv_T0r_3), \
                    &(grid_dv->vv_T0t_3), \
                    &(grid_dv->vv_T0p_3), \
                    &(grid_dv->vv_fun_3), \
                    &(grid_dv->vv_change_3), \
                    &(grid_dv->tau), \
                    &(grid_dv->loc_I_host), \
                    &(grid_dv->loc_J_host), \
                    &(grid_dv->loc_K_host), \
                    &(grid_dv->dr_host), \
                    &(grid_dv->dt_host), \
                    &(grid_dv->dp_host), \
                    &n_nodes_this_level, \
                    &i_node_offset \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_1st, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30003);

            } else if (iswp == 4){
                void* kernelArgs[]{\
                    &(grid_dv->vv_i__j__k___4), \
                    &(grid_dv->vv_ip1j__k___4), \
                    &(grid_dv->vv_im1j__k___4), \
                    &(grid_dv->vv_i__jp1k___4), \
                    &(grid_dv->vv_i__jm1k___4), \
                    &(grid_dv->vv_i__j__kp1_4), \
                    &(grid_dv->vv_i__j__km1_4), \
                    &(grid_dv->vv_fac_a_4), \
                    &(grid_dv->vv_fac_b_4), \
                    &(grid_dv->vv_fac_c_4), \
                    &(grid_dv->vv_fac_f_4), \
                    &(grid_dv->vv_T0v_4), \
                    &(grid_dv->vv_T0r_4), \
                    &(grid_dv->vv_T0t_4), \
                    &(grid_dv->vv_T0p_4), \
                    &(grid_dv->vv_fun_4), \
                    &(grid_dv->vv_change_4), \
                    &(grid_dv->tau), \
                    &(grid_dv->loc_I_host), \
                    &(grid_dv->loc_J_host), \
                    &(grid_dv->loc_K_host), \
                    &(grid_dv->dr_host), \
                    &(grid_dv->dt_host), \
                    &(grid_dv->dp_host), \
                    &n_nodes_this_level, \
                    &i_node_offset \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_1st, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30004);

            } else if (iswp == 5) {
                void* kernelArgs[]{\
                    &(grid_dv->vv_i__j__k___5), \
                    &(grid_dv->vv_ip1j__k___5), \
                    &(grid_dv->vv_im1j__k___5), \
                    &(grid_dv->vv_i__jp1k___5), \
                    &(grid_dv->vv_i__jm1k___5), \
                    &(grid_dv->vv_i__j__kp1_5), \
                    &(grid_dv->vv_i__j__km1_5), \
                    &(grid_dv->vv_fac_a_5), \
                    &(grid_dv->vv_fac_b_5), \
                    &(grid_dv->vv_fac_c_5), \
                    &(grid_dv->vv_fac_f_5), \
                    &(grid_dv->vv_T0v_5), \
                    &(grid_dv->vv_T0r_5), \
                    &(grid_dv->vv_T0t_5), \
                    &(grid_dv->vv_T0p_5), \
                    &(grid_dv->vv_fun_5), \
                    &(grid_dv->vv_change_5), \
                    &(grid_dv->tau), \
                    &(grid_dv->loc_I_host), \
                    &(grid_dv->loc_J_host), \
                    &(grid_dv->loc_K_host), \
                    &(grid_dv->dr_host), \
                    &(grid_dv->dt_host), \
                    &(grid_dv->dp_host), \
                    &n_nodes_this_level, \
                    &i_node_offset \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_1st, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30005);

            } else if (iswp == 6) {
                void* kernelArgs[]{\
                    &(grid_dv->vv_i__j__k___6), \
                    &(grid_dv->vv_ip1j__k___6), \
                    &(grid_dv->vv_im1j__k___6), \
                    &(grid_dv->vv_i__jp1k___6), \
                    &(grid_dv->vv_i__jm1k___6), \
                    &(grid_dv->vv_i__j__kp1_6), \
                    &(grid_dv->vv_i__j__km1_6), \
                    &(grid_dv->vv_fac_a_6), \
                    &(grid_dv->vv_fac_b_6), \
                    &(grid_dv->vv_fac_c_6), \
                    &(grid_dv->vv_fac_f_6), \
                    &(grid_dv->vv_T0v_6), \
                    &(grid_dv->vv_T0r_6), \
                    &(grid_dv->vv_T0t_6), \
                    &(grid_dv->vv_T0p_6), \
                    &(grid_dv->vv_fun_6), \
                    &(grid_dv->vv_change_6), \
                    &(grid_dv->tau), \
                    &(grid_dv->loc_I_host), \
                    &(grid_dv->loc_J_host), \
                    &(grid_dv->loc_K_host), \
                    &(grid_dv->dr_host), \
                    &(grid_dv->dt_host), \
                    &(grid_dv->dp_host), \
                    &n_nodes_this_level, \
                    &i_node_offset \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_1st, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30006);


            } else {
                void* kernelArgs[]{\
                    &(grid_dv->vv_i__j__k___7), \
                    &(grid_dv->vv_ip1j__k___7), \
                    &(grid_dv->vv_im1j__k___7), \
                    &(grid_dv->vv_i__jp1k___7), \
                    &(grid_dv->vv_i__jm1k___7), \
                    &(grid_dv->vv_i__j__kp1_7), \
                    &(grid_dv->vv_i__j__km1_7), \
                    &(grid_dv->vv_fac_a_7    ), \
                    &(grid_dv->vv_fac_b_7    ), \
                    &(grid_dv->vv_fac_c_7    ), \
                    &(grid_dv->vv_fac_f_7    ), \
                    &(grid_dv->vv_T0v_7      ), \
                    &(grid_dv->vv_T0r_7      ), \
                    &(grid_dv->vv_T0t_7      ), \
                    &(grid_dv->vv_T0p_7      ), \
                    &(grid_dv->vv_fun_7      ), \
                    &(grid_dv->vv_change_7   ), \
                    &(grid_dv->tau), \
                    &(grid_dv->loc_I_host), \
                    &(grid_dv->loc_J_host), \
                    &(grid_dv->loc_K_host), \
                    &(grid_dv->dr_host), \
                    &(grid_dv->dt_host), \
                    &(grid_dv->dp_host), \
                    &n_nodes_this_level, \
                    &i_node_offset \
                };

                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_1st, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30007);

            }
        }

        // synchronize all streams
        //print_CUDA_error_if_any(cudaStreamSynchronize(grid_dv->level_streams[id_stream]), 30008);
}


// this function calculate all levels of one single sweep direction
void cuda_run_iteration_forward(Grid_on_device* grid_dv, int const& iswp){

    initialize_sweep_params(grid_dv);

    int block_size = CUDA_SWEEPING_BLOCK_SIZE;
    int num_blocks_x, num_blocks_y;
    int actual_end_level = grid_dv->n_levels_host;
    int i_node_offset=0;

    for (size_t i_level = 0; i_level < actual_end_level; i_level++){
        get_block_xy(ceil(grid_dv->n_nodes_on_levels_host[i_level]/block_size+0.5), &num_blocks_x, &num_blocks_y);
        dim3 grid_each(num_blocks_x, num_blocks_y);
        dim3 threads_each(block_size, 1, 1);

        run_kernel(grid_dv, iswp, i_node_offset, i_level, grid_each, threads_each, grid_dv->n_nodes_on_levels_host[i_level]);
        //run_kernel(grid_dv, iswp, i_node_offset, i_level, grid_dv->grid_sweep_host, grid_dv->threads_sweep_host, grid_dv->n_nodes_on_levels_host[i_level]);

        i_node_offset += grid_dv->n_nodes_on_levels_host[i_level];
    }

    finalize_sweep_params(grid_dv);

    // check memory leak
    //print_memory_usage();

}