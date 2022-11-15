#include "iterator_wrapper.cuh"

#define PLUS 1.0
#define MINUS -1.0
#define v_eps 1e-12


__device__ CUSTOMREAL my_square_cu(CUSTOMREAL const& x) {
    return x*x;
}

__device__ CUSTOMREAL calc_stencil_1st(CUSTOMREAL const& a, CUSTOMREAL const& b, CUSTOMREAL const& Dinv){
    return Dinv*(a-b);
}

__device__ CUSTOMREAL calc_stencil_3rd(CUSTOMREAL const& a, CUSTOMREAL const& b, CUSTOMREAL const& c, CUSTOMREAL const& d, CUSTOMREAL const& Dinv_half, int const& sign){
    CUSTOMREAL tmp1 = v_eps + my_square_cu(a-2.0*b+c);
    CUSTOMREAL tmp2 = v_eps + my_square_cu(d-2.0*a+b);
    CUSTOMREAL ww   = 1.0/(1.0+2.0*my_square_cu(tmp1/tmp2));
    return sign*((1.0-ww)* (b-d)*Dinv_half + ww*(-3.0*a+4.0*b-c)*Dinv_half);
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
              fac_a_ * (T0r_ * tau_ + T0v_ * (pr1+pr2)/2.0) * (T0r_ * tau_ + T0v_ * (pr1+pr2)/2.0)\
    +         fac_b_ * (T0t_ * tau_ + T0v_ * (pt1+pt2)/2.0) * (T0t_ * tau_ + T0v_ * (pt1+pt2)/2.0)\
    +         fac_c_ * (T0p_ * tau_ + T0v_ * (pp1+pp2)/2.0) * (T0p_ * tau_ + T0v_ * (pp1+pp2)/2.0)\
    -     2.0*fac_f_ * (T0t_ * tau_ + T0v_ * (pt1+pt2)/2.0) \
                     * (T0p_ * tau_ + T0v_ * (pp1+pp2)/2.0) \
    );
}

__global__ void cuda_do_sweep_level_kernel_1st(\
    int* i__j__k__,\
    int* ip1j__k__,\
    int* im1j__k__,\
    int* i__jp1k__,\
    int* i__jm1k__,\
    int* i__j__kp1,\
    int* i__j__km1,\
    CUSTOMREAL* fac_a, \
    CUSTOMREAL* fac_b, \
    CUSTOMREAL* fac_c, \
    CUSTOMREAL* fac_f, \
    CUSTOMREAL* T0v, \
    CUSTOMREAL* T0r, \
    CUSTOMREAL* T0t, \
    CUSTOMREAL* T0p, \
    CUSTOMREAL* fun, \
    CUSTOMREAL* changed, \
    CUSTOMREAL* tau, \
    int loc_I, \
    int loc_J, \
    int loc_K, \
    CUSTOMREAL dr, \
    CUSTOMREAL dt, \
    CUSTOMREAL dp, \
    int n_nodes_this_level \
){
    unsigned int i_node = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;

    if (i_node > n_nodes_this_level) return;

    CUSTOMREAL sigr = 1.0*sqrt(fac_a[i_node])*T0v[i_node];
    CUSTOMREAL sigt = 1.0*sqrt(fac_b[i_node])*T0v[i_node];
    CUSTOMREAL sigp = 1.0*sqrt(fac_c[i_node])*T0v[i_node];
    CUSTOMREAL coe  = 1.0/((sigr/dr)+(sigt/dt)+(sigp/dp));

    CUSTOMREAL pp1 = calc_stencil_1st(tau[i__j__k__[i_node]],tau[im1j__k__[i_node]],1.0/dp);
    CUSTOMREAL pp2 = calc_stencil_1st(tau[ip1j__k__[i_node]],tau[i__j__k__[i_node]],1.0/dp);

    CUSTOMREAL pt1 = calc_stencil_1st(tau[i__j__k__[i_node]],tau[i__jm1k__[i_node]],1.0/dt);
    CUSTOMREAL pt2 = calc_stencil_1st(tau[i__jp1k__[i_node]],tau[i__j__k__[i_node]],1.0/dt);

    CUSTOMREAL pr1 = calc_stencil_1st(tau[i__j__k__[i_node]],tau[i__j__km1[i_node]],1.0/dr);
    CUSTOMREAL pr2 = calc_stencil_1st(tau[i__j__kp1[i_node]],tau[i__j__k__[i_node]],1.0/dr);

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

    tau[i__j__k__[i_node]] += coe*(fun[i_node] - Htau) \
                            + coe*(sigr*(pr2-pr1)/2.0 + sigt*(pt2-pt1)/2.0 + sigp*(pp2-pp1)/2.0);

}

__global__ void cuda_do_sweep_level_kernel_3rd(\
    int* i__j__k__,\
    int* ip1j__k__,\
    int* im1j__k__,\
    int* i__jp1k__,\
    int* i__jm1k__,\
    int* i__j__kp1,\
    int* i__j__km1,\
    int* ip2j__k__,\
    int* im2j__k__,\
    int* i__jp2k__,\
    int* i__jm2k__,\
    int* i__j__kp2,\
    int* i__j__km2,\
    CUSTOMREAL* fac_a, \
    CUSTOMREAL* fac_b, \
    CUSTOMREAL* fac_c, \
    CUSTOMREAL* fac_f, \
    CUSTOMREAL* T0v, \
    CUSTOMREAL* T0r, \
    CUSTOMREAL* T0t, \
    CUSTOMREAL* T0p, \
    CUSTOMREAL* fun, \
    CUSTOMREAL* changed, \
    CUSTOMREAL* tau, \
    int loc_I, \
    int loc_J, \
    int loc_K, \
    CUSTOMREAL dr, \
    CUSTOMREAL dt, \
    CUSTOMREAL dp,  \
    int n_nodes_this_level \
){

    CUSTOMREAL pp1, pp2, pt1, pt2, pr1, pr2;

    unsigned int i_node = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;

    if (i_node > n_nodes_this_level) return;

    int k = i_node/(loc_I*loc_J);
    int j = (i_node - k*loc_I*loc_J)/loc_I;
    int i = i_node - k*loc_I*loc_J - j*loc_I;
    CUSTOMREAL DRinv = 1.0/dr;
    CUSTOMREAL DTinv = 1.0/dt;
    CUSTOMREAL DPinv = 1.0/dp;
    CUSTOMREAL DRinv_half = DRinv*0.5;
    CUSTOMREAL DTinv_half = DTinv*0.5;
    CUSTOMREAL DPinv_half = DPinv*0.5;

    CUSTOMREAL sigr = 1.0*sqrt(fac_a[i__j__k__[i_node]])*T0v[i__j__k__[i_node]];
    CUSTOMREAL sigt = 1.0*sqrt(fac_b[i__j__k__[i_node]])*T0v[i__j__k__[i_node]];
    CUSTOMREAL sigp = 1.0*sqrt(fac_c[i__j__k__[i_node]])*T0v[i__j__k__[i_node]];
    CUSTOMREAL coe  = 1.0/((sigr/dr)+(sigt/dt)+(sigp/dp));

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
                                               fac_a[i__j__k__[i_node]], \
                                               fac_b[i__j__k__[i_node]], \
                                               fac_c[i__j__k__[i_node]], \
                                               fac_f[i__j__k__[i_node]], \
                                               T0r[i__j__k__[i_node]], \
                                               T0t[i__j__k__[i_node]], \
                                               T0p[i__j__k__[i_node]], \
                                               T0v[i__j__k__[i_node]], \
                                               tau[i__j__k__[i_node]], \
                                               pp1, pp2, pt1, pt2, pr1, pr2);

    tau[i__j__k__[i_node]] += coe*(fun[i__j__k__[i_node]] - Htau) \
                            + coe*(sigr*(pr2-pr1)/2.0 + sigt*(pt2-pt1)/2.0 + sigp*(pp2-pp1)/2.0);


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
    grid_dv->level_streams = (cudaStream_t*)malloc(CUDA_MAX_NUM_STREAMS*sizeof(cudaStream_t));
    for (int i = 0; i < CUDA_MAX_NUM_STREAMS; i++) {
        cudaStreamCreate(&(grid_dv->level_streams[i]));
    }


}

//        if (grid_dv->if_3rd_order) {
//            print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);
//        } else {
//            print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_1st, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30002);
//        }


void run_kernel(Grid_on_device* grid_dv, int const& iswp, int const& i_node_offset, int const& i_level, \
                dim3& grid_each, dim3& threads_each, int& n_nodes_this_level){

        int id_stream = i_level % CUDA_MAX_NUM_STREAMS;

        if (grid_dv->if_3rd_order) {
           if (iswp == 0){
                void *kernelArgs[]= {\
                    &(grid_dv->vv_i__j__k___0[i_node_offset]), \
                    &(grid_dv->vv_ip1j__k___0[i_node_offset]), \
                    &(grid_dv->vv_im1j__k___0[i_node_offset]), \
                    &(grid_dv->vv_i__jp1k___0[i_node_offset]), \
                    &(grid_dv->vv_i__jm1k___0[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp1_0[i_node_offset]), \
                    &(grid_dv->vv_i__j__km1_0[i_node_offset]), \
                    &(grid_dv->vv_ip2j__k___0[i_node_offset]), \
                    &(grid_dv->vv_im2j__k___0[i_node_offset]), \
                    &(grid_dv->vv_i__jp2k___0[i_node_offset]), \
                    &(grid_dv->vv_i__jm2k___0[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp2_0[i_node_offset]), \
                    &(grid_dv->vv_i__j__km2_0[i_node_offset]), \
                    &(grid_dv->vv_fac_a_0[i_node_offset]), \
                    &(grid_dv->vv_fac_b_0[i_node_offset]), \
                    &(grid_dv->vv_fac_c_0[i_node_offset]), \
                    &(grid_dv->vv_fac_f_0[i_node_offset]), \
                    &(grid_dv->vv_T0v_0[i_node_offset]), \
                    &(grid_dv->vv_T0r_0[i_node_offset]), \
                    &(grid_dv->vv_T0t_0[i_node_offset]), \
                    &(grid_dv->vv_T0p_0[i_node_offset]), \
                    &(grid_dv->vv_fun_0[i_node_offset]), \
                    &(grid_dv->vv_change_0[i_node_offset]), \
                    &(grid_dv->tau), \
                    &grid_dv->loc_I_host, \
                    &grid_dv->loc_J_host, \
                    &grid_dv->loc_K_host, \
                    &grid_dv->dr_host, \
                    &grid_dv->dt_host, \
                    &grid_dv->dp_host, \
                    &n_nodes_this_level \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);

            } else if (iswp == 1){
                void* kernelArgs[]= {\
                    &(grid_dv->vv_i__j__k___1[i_node_offset]), \
                    &(grid_dv->vv_i__jp1k___1[i_node_offset]), \
                    &(grid_dv->vv_i__jm1k___1[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp1_1[i_node_offset]), \
                    &(grid_dv->vv_i__j__km1_1[i_node_offset]), \
                    &(grid_dv->vv_ip1j__k___1[i_node_offset]), \
                    &(grid_dv->vv_im1j__k___1[i_node_offset]), \
                    &(grid_dv->vv_ip2j__k___1[i_node_offset]), \
                    &(grid_dv->vv_im2j__k___1[i_node_offset]), \
                    &(grid_dv->vv_i__jp2k___1[i_node_offset]), \
                    &(grid_dv->vv_i__jm2k___1[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp2_1[i_node_offset]), \
                    &(grid_dv->vv_i__j__km2_1[i_node_offset]), \
                    &(grid_dv->vv_fac_a_1[i_node_offset]), \
                    &(grid_dv->vv_fac_b_1[i_node_offset]), \
                    &(grid_dv->vv_fac_c_1[i_node_offset]), \
                    &(grid_dv->vv_fac_f_1[i_node_offset]), \
                    &(grid_dv->vv_T0v_1[i_node_offset]), \
                    &(grid_dv->vv_T0r_1[i_node_offset]), \
                    &(grid_dv->vv_T0t_1[i_node_offset]), \
                    &(grid_dv->vv_T0p_1[i_node_offset]), \
                    &(grid_dv->vv_fun_1[i_node_offset]), \
                    &(grid_dv->vv_change_1[i_node_offset]), \
                    &(grid_dv->tau), \
                    &grid_dv->loc_I_host, \
                    &grid_dv->loc_J_host, \
                    &grid_dv->loc_K_host, \
                    &grid_dv->dr_host, \
                    &grid_dv->dt_host, \
                    &grid_dv->dp_host, \
                    &n_nodes_this_level \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);

            } else if (iswp == 2){
                void* kernelArgs[]= {\
                    &(grid_dv->vv_i__j__k___2[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp1_2[i_node_offset]), \
                    &(grid_dv->vv_i__j__km1_2[i_node_offset]), \
                    &(grid_dv->vv_ip1j__k___2[i_node_offset]), \
                    &(grid_dv->vv_im1j__k___2[i_node_offset]), \
                    &(grid_dv->vv_i__jp1k___2[i_node_offset]), \
                    &(grid_dv->vv_i__jm1k___2[i_node_offset]), \
                    &(grid_dv->vv_ip2j__k___2[i_node_offset]), \
                    &(grid_dv->vv_im2j__k___2[i_node_offset]), \
                    &(grid_dv->vv_i__jp2k___2[i_node_offset]), \
                    &(grid_dv->vv_i__jm2k___2[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp2_2[i_node_offset]), \
                    &(grid_dv->vv_i__j__km2_2[i_node_offset]), \
                    &(grid_dv->vv_fac_a_2[i_node_offset]), \
                    &(grid_dv->vv_fac_b_2[i_node_offset]), \
                    &(grid_dv->vv_fac_c_2[i_node_offset]), \
                    &(grid_dv->vv_fac_f_2[i_node_offset]), \
                    &(grid_dv->vv_T0v_2[i_node_offset]), \
                    &(grid_dv->vv_T0r_2[i_node_offset]), \
                    &(grid_dv->vv_T0t_2[i_node_offset]), \
                    &(grid_dv->vv_T0p_2[i_node_offset]), \
                    &(grid_dv->vv_fun_2[i_node_offset]), \
                    &(grid_dv->vv_change_2[i_node_offset]), \
                    &(grid_dv->tau), \
                    &grid_dv->loc_I_host, \
                    &grid_dv->loc_J_host, \
                    &grid_dv->loc_K_host, \
                    &grid_dv->dr_host, \
                    &grid_dv->dt_host, \
                    &grid_dv->dp_host, \
                    &n_nodes_this_level \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);

            } else if (iswp == 3){
                void* kernelArgs[]= {\
                    &(grid_dv->vv_i__j__k___3[i_node_offset]), \
                    &(grid_dv->vv_ip1j__k___3[i_node_offset]), \
                    &(grid_dv->vv_im1j__k___3[i_node_offset]), \
                    &(grid_dv->vv_i__jp1k___3[i_node_offset]), \
                    &(grid_dv->vv_i__jm1k___3[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp1_3[i_node_offset]), \
                    &(grid_dv->vv_i__j__km1_3[i_node_offset]), \
                    &(grid_dv->vv_ip2j__k___3[i_node_offset]), \
                    &(grid_dv->vv_im2j__k___3[i_node_offset]), \
                    &(grid_dv->vv_i__jp2k___3[i_node_offset]), \
                    &(grid_dv->vv_i__jm2k___3[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp2_3[i_node_offset]), \
                    &(grid_dv->vv_i__j__km2_3[i_node_offset]), \
                    &(grid_dv->vv_fac_a_3[i_node_offset]), \
                    &(grid_dv->vv_fac_b_3[i_node_offset]), \
                    &(grid_dv->vv_fac_c_3[i_node_offset]), \
                    &(grid_dv->vv_fac_f_3[i_node_offset]), \
                    &(grid_dv->vv_T0v_3[i_node_offset]), \
                    &(grid_dv->vv_T0r_3[i_node_offset]), \
                    &(grid_dv->vv_T0t_3[i_node_offset]), \
                    &(grid_dv->vv_T0p_3[i_node_offset]), \
                    &(grid_dv->vv_fun_3[i_node_offset]), \
                    &(grid_dv->vv_change_3[i_node_offset]), \
                    &(grid_dv->tau), \
                    &grid_dv->loc_I_host, \
                    &grid_dv->loc_J_host, \
                    &grid_dv->loc_K_host, \
                    &grid_dv->dr_host, \
                    &grid_dv->dt_host, \
                    &grid_dv->dp_host, \
                    &n_nodes_this_level \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);

            } else if (iswp == 4){
                void* kernelArgs[]= {\
                    &(grid_dv->vv_i__j__k___4[i_node_offset]), \
                    &(grid_dv->vv_ip1j__k___4[i_node_offset]), \
                    &(grid_dv->vv_im1j__k___4[i_node_offset]), \
                    &(grid_dv->vv_i__jp1k___4[i_node_offset]), \
                    &(grid_dv->vv_i__jm1k___4[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp1_4[i_node_offset]), \
                    &(grid_dv->vv_i__j__km1_4[i_node_offset]), \
                    &(grid_dv->vv_ip2j__k___4[i_node_offset]), \
                    &(grid_dv->vv_im2j__k___4[i_node_offset]), \
                    &(grid_dv->vv_i__jp2k___4[i_node_offset]), \
                    &(grid_dv->vv_i__jm2k___4[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp2_4[i_node_offset]), \
                    &(grid_dv->vv_i__j__km2_4[i_node_offset]), \
                    &(grid_dv->vv_fac_a_4[i_node_offset]), \
                    &(grid_dv->vv_fac_b_4[i_node_offset]), \
                    &(grid_dv->vv_fac_c_4[i_node_offset]), \
                    &(grid_dv->vv_fac_f_4[i_node_offset]), \
                    &(grid_dv->vv_T0v_4[i_node_offset]), \
                    &(grid_dv->vv_T0r_4[i_node_offset]), \
                    &(grid_dv->vv_T0t_4[i_node_offset]), \
                    &(grid_dv->vv_T0p_4[i_node_offset]), \
                    &(grid_dv->vv_fun_4[i_node_offset]), \
                    &(grid_dv->vv_change_4[i_node_offset]), \
                    &(grid_dv->tau), \
                    &grid_dv->loc_I_host, \
                    &grid_dv->loc_J_host, \
                    &grid_dv->loc_K_host, \
                    &grid_dv->dr_host, \
                    &grid_dv->dt_host, \
                    &grid_dv->dp_host, \
                    &n_nodes_this_level \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);

            } else if (iswp == 5) {
                void* kernelArgs[]= {\
                    &(grid_dv->vv_i__j__k___5[i_node_offset]), \
                    &(grid_dv->vv_ip1j__k___5[i_node_offset]), \
                    &(grid_dv->vv_im1j__k___5[i_node_offset]), \
                    &(grid_dv->vv_i__jp1k___5[i_node_offset]), \
                    &(grid_dv->vv_i__jm1k___5[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp1_5[i_node_offset]), \
                    &(grid_dv->vv_i__j__km1_5[i_node_offset]), \
                    &(grid_dv->vv_ip2j__k___5[i_node_offset]), \
                    &(grid_dv->vv_im2j__k___5[i_node_offset]), \
                    &(grid_dv->vv_i__jp2k___5[i_node_offset]), \
                    &(grid_dv->vv_i__jm2k___5[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp2_5[i_node_offset]), \
                    &(grid_dv->vv_i__j__km2_5[i_node_offset]), \
                    &(grid_dv->vv_fac_a_5[i_node_offset]), \
                    &(grid_dv->vv_fac_b_5[i_node_offset]), \
                    &(grid_dv->vv_fac_c_5[i_node_offset]), \
                    &(grid_dv->vv_fac_f_5[i_node_offset]), \
                    &(grid_dv->vv_T0v_5[i_node_offset]), \
                    &(grid_dv->vv_T0r_5[i_node_offset]), \
                    &(grid_dv->vv_T0t_5[i_node_offset]), \
                    &(grid_dv->vv_T0p_5[i_node_offset]), \
                    &(grid_dv->vv_fun_5[i_node_offset]), \
                    &(grid_dv->vv_change_5[i_node_offset]), \
                    &(grid_dv->tau), \
                    &grid_dv->loc_I_host, \
                    &grid_dv->loc_J_host, \
                    &grid_dv->loc_K_host, \
                    &grid_dv->dr_host, \
                    &grid_dv->dt_host, \
                    &grid_dv->dp_host, \
                    &n_nodes_this_level \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);

            } else if (iswp == 6) {
                void* kernelArgs[]= {\
                    &(grid_dv->vv_i__j__k___6[i_node_offset]), \
                    &(grid_dv->vv_ip1j__k___6[i_node_offset]), \
                    &(grid_dv->vv_im1j__k___6[i_node_offset]), \
                    &(grid_dv->vv_i__jp1k___6[i_node_offset]), \
                    &(grid_dv->vv_i__jm1k___6[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp1_6[i_node_offset]), \
                    &(grid_dv->vv_i__j__km1_6[i_node_offset]), \
                    &(grid_dv->vv_ip2j__k___6[i_node_offset]), \
                    &(grid_dv->vv_im2j__k___6[i_node_offset]), \
                    &(grid_dv->vv_i__jp2k___6[i_node_offset]), \
                    &(grid_dv->vv_i__jm2k___6[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp2_6[i_node_offset]), \
                    &(grid_dv->vv_i__j__km2_6[i_node_offset]), \
                    &(grid_dv->vv_fac_a_6[i_node_offset]), \
                    &(grid_dv->vv_fac_b_6[i_node_offset]), \
                    &(grid_dv->vv_fac_c_6[i_node_offset]), \
                    &(grid_dv->vv_fac_f_6[i_node_offset]), \
                    &(grid_dv->vv_T0v_6[i_node_offset]), \
                    &(grid_dv->vv_T0r_6[i_node_offset]), \
                    &(grid_dv->vv_T0t_6[i_node_offset]), \
                    &(grid_dv->vv_T0p_6[i_node_offset]), \
                    &(grid_dv->vv_fun_6[i_node_offset]), \
                    &(grid_dv->vv_change_6[i_node_offset]), \
                    &(grid_dv->tau), \
                    &grid_dv->loc_I_host, \
                    &grid_dv->loc_J_host, \
                    &grid_dv->loc_K_host, \
                    &grid_dv->dr_host, \
                    &grid_dv->dt_host, \
                    &grid_dv->dp_host, \
                    &n_nodes_this_level \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);

            } else {
                void* kernelArgs[]= {\
                    &(grid_dv->vv_i__j__k___7[i_node_offset]), \
                    &(grid_dv->vv_ip1j__k___7[i_node_offset]), \
                    &(grid_dv->vv_im1j__k___7[i_node_offset]), \
                    &(grid_dv->vv_i__jp1k___7[i_node_offset]), \
                    &(grid_dv->vv_i__jm1k___7[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp1_7[i_node_offset]), \
                    &(grid_dv->vv_i__j__km1_7[i_node_offset]), \
                    &(grid_dv->vv_ip2j__k___7[i_node_offset]), \
                    &(grid_dv->vv_im2j__k___7[i_node_offset]), \
                    &(grid_dv->vv_i__jp2k___7[i_node_offset]), \
                    &(grid_dv->vv_i__jm2k___7[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp2_7[i_node_offset]), \
                    &(grid_dv->vv_i__j__km2_7[i_node_offset]), \
                    &(grid_dv->vv_fac_a_7[i_node_offset]), \
                    &(grid_dv->vv_fac_b_7[i_node_offset]), \
                    &(grid_dv->vv_fac_c_7[i_node_offset]), \
                    &(grid_dv->vv_fac_f_7[i_node_offset]), \
                    &(grid_dv->vv_T0v_7[i_node_offset]), \
                    &(grid_dv->vv_T0r_7[i_node_offset]), \
                    &(grid_dv->vv_T0t_7[i_node_offset]), \
                    &(grid_dv->vv_T0p_7[i_node_offset]), \
                    &(grid_dv->vv_fun_7[i_node_offset]), \
                    &(grid_dv->vv_change_7[i_node_offset]), \
                    &(grid_dv->tau), \
                    &grid_dv->loc_I_host, \
                    &grid_dv->loc_J_host, \
                    &grid_dv->loc_K_host, \
                    &grid_dv->dr_host, \
                    &grid_dv->dt_host, \
                    &grid_dv->dp_host, \
                    &n_nodes_this_level \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);

            }
        } else { // 1st order
            if (iswp == 0){
                void* kernelArgs[]= {\
                    &(grid_dv->vv_i__j__k___0[i_node_offset]), \
                    &(grid_dv->vv_ip1j__k___0[i_node_offset]), \
                    &(grid_dv->vv_im1j__k___0[i_node_offset]), \
                    &(grid_dv->vv_i__jp1k___0[i_node_offset]), \
                    &(grid_dv->vv_i__jm1k___0[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp1_0[i_node_offset]), \
                    &(grid_dv->vv_i__j__km1_0[i_node_offset]), \
                    &(grid_dv->vv_fac_a_0[i_node_offset]), \
                    &(grid_dv->vv_fac_b_0[i_node_offset]), \
                    &(grid_dv->vv_fac_c_0[i_node_offset]), \
                    &(grid_dv->vv_fac_f_0[i_node_offset]), \
                    &(grid_dv->vv_T0v_0[i_node_offset]), \
                    &(grid_dv->vv_T0r_0[i_node_offset]), \
                    &(grid_dv->vv_T0t_0[i_node_offset]), \
                    &(grid_dv->vv_T0p_0[i_node_offset]), \
                    &(grid_dv->vv_fun_0[i_node_offset]), \
                    &(grid_dv->vv_change_0[i_node_offset]), \
                    &(grid_dv->tau), \
                    &grid_dv->loc_I_host, \
                    &grid_dv->loc_J_host, \
                    &grid_dv->loc_K_host, \
                    &grid_dv->dr_host, \
                    &grid_dv->dt_host, \
                    &grid_dv->dp_host, \
                    &n_nodes_this_level \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_1st, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30000);

            } else if (iswp == 1){
                void* kernelArgs[]= {\
                    &(grid_dv->vv_i__j__k___1[i_node_offset]), \
                    &(grid_dv->vv_i__jp1k___1[i_node_offset]), \
                    &(grid_dv->vv_i__jm1k___1[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp1_1[i_node_offset]), \
                    &(grid_dv->vv_i__j__km1_1[i_node_offset]), \
                    &(grid_dv->vv_ip1j__k___1[i_node_offset]), \
                    &(grid_dv->vv_im1j__k___1[i_node_offset]), \
                    &(grid_dv->vv_fac_a_1[i_node_offset]), \
                    &(grid_dv->vv_fac_b_1[i_node_offset]), \
                    &(grid_dv->vv_fac_c_1[i_node_offset]), \
                    &(grid_dv->vv_fac_f_1[i_node_offset]), \
                    &(grid_dv->vv_T0v_1[i_node_offset]), \
                    &(grid_dv->vv_T0r_1[i_node_offset]), \
                    &(grid_dv->vv_T0t_1[i_node_offset]), \
                    &(grid_dv->vv_T0p_1[i_node_offset]), \
                    &(grid_dv->vv_fun_1[i_node_offset]), \
                    &(grid_dv->vv_change_1[i_node_offset]), \
                    &(grid_dv->tau), \
                    &grid_dv->loc_I_host, \
                    &grid_dv->loc_J_host, \
                    &grid_dv->loc_K_host, \
                    &grid_dv->dr_host, \
                    &grid_dv->dt_host, \
                    &grid_dv->dp_host, \
                    &n_nodes_this_level \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_1st, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30001);

            } else if (iswp == 2){
                void* kernelArgs[]= {\
                    &(grid_dv->vv_i__j__k___2[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp1_2[i_node_offset]), \
                    &(grid_dv->vv_i__j__km1_2[i_node_offset]), \
                    &(grid_dv->vv_ip1j__k___2[i_node_offset]), \
                    &(grid_dv->vv_im1j__k___2[i_node_offset]), \
                    &(grid_dv->vv_i__jp1k___2[i_node_offset]), \
                    &(grid_dv->vv_i__jm1k___2[i_node_offset]), \
                    &(grid_dv->vv_fac_a_2[i_node_offset]), \
                    &(grid_dv->vv_fac_b_2[i_node_offset]), \
                    &(grid_dv->vv_fac_c_2[i_node_offset]), \
                    &(grid_dv->vv_fac_f_2[i_node_offset]), \
                    &(grid_dv->vv_T0v_2[i_node_offset]), \
                    &(grid_dv->vv_T0r_2[i_node_offset]), \
                    &(grid_dv->vv_T0t_2[i_node_offset]), \
                    &(grid_dv->vv_T0p_2[i_node_offset]), \
                    &(grid_dv->vv_fun_2[i_node_offset]), \
                    &(grid_dv->vv_change_2[i_node_offset]), \
                    &(grid_dv->tau), \
                    &grid_dv->loc_I_host, \
                    &grid_dv->loc_J_host, \
                    &grid_dv->loc_K_host, \
                    &grid_dv->dr_host, \
                    &grid_dv->dt_host, \
                    &grid_dv->dp_host, \
                    &n_nodes_this_level \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_1st, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30002);

            } else if (iswp == 3){
                void* kernelArgs[]= {\
                    &(grid_dv->vv_i__j__k___3[i_node_offset]), \
                    &(grid_dv->vv_ip1j__k___3[i_node_offset]), \
                    &(grid_dv->vv_im1j__k___3[i_node_offset]), \
                    &(grid_dv->vv_i__jp1k___3[i_node_offset]), \
                    &(grid_dv->vv_i__jm1k___3[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp1_3[i_node_offset]), \
                    &(grid_dv->vv_i__j__km1_3[i_node_offset]), \
                    &(grid_dv->vv_fac_a_3[i_node_offset]), \
                    &(grid_dv->vv_fac_b_3[i_node_offset]), \
                    &(grid_dv->vv_fac_c_3[i_node_offset]), \
                    &(grid_dv->vv_fac_f_3[i_node_offset]), \
                    &(grid_dv->vv_T0v_3[i_node_offset]), \
                    &(grid_dv->vv_T0r_3[i_node_offset]), \
                    &(grid_dv->vv_T0t_3[i_node_offset]), \
                    &(grid_dv->vv_T0p_3[i_node_offset]), \
                    &(grid_dv->vv_fun_3[i_node_offset]), \
                    &(grid_dv->vv_change_3[i_node_offset]), \
                    &(grid_dv->tau), \
                    &grid_dv->loc_I_host, \
                    &grid_dv->loc_J_host, \
                    &grid_dv->loc_K_host, \
                    &grid_dv->dr_host, \
                    &grid_dv->dt_host, \
                    &grid_dv->dp_host, \
                    &n_nodes_this_level \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_1st, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30003);

            } else if (iswp == 4){
                void* kernelArgs[]= {\
                    &(grid_dv->vv_i__j__k___4[i_node_offset]), \
                    &(grid_dv->vv_ip1j__k___4[i_node_offset]), \
                    &(grid_dv->vv_im1j__k___4[i_node_offset]), \
                    &(grid_dv->vv_i__jp1k___4[i_node_offset]), \
                    &(grid_dv->vv_i__jm1k___4[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp1_4[i_node_offset]), \
                    &(grid_dv->vv_i__j__km1_4[i_node_offset]), \
                    &(grid_dv->vv_fac_a_4[i_node_offset]), \
                    &(grid_dv->vv_fac_b_4[i_node_offset]), \
                    &(grid_dv->vv_fac_c_4[i_node_offset]), \
                    &(grid_dv->vv_fac_f_4[i_node_offset]), \
                    &(grid_dv->vv_T0v_4[i_node_offset]), \
                    &(grid_dv->vv_T0r_4[i_node_offset]), \
                    &(grid_dv->vv_T0t_4[i_node_offset]), \
                    &(grid_dv->vv_T0p_4[i_node_offset]), \
                    &(grid_dv->vv_fun_4[i_node_offset]), \
                    &(grid_dv->vv_change_4[i_node_offset]), \
                    &(grid_dv->tau), \
                    &grid_dv->loc_I_host, \
                    &grid_dv->loc_J_host, \
                    &grid_dv->loc_K_host, \
                    &grid_dv->dr_host, \
                    &grid_dv->dt_host, \
                    &grid_dv->dp_host, \
                    &n_nodes_this_level \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_1st, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30004);

            } else if (iswp == 5) {
                void* kernelArgs[]= {\
                    &(grid_dv->vv_i__j__k___5[i_node_offset]), \
                    &(grid_dv->vv_ip1j__k___5[i_node_offset]), \
                    &(grid_dv->vv_im1j__k___5[i_node_offset]), \
                    &(grid_dv->vv_i__jp1k___5[i_node_offset]), \
                    &(grid_dv->vv_i__jm1k___5[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp1_5[i_node_offset]), \
                    &(grid_dv->vv_i__j__km1_5[i_node_offset]), \
                    &(grid_dv->vv_fac_a_5[i_node_offset]), \
                    &(grid_dv->vv_fac_b_5[i_node_offset]), \
                    &(grid_dv->vv_fac_c_5[i_node_offset]), \
                    &(grid_dv->vv_fac_f_5[i_node_offset]), \
                    &(grid_dv->vv_T0v_5[i_node_offset]), \
                    &(grid_dv->vv_T0r_5[i_node_offset]), \
                    &(grid_dv->vv_T0t_5[i_node_offset]), \
                    &(grid_dv->vv_T0p_5[i_node_offset]), \
                    &(grid_dv->vv_fun_5[i_node_offset]), \
                    &(grid_dv->vv_change_5[i_node_offset]), \
                    &(grid_dv->tau), \
                    &grid_dv->loc_I_host, \
                    &grid_dv->loc_J_host, \
                    &grid_dv->loc_K_host, \
                    &grid_dv->dr_host, \
                    &grid_dv->dt_host, \
                    &grid_dv->dp_host, \
                    &n_nodes_this_level \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_1st, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30005);

            } else if (iswp == 6) {
                void* kernelArgs[]= {\
                    &(grid_dv->vv_i__j__k___6[i_node_offset]), \
                    &(grid_dv->vv_ip1j__k___6[i_node_offset]), \
                    &(grid_dv->vv_im1j__k___6[i_node_offset]), \
                    &(grid_dv->vv_i__jp1k___6[i_node_offset]), \
                    &(grid_dv->vv_i__jm1k___6[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp1_6[i_node_offset]), \
                    &(grid_dv->vv_i__j__km1_6[i_node_offset]), \
                    &(grid_dv->vv_fac_a_6[i_node_offset]), \
                    &(grid_dv->vv_fac_b_6[i_node_offset]), \
                    &(grid_dv->vv_fac_c_6[i_node_offset]), \
                    &(grid_dv->vv_fac_f_6[i_node_offset]), \
                    &(grid_dv->vv_T0v_6[i_node_offset]), \
                    &(grid_dv->vv_T0r_6[i_node_offset]), \
                    &(grid_dv->vv_T0t_6[i_node_offset]), \
                    &(grid_dv->vv_T0p_6[i_node_offset]), \
                    &(grid_dv->vv_fun_6[i_node_offset]), \
                    &(grid_dv->vv_change_6[i_node_offset]), \
                    &(grid_dv->tau), \
                    &grid_dv->loc_I_host, \
                    &grid_dv->loc_J_host, \
                    &grid_dv->loc_K_host, \
                    &grid_dv->dr_host, \
                    &grid_dv->dt_host, \
                    &grid_dv->dp_host, \
                    &n_nodes_this_level \
                };
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_1st, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30006);


            } else {
                void* kernelArgs[]= {\
                    &(grid_dv->vv_i__j__k___7[i_node_offset]), \
                    &(grid_dv->vv_ip1j__k___7[i_node_offset]), \
                    &(grid_dv->vv_im1j__k___7[i_node_offset]), \
                    &(grid_dv->vv_i__jp1k___7[i_node_offset]), \
                    &(grid_dv->vv_i__jm1k___7[i_node_offset]), \
                    &(grid_dv->vv_i__j__kp1_7[i_node_offset]), \
                    &(grid_dv->vv_i__j__km1_7[i_node_offset]), \
                    &(grid_dv->vv_fac_a_7[i_node_offset]), \
                    &(grid_dv->vv_fac_b_7[i_node_offset]), \
                    &(grid_dv->vv_fac_c_7[i_node_offset]), \
                    &(grid_dv->vv_fac_f_7[i_node_offset]), \
                    &(grid_dv->vv_T0v_7[i_node_offset]), \
                    &(grid_dv->vv_T0r_7[i_node_offset]), \
                    &(grid_dv->vv_T0t_7[i_node_offset]), \
                    &(grid_dv->vv_T0p_7[i_node_offset]), \
                    &(grid_dv->vv_fun_7[i_node_offset]), \
                    &(grid_dv->vv_change_7[i_node_offset]), \
                    &(grid_dv->tau), \
                    &grid_dv->loc_I_host, \
                    &grid_dv->loc_J_host, \
                    &grid_dv->loc_K_host, \
                    &grid_dv->dr_host, \
                    &grid_dv->dt_host, \
                    &grid_dv->dp_host, \
                    &n_nodes_this_level \
                };

                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_1st, grid_each, threads_each, kernelArgs, 0, grid_dv->level_streams[id_stream]), 30007);
            }
        }
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

        i_node_offset += grid_dv->n_nodes_on_levels_host[i_level];
    }
std::cout <<"aaa3"<<std::endl;

}