#include "iterator_wrapper.cuh"

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
    int& loc_I, \
    int& loc_J, \
    int& loc_K, \
    CUSTOMREAL& dr, \
    CUSTOMREAL& dt, \
    CUSTOMREAL& dp  \
){
    unsigned int i_node = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;

    if (i_node > loc_I*loc_J*loc_K) return;

    CUSTOMREAL sigr = 1.0*sqrt(fac_a[i__j__k__[i_node]])*T0v[i__j__k__[i_node]];
    CUSTOMREAL sigt = 1.0*sqrt(fac_b[i__j__k__[i_node]])*T0v[i__j__k__[i_node]];
    CUSTOMREAL sigp = 1.0*sqrt(fac_c[i__j__k__[i_node]])*T0v[i__j__k__[i_node]];
    CUSTOMREAL coe  = 1.0/((sigr/dr)+(sigt/dt)+(sigp/dp));

    CUSTOMREAL pp1 = (tau[i__j__k__[i_node]] - tau[im1j__k__[i_node]])/dp;
    CUSTOMREAL pp2 = (tau[ip1j__k__[i_node]] - tau[i__j__k__[i_node]])/dp;

    CUSTOMREAL pt1 = (tau[i__j__k__[i_node]] - tau[i__jm1k__[i_node]])/dt;
    CUSTOMREAL pt2 = (tau[i__jp1k__[i_node]] - tau[i__j__k__[i_node]])/dt;

    CUSTOMREAL pr1 = (tau[i__j__k__[i_node]] - tau[i__j__km1[i_node]])/dr;
    CUSTOMREAL pr2 = (tau[i__j__kp1[i_node]] - tau[i__j__k__[i_node]])/dr;

    // LF Hamiltonian
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
    int& loc_I, \
    int& loc_J, \
    int& loc_K, \
    CUSTOMREAL& dr, \
    CUSTOMREAL& dt, \
    CUSTOMREAL& dp  \
){

    CUSTOMREAL pp1, pp2, pt1, pt2, pr1, pr2;
    CUSTOMREAL wp1, wp2, wt1, wt2, wr1, wr2;
    const CUSTOMREAL eps = 1.0e-12;

    unsigned int i_node = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;

    if (i_node > loc_I*loc_J*loc_K) return;

    int k = i_node/(loc_I*loc_J);
    int j = (i_node - k*loc_I*loc_J)/loc_I;
    int i = i_node - k*loc_I*loc_J - j*loc_I;

    // direction p
    if (i == 1) {
        pp1 = (tau[i__j__k__[i_node]] \
             - tau[im1j__k__[i_node]]) / dp;

        wp2 = 1.0/(1.0+2.0*pow((eps + pow(        tau[i__j__k__] \
                                             -2.0*tau[ipj_k_] \
                                             +    tau[ip2j_k_],2.0) ) \
                               / (eps + pow(      tau[imj_k_] \
                                             -2.0*tau[i_j_k_] \
                                             +    tau[ipj_k_],2.0) ),2.0) );

        pp2 = (1.0 - wp2) * (         tau_loc[ipj_k_] \
                              -       tau_loc[imj_k_]) / 2.0 / dp \
                   + wp2  * ( - 3.0 * tau_loc[i_j_k_] \
                              + 4.0 * tau_loc[ipj_k_] \
                              -       tau_loc[ip2j_k_] ) / 2.0 / dp;

    } else if (i == loc_I-2) {
        wp1 = 1.0/(1.0+2.0*pow((eps + pow(            tau_loc[i_j_k_] \
                                                 -2.0*tau_loc[imj_k_] \
                                                 +    tau_loc[im2j_k_],2.0) ) \
                                   / (eps + pow(      tau_loc[ipj_k_] \
                                                 -2.0*tau_loc[i_j_k_] \
                                                 +    tau_loc[imj_k_],2.0) ),2.0) );

        pp1 = (1.0 - wp1) * (           tau_loc[ipj_k_] \
                                -       tau_loc[imj_k_]) / 2.0 / dp \
                     + wp1  * ( + 3.0 * tau_loc[i_j_k_] \
                                - 4.0 * tau_loc[imj_k_] \
                                +       tau_loc[im2j_k_] ) / 2.0 / dp;

        pp2 = (tau_loc[ipj_k_] - tau_loc[i_j_k_]) / dp;

    } else {
        wp1 = 1.0/(1.0+2.0*pow((eps + pow(        tau_loc[i_j_k_] \
                                             -2.0*tau_loc[imj_k_] \
                                             +    tau_loc[im2j_k_],2.0) ) \
                               / (eps + pow(      tau_loc[ipj_k_] \
                                             -2.0*tau_loc[i_j_k_] \
                                             +    tau_loc[imj_k_],2.0) ),2.0) );

        pp1 = (1.0 - wp1) * (         tau_loc[ipj_k_] \
                              -       tau_loc[imj_k_]) / 2.0 / dp \
                   + wp1  * ( + 3.0 * tau_loc[i_j_k_] \
                              - 4.0 * tau_loc[imj_k_] \
                              +       tau_loc[im2j_k_] ) / 2.0 / dp;

        wp2 = 1.0/(1.0+2.0*pow((eps + pow(        tau_loc[i_j_k_] \
                                             -2.0*tau_loc[ipj_k_] \
                                             +    tau_loc[ip2j_k_],2.0) ) \
                             / (eps + pow(        tau_loc[imj_k_] \
                                             -2.0*tau_loc[i_j_k_] \
                                             +    tau_loc[ipj_k_],2.0) ),2.0) );

        pp2 = (1.0 - wp2) * (           tau_loc[ipj_k_] \
                                -       tau_loc[imj_k_]) / 2.0 / dp \
                     + wp2  * ( - 3.0 * tau_loc[i_j_k_] \
                                + 4.0 * tau_loc[ipj_k_] \
                                -       tau_loc[ip2j_k_] ) / 2.0 / dp;

    }

    // direction t
    if (j == 1) {
        pt1 = (tau_loc[i_j_k_] \
             - tau_loc[i_jmk_]) / dt;

        wt2 = 1.0/(1.0+2.0*pow((eps + pow(        tau_loc[i_j_k_] \
                                             -2.0*tau_loc[i_jpk_] \
                                             +    tau_loc[i_jp2k_],2.0) ) \
                             / (eps + pow(        tau_loc[i_jmk_] \
                                             -2.0*tau_loc[i_j_k_] \
                                           +      tau_loc[i_jpk_],2.0) ),2.0) );

        pt2 = (1.0 - wt2) * (         tau_loc[i_jpk_] \
                            -         tau_loc[i_jmk_]) / 2.0 / dt \
                   + wt2  * ( - 3.0 * tau_loc[i_j_k_] \
                              + 4.0 * tau_loc[i_jpk_] \
                            -         tau_loc[i_jp2k_] ) / 2.0 / dt;

    } else if (j == loc_J-2) {
        wt1 = 1.0/(1.0+2.0*pow((eps + pow(        tau_loc[i_j_k_] \
                                             -2.0*tau_loc[i_jmk_] \
                                           +      tau_loc[i_jm2k_],2.0) ) \
                             / (eps + pow(        tau_loc[i_jpk_] \
                                             -2.0*tau_loc[i_j_k_] \
                                           +      tau_loc[i_jmk_],2.0) ),2.0) );

        pt1 = (1.0 - wt1) * (           tau_loc[i_jpk_] \
                                -       tau_loc[i_jmk_]) / 2.0 / dt \
                     + wt1  * ( + 3.0 * tau_loc[i_j_k_] \
                                - 4.0 * tau_loc[i_jmk_] \
                                +       tau_loc[i_jm2k_] ) / 2.0 / dt;

        pt2 = (tau_loc[i_jpk_] - tau_loc[i_j_k_]) / dt;

    } else {
        wt1 = 1.0/(1.0+2.0*pow((eps + pow(        tau_loc[i_j_k_] \
                                             -2.0*tau_loc[i_jmk_] \
                                             +    tau_loc[i_jm2k_],2.0) ) \
                             / (eps + pow(        tau_loc[i_jpk_] \
                                             -2.0*tau_loc[i_j_k_] \
                                             +    tau_loc[i_jmk_],2.0) ),2.0) );

        pt1 = (1.0 - wt1) * (           tau_loc[i_jpk_] \
                                -       tau_loc[i_jmk_]) / 2.0 / dt \
                     + wt1  * ( + 3.0 * tau_loc[i_j_k_] \
                                - 4.0 * tau_loc[i_jmk_] \
                                +       tau_loc[i_jm2k_] ) / 2.0 / dt;

        wt2 = 1.0/(1.0+2.0*pow((eps + pow(        tau_loc[i_j_k_] \
                                             -2.0*tau_loc[i_jpk_] \
                                           +      tau_loc[i_jp2k_],2.0) ) \
                             / (eps + pow(        tau_loc[i_jmk_] \
                                             -2.0*tau_loc[i_j_k_] \
                                           +      tau_loc[i_jpk_],2.0) ),2.0) );

        pt2 = (1.0 - wt2) * (           tau_loc[i_jpk_] \
                                -       tau_loc[i_jmk_]) / 2.0 / dt \
                     + wt2  * ( - 3.0 * tau_loc[i_j_k_] \
                                + 4.0 * tau_loc[i_jpk_] \
                                -       tau_loc[i_jp2k_] ) / 2.0 / dt;

    }

    // direction r
    if (k == 1) {
        pr1 = (tau_loc[i_j_k_] \
             - tau_loc[i_j_km]) / dr;

        wr2 = 1.0/(1.0+2.0*pow((eps + pow(        tau_loc[i_j_k_] \
                                             -2.0*tau_loc[i_j_kp] \
                                           +      tau_loc[i_j_kp2],2.0) ) \
                             / (eps + pow(        tau_loc[i_j_km] \
                                             -2.0*tau_loc[i_j_k_] \
                                           +      tau_loc[i_j_kp],2.0) ),2.0) );

        pr2 = (1.0 - wr2) * (           tau_loc[i_j_kp] \
                              -         tau_loc[i_j_km]) / 2.0 / dr \
                     + wr2  * ( - 3.0 * tau_loc[i_j_k_] \
                                + 4.0 * tau_loc[i_j_kp] \
                              -         tau_loc[i_j_kp2] ) / 2.0 / dr;

    } else if (k == loc_K - 2) {
        wr1 = 1.0/(1.0+2.0*pow((eps + pow(        tau_loc[i_j_k_] \
                                             -2.0*tau_loc[i_j_km] \
                                           +      tau_loc[i_j_km2],2.0) ) \
                             / (eps + pow(        tau_loc[i_j_kp] \
                                             -2.0*tau_loc[i_j_k_] \
                                           +      tau_loc[i_j_km],2.0) ),2.0) );

        pr1 = (1.0 - wr1) * (            tau_loc[i_j_kp] \
                              -          tau_loc[i_j_km]) / 2.0 / dr \
                      + wr1  * ( + 3.0 * tau_loc[i_j_k_] \
                                 - 4.0 * tau_loc[i_j_km] \
                               +         tau_loc[i_j_km2] ) / 2.0 / dr;

        pr2 = (tau_loc[i_j_kp] \
             - tau_loc[i_j_k_]) / dr;

    } else {
        wr1 = 1.0/(1.0+2.0*pow((eps + pow(        tau_loc[i_j_k_] \
                                             -2.0*tau_loc[i_j_km] \
                                           +      tau_loc[i_j_km2],2.0) ) \
                             / (eps + pow(        tau_loc[i_j_kp] \
                                             -2.0*tau_loc[i_j_k_] \
                                           +      tau_loc[i_j_km],2.0) ),2.0) );

        pr1 = (1.0 - wr1) * (           tau_loc[i_j_kp] \
                              -         tau_loc[i_j_km]) / 2.0 / dr \
                     + wr1  * ( + 3.0 * tau_loc[i_j_k_] \
                                - 4.0 * tau_loc[i_j_km] \
                              +         tau_loc[i_j_km2] ) / 2.0 / dr;

        wr2 = 1.0/(1.0+2.0*pow((eps + pow(        tau_loc[i_j_k_] \
                                             -2.0*tau_loc[i_j_kp] \
                                           +      tau_loc[i_j_kp2],2.0) ) \
                             / (eps + pow(        tau_loc[i_j_km] \
                                             -2.0*tau_loc[i_j_k_] \
                                           +      tau_loc[i_j_kp],2.0) ),2.0) );

        pr2 = (1.0 - wr2) * (          tau_loc[i_j_kp] \
                              -        tau_loc[i_j_km]) / 2.0 / dr \
                     + wr2  * ( - 3.0 *tau_loc[i_j_k_] \
                                + 4.0 *tau_loc[i_j_kp] \
                              -        tau_loc[i_j_kp2] ) / 2.0 / dr;

    }



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


// this function calculate all levels of one single sweep direction
void cuda_run_iterate_forward(Grid_on_device* grid_dv, int const& iswp){

    int block_size = CUDA_SWEEPING_BLOCK_SIZE;
    int num_blocks_x, num_blocks_y;
    int actual_end_level = grid_dv->n_levels_host;

    for (size_t i_level = 0; i_level < actual_end_level; i_level++){
        get_block_xy(ceil(grid_dv->n_nodes_on_levels[i_level]/block_size+0.5), &num_blocks_x, &num_blocks_y);
        dim3 grid_each(num_blocks_x, num_blocks_y);
        dim3 block_each(block_size, 1, 1);

        void* kernelArgs[]= {\


        }


    }

}