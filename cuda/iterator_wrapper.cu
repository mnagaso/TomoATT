#include "iterator_wrapper.cuh"


//__device__ void warpReduce(volatile CUSTOMREAL *sdata, unsigned int tid) {
//    if (CUDA_L1_BLOCK_SIZE >= 64) sdata[tid] += sdata[tid + 32];
//    if (CUDA_L1_BLOCK_SIZE >= 32) sdata[tid] += sdata[tid + 16];
//    if (CUDA_L1_BLOCK_SIZE >= 16) sdata[tid] += sdata[tid + 8];
//    if (CUDA_L1_BLOCK_SIZE >= 8) sdata[tid] += sdata[tid + 4];
//    if (CUDA_L1_BLOCK_SIZE >= 4) sdata[tid] += sdata[tid + 2];
//    if (CUDA_L1_BLOCK_SIZE >= 2) sdata[tid] += sdata[tid + 1];
//}
//
//__global__ void L1_reduce_kernel(CUSTOMREAL* g_idata, CUSTOMREAL* g_odata, int n)
//{
//    extern __shared__ CUSTOMREAL sdata[];
//    unsigned int tid = threadIdx.x;
//    unsigned int i = blockIdx.x*(CUDA_L1_BLOCK_SIZE*2) + tid;
//    unsigned int gridSize = CUDA_L1_BLOCK_SIZE*2*gridDim.x;
//    sdata[tid] = 0;
//    while (i < n) { sdata[tid] += g_idata[i] + g_idata[i+CUDA_L1_BLOCK_SIZE]; i += gridSize; }
//    __syncthreads();
//    if (CUDA_L1_BLOCK_SIZE >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
//    if (CUDA_L1_BLOCK_SIZE >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
//    if (CUDA_L1_BLOCK_SIZE >= 128) { if (tid < 64)  { sdata[tid] += sdata[tid + 64]; } __syncthreads(); }
//    if (tid < 32) warpReduce(sdata, tid);
//    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
//}

__global__ void L1_reduce_kernel(CUSTOMREAL* d_in, CUSTOMREAL* g_out, int n)
{

    extern __shared__ CUSTOMREAL diff_shared[]; // use shared memory accumulation for reduction
    //CUSTOMREAL *diff_shared = SharedMemory<CUSTOMREAL>();
    unsigned int i_block = blockIdx.x + blockIdx.y*gridDim.x;
    //unsigned int n_blocks = gridDim.x*gridDim.y;

    unsigned int ithread = threadIdx.x;
    // x y grid
    unsigned int iglobal = threadIdx.x + (i_block)*blockDim.x;

    if (iglobal >= n || ithread >= blockDim.x) return;

    //printf("%d %d %d\n", iglobal, ithread, blockDim.x);
    //printf("iglobal %d\n", iglobal);
    //printf("d_in[iglobal] %f\n", d_in[iglobal]);
    diff_shared[ithread] = d_in[iglobal];

    //printf("diff_shared[ithread] %f\n", diff_shared[ithread]);

    __syncthreads();


   // do reduction in shared mem
//    for(unsigned int s=1; s < blockDim.x; s *= 2) {
//        if (ithread % (2*s) == 0) {
//            diff_shared[ithread] += diff_shared[ithread + s];
//        }
//        __syncthreads();
//    }

    for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
        if (ithread < s) {
            diff_shared[ithread] += diff_shared[ithread + s];
        }
        __syncthreads();
    }

    // write result for this block to global mem
    if (ithread == 0) g_out[i_block] = diff_shared[0];

    __syncthreads();

}

__global__ void calculate_L1_kernel(CUSTOMREAL* tau_loc, CUSTOMREAL* tau_old_loc, CUSTOMREAL* T0v_loc, \
                                    CUSTOMREAL* L1_tmp_dev, \
                                    int loc_I, int loc_J, int loc_K, int n_ghost_layer, CUSTOMREAL* L1){


    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;


    // skip if the assigned ijk is greater than the grid size
    if (i >= loc_I-n_ghost_layer || j >= loc_J-n_ghost_layer || k >= loc_K-n_ghost_layer || i < n_ghost_layer || j < n_ghost_layer || k < n_ghost_layer){
        return;
        //tmp = 0.0;
        //unsigned int ijk = I2V_cuda(i,j,k, loc_I, loc_J);
        //printf("i j k ijk %d %d %d %d\n", i, j, k, ijk);

        //printf("ijk %d ithread %d\n", ijk, ithread);
        //L1_tmp_dev[ijk] = 0.0;

    }
    else{
        unsigned int ijk = I2V_cuda(i,j,k, loc_I, loc_J);
        //printf("i j k ijk %d %d %d %d\n", i, j, k, ijk);

        //printf("ijk %d ithread %d\n", ijk, ithread);
        L1_tmp_dev[ijk] = (CUSTOMREAL) fabs(tau_loc[ijk] - tau_old_loc[ijk])*T0v_loc[ijk];
    }

}


void cuda_calculate_L1(Grid_on_device* grid_on_dv){
    // calculate L1

    // set tau_old_loc to 0.0
    print_CUDA_error_if_any(cudaMemset(grid_on_dv->L1_tmp_dev, 0, grid_on_dv->loc_I_host*grid_on_dv->loc_J_host*grid_on_dv->loc_K_host*sizeof(CUSTOMREAL)), 44444);

    // calculate the initial L1 and Linf values

    calculate_L1_kernel<<<grid_on_dv->grid_3d_full_incr_host, grid_on_dv->threads_3d_full_incr_host>>>\
                                             (grid_on_dv->tau_loc_dev, grid_on_dv->tau_old_loc_dev, grid_on_dv->T0v_loc_dev, \
                                              grid_on_dv->L1_tmp_dev, \
                                              grid_on_dv->loc_I_host, grid_on_dv->loc_J_host, grid_on_dv->loc_K_host, \
                                              grid_on_dv->n_ghost_layers_host, \
                                              &(grid_on_dv->L1_host));
    //cudaDeviceSynchronize(); // not needed

    //printf("nblocks %d\n", nblocks);

    // initialize
    grid_on_dv->L1_host = 0.0;
    for (int i = 0; i < grid_on_dv->n_blocks_L1_host; i++){
        grid_on_dv->L1_arr_unified[i] = 0.0;
    }


    int n_total_pt = grid_on_dv->loc_I_host*grid_on_dv->loc_J_host*grid_on_dv->loc_K_host;

    L1_reduce_kernel<<<grid_on_dv->grid_L1_host, grid_on_dv->threads_L1_host, CUDA_L1_BLOCK_SIZE*sizeof(CUSTOMREAL)>>> \
                        (grid_on_dv->L1_tmp_dev, grid_on_dv->L1_arr_unified, n_total_pt);
    cudaDeviceSynchronize();

    // reduce among the blocks
    for (int i = 0; i < grid_on_dv->n_blocks_L1_host; ++i) {
        //printf("L1_arr_unified[%d] %f\n", i, grid_on_dv->L1_arr_unified[i]);
        grid_on_dv->L1_host += grid_on_dv->L1_arr_unified[i];
    }


    // copy L1_dev to L1_host
    //printf("L1: %3.16f\n", grid_on_dv->L1_host);

    // mpi reduce summation of L1
    CUSTOMREAL L1tmp=grid_on_dv->L1_host;
    MPI_Allreduce(&L1tmp, &(grid_on_dv->L1_host), 1, MPI_CR, MPI_SUM, grid_on_dv->inter_sub_comm_host);

    grid_on_dv->L1_host = grid_on_dv->L1_host/((grid_on_dv->ngrid_i_host-2) \
                                             * (grid_on_dv->ngrid_j_host-2) \
                                             * (grid_on_dv->ngrid_k_host-2));
    //printf("L1: %3.16f\n", grid_on_dv->L1_host);


}


__global__ void copy_tau_to_tau_old_kernel( \
    int loc_I, int loc_J, int loc_K, \
    CUSTOMREAL* tau_loc, CUSTOMREAL* tau_old_loc){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    // skip if the assigned ijk is greater than the grid size
    if (i >= loc_I || j >= loc_J || k >= loc_K)
        return;

    int ijk = I2V_cuda(i,j,k, loc_I, loc_J);
    tau_old_loc[ijk] = tau_loc[ijk];

}


void cuda_tau2old_tau(Grid_on_device* grid_on_dv){
    // copy tau to tau_old
    //copy_tau_to_tau_old_kernel<<<grid_on_dv->grid_3d_full_incr_host, grid_on_dv->threads_3d_full_incr_host>>>\
    //                                               (grid_on_dv->loc_I_host,grid_on_dv->loc_J_host, grid_on_dv->loc_K_host, \
    //                                                grid_on_dv->tau_loc_dev, grid_on_dv->tau_old_loc_dev);

    print_CUDA_error_if_any( cudaMemcpy(grid_on_dv->tau_old_loc_dev, grid_on_dv->tau_loc_dev, grid_on_dv->loc_I_host*grid_on_dv->loc_J_host*grid_on_dv->loc_K_host*sizeof(CUSTOMREAL), cudaMemcpyDeviceToDevice), 33333);
}



__device__ CUSTOMREAL cuda_calc_LF_Hamiltonian( \
                                            CUSTOMREAL* fac_a_loc, \
                                            CUSTOMREAL* fac_b_loc, \
                                            CUSTOMREAL* fac_c_loc, \
                                            CUSTOMREAL* fac_f_loc, \
                                            CUSTOMREAL* T0r_loc, \
                                            CUSTOMREAL* T0t_loc, \
                                            CUSTOMREAL* T0p_loc, \
                                            CUSTOMREAL* T0v_loc, \
                                            CUSTOMREAL* tau_loc, \
                                            CUSTOMREAL& pp1, CUSTOMREAL& pp2, \
                                            CUSTOMREAL& pt1, CUSTOMREAL& pt2, \
                                            CUSTOMREAL& pr1, CUSTOMREAL& pr2, \
                                            int& i_j_k_) {
    // LF Hamiltonian for T = T0 + tau
    /*
    return sqrt(
              grid.fac_a_loc[i_j_k_] * (pow((pr1+pr2)/2.0,2.0)) \
      +       grid.fac_b_loc[i_j_k_] * (pow((pt1+pt2)/2.0,2.0)) \
      +       grid.fac_c_loc[i_j_k_] * (pow((pp1+pp2)/2.0,2.0)) \
      - 2.0*grid.fac_f_loc[i_j_k_] * (pt1+pt2)/2.0 * (pp1+pp2)/2.0 \

      +       grid.fac_a_loc[i_j_k_] * (pr1+pr2) * grid.T0r_loc[i_j_k_] \
      +       grid.fac_b_loc[i_j_k_] * (pt1+pt2) * grid.T0t_loc[i_j_k_] - grid.fac_f_loc[i_j_k_] * (pt1+pt2) * grid.T0p_loc[i_j_k_] \
      +       grid.fac_c_loc[i_j_k_] * (pp1+pp2) * grid.T0p_loc[i_j_k_] - grid.fac_f_loc[i_j_k_] * (pp1+pp2) * grid.T0t_loc[i_j_k_] \

      +       grid.fac_a_loc[i_j_k_] * pow(grid.T0r_loc[i_j_k_],2.0) \
      +       grid.fac_b_loc[i_j_k_] * pow(grid.T0t_loc[i_j_k_],2.0) \
      +       grid.fac_c_loc[i_j_k_] * pow(grid.T0p_loc[i_j_k_],2.0) \
      - 2.0*grid.fac_f_loc[i_j_k_] * grid.T0t_loc[i_j_k_] * grid.T0p_loc[i_j_k_] \
    );
    */

    // LF Hamiltonian for T = T0 * tau
    return sqrt(
              fac_a_loc[i_j_k_] * (T0r_loc[i_j_k_] * tau_loc[i_j_k_] + T0v_loc[i_j_k_] * (pr1+pr2)/2.0) * (T0r_loc[i_j_k_] * tau_loc[i_j_k_] + T0v_loc[i_j_k_] * (pr1+pr2)/2.0)\
    +         fac_b_loc[i_j_k_] * (T0t_loc[i_j_k_] * tau_loc[i_j_k_] + T0v_loc[i_j_k_] * (pt1+pt2)/2.0) * (T0t_loc[i_j_k_] * tau_loc[i_j_k_] + T0v_loc[i_j_k_] * (pt1+pt2)/2.0)\
    +         fac_c_loc[i_j_k_] * (T0p_loc[i_j_k_] * tau_loc[i_j_k_] + T0v_loc[i_j_k_] * (pp1+pp2)/2.0) * (T0p_loc[i_j_k_] * tau_loc[i_j_k_] + T0v_loc[i_j_k_] * (pp1+pp2)/2.0)\
    -     2.0*fac_f_loc[i_j_k_] * (T0t_loc[i_j_k_] * tau_loc[i_j_k_] + T0v_loc[i_j_k_] * (pt1+pt2)/2.0) \
                                * (T0p_loc[i_j_k_] * tau_loc[i_j_k_] + T0v_loc[i_j_k_] * (pp1+pp2)/2.0) \
    );


}


__device__ void cuda_calculate_stencil_1st_order( \
                                            CUSTOMREAL* fac_a_loc, \
                                            CUSTOMREAL* fac_b_loc, \
                                            CUSTOMREAL* fac_c_loc, \
                                            CUSTOMREAL* fac_f_loc, \
                                            CUSTOMREAL* T0r_loc, \
                                            CUSTOMREAL* T0t_loc, \
                                            CUSTOMREAL* T0p_loc, \
                                            CUSTOMREAL* T0v_loc, \
                                            CUSTOMREAL* tau_loc, \
                                            CUSTOMREAL* fun_loc, \
                                            CUSTOMREAL dr, CUSTOMREAL dt, CUSTOMREAL dp, \
                                            int i_j_k_, \
                                            int ipj_k_, \
                                            int imj_k_, \
                                            int i_jpk_, \
                                            int i_jmk_, \
                                            int i_j_kp, \
                                            int i_j_km){

    CUSTOMREAL sigr = 1.0*sqrt(fac_a_loc[i_j_k_])*T0v_loc[i_j_k_];
    CUSTOMREAL sigt = 1.0*sqrt(fac_b_loc[i_j_k_])*T0v_loc[i_j_k_];
    CUSTOMREAL sigp = 1.0*sqrt(fac_c_loc[i_j_k_])*T0v_loc[i_j_k_];
    CUSTOMREAL coe  = 1.0/((sigr/dr)+(sigt/dt)+(sigp/dp));

    //printf("i_j_k_ , imj_k_ , ipj_k_ , i_jpk_ , i_jmk_ , i_j_kp , i_j_km, dp, dt, dr :  %d , %d , %d , %d , %d , %d , %d , %f , %f, %f \n", i_j_k_, imj_k_, ipj_k_, i_jpk_, i_jmk_, i_j_kp, i_j_km, dp, dt, dr);

    CUSTOMREAL pp1 = (tau_loc[i_j_k_] - tau_loc[imj_k_])/dp;
    CUSTOMREAL pp2 = (tau_loc[ipj_k_] - tau_loc[i_j_k_])/dp;

    CUSTOMREAL pt1 = (tau_loc[i_j_k_] - tau_loc[i_jmk_])/dt;
    CUSTOMREAL pt2 = (tau_loc[i_jpk_] - tau_loc[i_j_k_])/dt;

    CUSTOMREAL pr1 = (tau_loc[i_j_k_] - tau_loc[i_j_km])/dr;
    CUSTOMREAL pr2 = (tau_loc[i_j_kp] - tau_loc[i_j_k_])/dr;

    // LF Hamiltonian
    CUSTOMREAL Htau = cuda_calc_LF_Hamiltonian(\
                                               fac_a_loc, \
                                               fac_b_loc, \
                                               fac_c_loc, \
                                               fac_f_loc, \
                                               T0r_loc, \
                                               T0t_loc, \
                                               T0p_loc, \
                                               T0v_loc, \
                                               tau_loc, \
                                               pp1, pp2, pt1, pt2, pr1, pr2, i_j_k_);

    tau_loc[i_j_k_] += coe*(fun_loc[i_j_k_] - Htau) \
                     + coe*(sigr*(pr2-pr1)/2.0 + sigt*(pt2-pt1)/2.0 + sigp*(pp2-pp1)/2.0);

}


__device__ void cuda_calculate_stencil_3rd_order( \
                                            CUSTOMREAL* fac_a_loc, \
                                            CUSTOMREAL* fac_b_loc, \
                                            CUSTOMREAL* fac_c_loc, \
                                            CUSTOMREAL* fac_f_loc, \
                                            CUSTOMREAL* T0r_loc, \
                                            CUSTOMREAL* T0t_loc, \
                                            CUSTOMREAL* T0p_loc, \
                                            CUSTOMREAL* T0v_loc, \
                                            CUSTOMREAL* tau_loc, \
                                            CUSTOMREAL* fun_loc, \
                                            CUSTOMREAL& dr, CUSTOMREAL& dt, CUSTOMREAL& dp, \
                                            int& i_j_k_, \
                                            int& ipj_k_, \
                                            int& imj_k_, \
                                            int& i_jpk_, \
                                            int& i_jmk_, \
                                            int& i_j_kp, \
                                            int& i_j_km, \
                                            int& ip2j_k_, \
                                            int& im2j_k_, \
                                            int& i_jp2k_, \
                                            int& i_jm2k_, \
                                            int& i_j_kp2, \
                                            int& i_j_km2, \
                                            int& i, int& j, int& k, \
                                            int& loc_I, int& loc_J, int& loc_K){

    CUSTOMREAL sigr = 1.0*sqrt(fac_a_loc[i_j_k_])*T0v_loc[i_j_k_];
    CUSTOMREAL sigt = 1.0*sqrt(fac_b_loc[i_j_k_])*T0v_loc[i_j_k_];
    CUSTOMREAL sigp = 1.0*sqrt(fac_c_loc[i_j_k_])*T0v_loc[i_j_k_];
    CUSTOMREAL coe  = 1.0/((sigr/dr)+(sigt/dt)+(sigp/dp));

    CUSTOMREAL pp1, pp2, pt1, pt2, pr1, pr2;
    CUSTOMREAL wp1, wp2, wt1, wt2, wr1, wr2;

    const CUSTOMREAL eps = 1.0e-12;

    // direction p
    if (i == 1) {
        pp1 = (tau_loc[i_j_k_] \
             - tau_loc[imj_k_]) / dp;

        wp2 = 1.0/(1.0+2.0*pow((eps + pow(        tau_loc[i_j_k_] \
                                             -2.0*tau_loc[ipj_k_] \
                                             +    tau_loc[ip2j_k_],2.0) ) \
                               / (eps + pow(      tau_loc[imj_k_] \
                                             -2.0*tau_loc[i_j_k_] \
                                             +    tau_loc[ipj_k_],2.0) ),2.0) );

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

    // LF Hamiltonian
    CUSTOMREAL Htau = cuda_calc_LF_Hamiltonian(\
                                               fac_a_loc, \
                                               fac_b_loc, \
                                               fac_c_loc, \
                                               fac_f_loc, \
                                               T0r_loc, \
                                               T0t_loc, \
                                               T0p_loc, \
                                               T0v_loc, \
                                               tau_loc, \
                                               pp1, pp2, pt1, pt2, pr1, pr2, i_j_k_);


    // update tau
    tau_loc[i_j_k_] += coe * (fun_loc[i_j_k_] - Htau) \
                     + coe * (sigr*(pr2-pr1)/2.0 + sigt*(pt2-pt1)/2.0 + sigp*(pp2-pp1)/2.0);


}


__global__ void cuda_do_sweep_level_kernel_1st(\
                                           int i_level_in, \
                                           int id_start_in, \
                                           int actual_end_level, \
                                           int* n_points_each_level, \
                                           int r_dirc, int t_dirc, int p_dirc, \
                                           int* i_id_all_level, \
                                           int* j_id_all_level, \
                                           int* k_id_all_level, \
                                           int loc_I, int loc_J, int loc_K, \
                                           int stencil_order, \
                                           CUSTOMREAL* fac_a_loc, \
                                           CUSTOMREAL* fac_b_loc, \
                                           CUSTOMREAL* fac_c_loc, \
                                           CUSTOMREAL* fun_loc, \
                                           CUSTOMREAL* T0v_loc, \
                                           CUSTOMREAL* tau_loc, \
                                           CUSTOMREAL dr, CUSTOMREAL dt, CUSTOMREAL dp, \
                                           CUSTOMREAL* fac_f_loc, \
                                           CUSTOMREAL* T0r_loc, \
                                           CUSTOMREAL* T0t_loc, \
                                           CUSTOMREAL* T0p_loc, \
                                           bool* is_changed){


    // prepare for grid synchronization (60 Pascal or later version of CUDA device is required)
    cooperative_groups::grid_group   grid = cooperative_groups::this_grid(); // use for grid synchronization
    cooperative_groups::thread_block block = cooperative_groups::this_thread_block(); // use for thread synchronization

    int id_start = 0; // offset
    int i_level_start = 0;
    bool unified_kernel = true;

    if (i_level_in != 9999) {
        actual_end_level = i_level_in;
        i_level_start = i_level_in;
        id_start = id_start_in;
        unified_kernel = false;
    }

    for (int i_level = i_level_start; i_level <= actual_end_level; i_level++) {

        if (unified_kernel) {
            block.sync();
            grid.sync();
        }

        int n_point_this_level = n_points_each_level[i_level];

        int id_point = (blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x + threadIdx.x;

        //grid.sync(); // synchronize grid
        //printf("i_level id_point %d %d\n", i_level, id_point);
        // continue if the point is out of the range of this level
        if (id_point < n_point_this_level) {

            int i = i_id_all_level[id_start+id_point];
            int j = j_id_all_level[id_start+id_point];
            int k = k_id_all_level[id_start+id_point];

            //printf("id_point, id_start  , i, j, k, loc_I, loc_J, loc_K : %d %d %d, %d, %d, %d, %d, %d \n", id_point,id_start,i, j, k, loc_I, loc_J, loc_K);
            //if (i >= loc_I-1|| j >= loc_J-1 || k >= loc_K-1 || i <= 1 || j <= 1 || k <= 1) return;

            if (r_dirc < 0) k = loc_K-k; //kk-1;
            else            k = k-1;  //nr-kk;
            if (t_dirc < 0) j = loc_J-j; //jj-1;
            else            j = j-1;  //nt-jj;
            if (p_dirc < 0) i = loc_I-i; //ii-1;
            else            i = i-1; //np-ii;

            // indices
            int i_j_k_ = I2V_cuda(i, j, k, loc_I, loc_J);
            int ipj_k_ = I2V_cuda(i+1,j,k, loc_I, loc_J);
            int imj_k_ = I2V_cuda(i-1,j,k, loc_I, loc_J);
            int i_jpk_ = I2V_cuda(i,j+1,k, loc_I, loc_J);
            int i_jmk_ = I2V_cuda(i,j-1,k, loc_I, loc_J);
            int i_j_kp = I2V_cuda(i,j,k+1, loc_I, loc_J);
            int i_j_km = I2V_cuda(i,j,k-1, loc_I, loc_J);

            // return if is_changed is false
            if (is_changed[i_j_k_]) {
                // calculate stencil
                // first order
                cuda_calculate_stencil_1st_order( \
                                                fac_a_loc, \
                                                fac_b_loc, \
                                                fac_c_loc, \
                                                fac_f_loc, \
                                                T0r_loc, \
                                                T0t_loc, \
                                                T0p_loc, \
                                                T0v_loc, \
                                                tau_loc, \
                                                fun_loc, \
                                                dr, dt, dp, \
                i_j_k_, ipj_k_, imj_k_, i_jpk_, i_jmk_, i_j_kp, i_j_km);

            } // end if is_changed

            if(unified_kernel) id_start += n_point_this_level;

        } // end if id_point < n_point_this_level
        // synchronize threads
        // synchronize grids

    } // end for i_level


}


__global__ void cuda_do_sweep_level_kernel_3rd(\
                                           int i_level_in, \
                                           int id_start_in, \
                                           int actual_end_level, \
                                           int* n_points_each_level, \
                                           int r_dirc, int t_dirc, int p_dirc, \
                                           int* i_id_all_level, \
                                           int* j_id_all_level, \
                                           int* k_id_all_level, \
                                           int loc_I, int loc_J, int loc_K, \
                                           int stencil_order, \
                                           CUSTOMREAL* fac_a_loc, \
                                           CUSTOMREAL* fac_b_loc, \
                                           CUSTOMREAL* fac_c_loc, \
                                           CUSTOMREAL* fun_loc, \
                                           CUSTOMREAL* T0v_loc, \
                                           CUSTOMREAL* tau_loc, \
                                           CUSTOMREAL dr, CUSTOMREAL dt, CUSTOMREAL dp, \
                                           CUSTOMREAL* fac_f_loc, \
                                           CUSTOMREAL* T0r_loc, \
                                           CUSTOMREAL* T0t_loc, \
                                           CUSTOMREAL* T0p_loc, \
                                           bool* is_changed){


    // prepare for grid synchronization (60 Pascal or later version of CUDA device is required)
    cooperative_groups::grid_group        grid  = cooperative_groups::this_grid(); // use for grid synchronization
    cooperative_groups::thread_block      block = cooperative_groups::this_thread_block(); // use for thread synchronization
    //cooperative_groups::coalesced_group active = cooperative_groups::coalesced_threads();

    int id_start = 0; // offset
    int i_level_start = 0;
    bool unified_kernel = true;

    if (i_level_in != 9999) {
        actual_end_level = i_level_in;
        i_level_start = i_level_in;
        id_start = id_start_in;
        unified_kernel = false;
    }

    for (int i_level = i_level_start; i_level <= actual_end_level; i_level++) {

        if (unified_kernel) {
            block.sync();
            grid.sync();
        }


        int n_point_this_level = n_points_each_level[i_level];

        int id_point = (blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x + threadIdx.x;

        //printf("i_level id_point %d %d\n", i_level, id_point);
        // continue if the point is out of the range of this level
        if (id_point < n_point_this_level) {

            int i = i_id_all_level[id_start+id_point];
            int j = j_id_all_level[id_start+id_point];
            int k = k_id_all_level[id_start+id_point];

            //printf("id_point, id_start  , i, j, k, loc_I, loc_J, loc_K : %d %d %d, %d, %d, %d, %d, %d \n", id_point,id_start,i, j, k, loc_I, loc_J, loc_K);
            //if (i >= loc_I-1|| j >= loc_J-1 || k >= loc_K-1 || i <= 1 || j <= 1 || k <= 1) continue;
            //if (i >= loc_I - 1|| j >= loc_J -1|| k >= loc_K -1 || i <= 0 || j <= 0 || k <= 0) continue;

            if (r_dirc < 0) k = loc_K-k; //kk-1;
            else            k = k-1;  //nr-kk;
            if (t_dirc < 0) j = loc_J-j; //jj-1;
            else            j = j-1;  //nt-jj;
            if (p_dirc < 0) i = loc_I-i; //ii-1;
            else            i = i-1; //np-ii;

            // indices
            int i_j_k_ = I2V_cuda(i, j, k, loc_I, loc_J);
            int ipj_k_ = I2V_cuda(i+1,j,k, loc_I, loc_J);
            int imj_k_ = I2V_cuda(i-1,j,k, loc_I, loc_J);
            int i_jpk_ = I2V_cuda(i,j+1,k, loc_I, loc_J);
            int i_jmk_ = I2V_cuda(i,j-1,k, loc_I, loc_J);
            int i_j_kp = I2V_cuda(i,j,k+1, loc_I, loc_J);
            int i_j_km = I2V_cuda(i,j,k-1, loc_I, loc_J);
            int ip2j_k_ = I2V_cuda(i+2,j,k, loc_I, loc_J);
            int im2j_k_ = I2V_cuda(i-2,j,k, loc_I, loc_J);
            int i_jp2k_ = I2V_cuda(i,j+2,k, loc_I, loc_J);
            int i_jm2k_ = I2V_cuda(i,j-2,k, loc_I, loc_J);
            int i_j_kp2 = I2V_cuda(i,j,k+2, loc_I, loc_J);
            int i_j_km2 = I2V_cuda(i,j,k-2, loc_I, loc_J);


            // return if is_changed is false
            if (is_changed[i_j_k_]) {
                // calculate stencil
                // third order
                cuda_calculate_stencil_3rd_order( \
                                                fac_a_loc, \
                                                fac_b_loc, \
                                                fac_c_loc, \
                                                fac_f_loc, \
                                                T0r_loc, \
                                                T0t_loc, \
                                                T0p_loc, \
                                                T0v_loc, \
                                                tau_loc, \
                                                fun_loc, \
                                                dr, dt, dp, \
                i_j_k_, ipj_k_, imj_k_, i_jpk_, i_jmk_, i_j_kp, i_j_km, \
                ip2j_k_, im2j_k_, i_jp2k_, i_jm2k_, i_j_kp2, i_j_km2, \
                                                i,j,k, \
                                                loc_I, loc_J, loc_K);
            } // end if is_changed

            if(unified_kernel) id_start += n_point_this_level;

        } // end if id_point < n_point_this_level

        // synchronize threads
        // synchronize grids
        //if (unified_kernel) {
        //    block.sync();
        //    grid.sync();
        //}

    } // end for i_level


}




__global__ void cuda_calc_boundary_nodes_kernel_ij( \
    int loc_I, int loc_J, int loc_K, \
    CUSTOMREAL* tau_loc) {

    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int j = blockIdx.y*blockDim.y+threadIdx.y;

    if (i >= loc_I || j >= loc_J || i < 0 || j < 0) return;

    int i_j_0_ = I2V_cuda(i, j, 0, loc_I, loc_J);
    int i_j_1_ = I2V_cuda(i, j, 1, loc_I, loc_J);
    int i_j_2_ = I2V_cuda(i, j, 2, loc_I, loc_J);
    int i_j_n1 = I2V_cuda(i, j, loc_K-1, loc_I, loc_J);
    int i_j_n2 = I2V_cuda(i, j, loc_K-2, loc_I, loc_J);
    int i_j_n3 = I2V_cuda(i, j, loc_K-3, loc_I, loc_J);

    CUSTOMREAL v0 = 2.0 * tau_loc[i_j_1_] - tau_loc[i_j_2_];
    CUSTOMREAL v1 = tau_loc[i_j_2_];
    if (v0 > v1)
        tau_loc[i_j_0_] = v0;
    else
        tau_loc[i_j_0_] = v1;


    v0 = 2.0 * tau_loc[i_j_n2] - tau_loc[i_j_n3];
    v1 = tau_loc[i_j_n3];
    if (v0 > v1)
        tau_loc[i_j_n1] = v0;
    else
        tau_loc[i_j_n1] = v1;
}


__global__ void cuda_calc_boundary_nodes_kernel_ik( \
    int loc_I, int loc_J, int loc_K, \
    CUSTOMREAL* tau_loc) {

    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int k = blockIdx.z*blockDim.z+threadIdx.z;

    if (i >= loc_I || k >= loc_K || i < 0 || k < 0) return;

    int i_0_j_ = I2V_cuda(i, 0, k, loc_I, loc_J);
    int i_1_j_ = I2V_cuda(i, 1, k, loc_I, loc_J);
    int i_2_j_ = I2V_cuda(i, 2, k, loc_I, loc_J);
    int i_n1_j_ = I2V_cuda(i, loc_J-1, k, loc_I, loc_J);
    int i_n2_j_ = I2V_cuda(i, loc_J-2, k, loc_I, loc_J);
    int i_n3_j_ = I2V_cuda(i, loc_J-3, k, loc_I, loc_J);

    CUSTOMREAL v0 = 2.0 * tau_loc[i_1_j_] - tau_loc[i_2_j_];
    CUSTOMREAL v1 = tau_loc[i_2_j_];
    if (v0 > v1)
        tau_loc[i_0_j_] = v0;
    else
        tau_loc[i_0_j_] = v1;

    v0 = 2.0 * tau_loc[i_n2_j_] - tau_loc[i_n3_j_];
    v1 = tau_loc[i_n3_j_];
    if (v0 > v1)
        tau_loc[i_n1_j_] = v0;
    else
        tau_loc[i_n1_j_] = v1;
}


__global__ void cuda_calc_boundary_nodes_kernel_jk(\
    int loc_I, int loc_J, int loc_K, \
    CUSTOMREAL* tau_loc) {

    int j = blockIdx.y*blockDim.y+threadIdx.y;
    int k = blockIdx.z*blockDim.z+threadIdx.z;

    if (j >= loc_J || k >= loc_K || j < 0 || k < 0) return;

    int _0_j_k_ = I2V_cuda(0, j, k, loc_I, loc_J);
    int _1_j_k_ = I2V_cuda(1, j, k, loc_I, loc_J);
    int _2_j_k_ = I2V_cuda(2, j, k, loc_I, loc_J);
    int n1_j_k_ = I2V_cuda(loc_I-1, j, k, loc_I, loc_J);
    int n2_j_k_ = I2V_cuda(loc_I-2, j, k, loc_I, loc_J);
    int n3_j_k_ = I2V_cuda(loc_I-3, j, k, loc_I, loc_J);

    CUSTOMREAL v0 = 2.0 * tau_loc[_1_j_k_] \
                        - tau_loc[_2_j_k_];
    CUSTOMREAL v1 = tau_loc[_2_j_k_];
    if (v0 > v1)
        tau_loc[_0_j_k_] = v0;
    else
        tau_loc[_0_j_k_] = v1;

    v0 = 2.0 * tau_loc[n2_j_k_] \
             - tau_loc[n3_j_k_];
    v1 = tau_loc[n3_j_k_];
    if (v0 > v1)
        tau_loc[n1_j_k_] = v0;
    else
        tau_loc[n1_j_k_] = v1;
}


void cuda_calculate_boundary_nodes(Grid_on_device* grid) {
    dim3 blocks, threads;

    // i-j plane
    get_thread_block_for_kbound(grid->loc_I_host, grid->loc_J_host, 1, &threads, &blocks);
    cuda_calc_boundary_nodes_kernel_ij<<<blocks, threads>>>( \
                                    grid->loc_I_host, grid->loc_J_host, grid->loc_K_host, \
                                    grid->tau_loc_dev);

    // i-k plane
    get_thread_block_for_jbound(grid->loc_I_host, 1, grid->loc_K_host, &threads, &blocks);
    cuda_calc_boundary_nodes_kernel_ik<<<blocks, threads>>>( \
                                    grid->loc_I_host, grid->loc_J_host, grid->loc_K_host, \
                                    grid->tau_loc_dev);

    // j-k plane
    get_thread_block_for_ibound(1, grid->loc_J_host, grid->loc_K_host, &threads, &blocks);
    cuda_calc_boundary_nodes_kernel_jk<<<blocks, threads>>>( \
                                    grid->loc_I_host, grid->loc_J_host, grid->loc_K_host, \
                                    grid->tau_loc_dev);

}


void cuda_do_sweep_level(int iswp, Grid_on_device* grid_on_dv) {

    int actual_end_level = grid_on_dv->ed_level_host - grid_on_dv->st_level_host;
    //int start_point = 0; // count up for checking the start id of linearized level-point collection
    int r_dirc, t_dirc, p_dirc;

    // set sweep direction
    if (iswp == 0) {
        r_dirc = -1;
        t_dirc = -1;
        p_dirc = -1;
    } else if (iswp == 1) {
        r_dirc = -1;
        t_dirc = -1;
        p_dirc =  1;
    } else if (iswp == 2) {
        r_dirc = -1;
        t_dirc =  1;
        p_dirc = -1;
    } else if (iswp == 3) {
        r_dirc = -1;
        t_dirc =  1;
        p_dirc =  1;
    } else if (iswp == 4) {
        r_dirc =  1;
        t_dirc = -1;
        p_dirc = -1;
    } else if (iswp == 5) {
        r_dirc =  1;
        t_dirc = -1;
        p_dirc =  1;
    } else if (iswp == 6) {
        r_dirc =  1;
        t_dirc =  1;
        p_dirc = -1;
    } else if (iswp == 7) {
        r_dirc =  1;
        t_dirc =  1;
        p_dirc =  1;
    }

//
// slow method : definition of threads grid
//
//
// kernel launching by each level (slow by overheads for kernel launch)
//
//
// launch kernel for each level
// (this may be slower for an overhead to launch kernel each time,)
// (but at the moment use this, because grid synchronization is not working well)

    if (grid_on_dv->use_unified_kernel_host == false) {
        // define thread and block size for this level
        int block_size = CUDA_SWEEPING_BLOCK_SIZE;
        int num_blocks_x, num_blocks_y;

        get_block_xy(ceil(grid_on_dv->n_max_points_level_host/block_size+0.5), &num_blocks_x, &num_blocks_y);
        dim3 grid_each(num_blocks_x,num_blocks_y);
        dim3 threads_each(block_size,1,1);

        int id_start = 0;

        for (int i_level = 0; i_level <= actual_end_level; i_level++) {

            //get_block_xy(ceil(grid_on_dv->n_points_per_level_host[i_level]/block_size+0.5), &num_blocks_x, &num_blocks_y);
            //dim3 grid_each(num_blocks_x,num_blocks_y);
            //dim3 threads_each(block_size,1,1);

            void *kernelArgs[] = {\
                                    &i_level, \
                                    &id_start, \
                                    &actual_end_level, \
                                    &(grid_on_dv->n_points_per_level_dev), \
                                    &r_dirc, &t_dirc, &p_dirc, \
                                    &(grid_on_dv->i_id_all_level), \
                                    &(grid_on_dv->j_id_all_level), \
                                    &(grid_on_dv->k_id_all_level), \
                                    &grid_on_dv->loc_I_host, &grid_on_dv->loc_J_host, &grid_on_dv->loc_K_host, \
                                    &grid_on_dv->stencil_order_host, \
                                    &(grid_on_dv->fac_a_loc_dev), \
                                    &(grid_on_dv->fac_b_loc_dev), \
                                    &(grid_on_dv->fac_c_loc_dev), \
                                    &(grid_on_dv->fun_loc_dev), \
                                    &(grid_on_dv->T0v_loc_dev), \
                                    &(grid_on_dv->tau_loc_dev), \
                                    &grid_on_dv->dr_host, &grid_on_dv->dt_host, &grid_on_dv->dp_host, \
                                    &(grid_on_dv->fac_f_loc_dev), \
                                    &(grid_on_dv->T0r_loc_dev  ), \
                                    &(grid_on_dv->T0t_loc_dev  ), \
                                    &(grid_on_dv->T0p_loc_dev  ), \
                                    &(grid_on_dv->is_changed_dev) \
                                    };

            int id_stream = i_level % CUDA_MAX_NUM_STREAMS;
            if (grid_on_dv->stencil_order_host == 3)
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_each, threads_each, kernelArgs, 0, grid_on_dv->level_streams[id_stream]), 30001);
            else if (grid_on_dv->stencil_order_host == 1)
                print_CUDA_error_if_any(cudaLaunchKernel((void*) cuda_do_sweep_level_kernel_1st, grid_each, threads_each, kernelArgs, 0, grid_on_dv->level_streams[id_stream]), 30002);

            id_start += grid_on_dv->n_points_per_level_host[i_level];

            // sychronize the stream for this level
            //cudaStreamSynchronize(grid_on_dv->level_streams[i_level]);

        }

    } else {

        //
        // launch kernel only once for all levels
        //
        int i_level_dummy = 9999;
        int id_start_dummy = 9999;

        void *kernelArgs[] = {\
                                &i_level_dummy, \
                                &id_start_dummy, \
                                &actual_end_level, \
                                &(grid_on_dv->n_points_per_level_dev), \
                                &r_dirc, &t_dirc, &p_dirc, \
                                &(grid_on_dv->i_id_all_level), \
                                &(grid_on_dv->j_id_all_level), \
                                &(grid_on_dv->k_id_all_level), \
                                &grid_on_dv->loc_I_host, &grid_on_dv->loc_J_host, &grid_on_dv->loc_K_host, \
                                &grid_on_dv->stencil_order_host, \
                                &(grid_on_dv->fac_a_loc_dev), \
                                &(grid_on_dv->fac_b_loc_dev), \
                                &(grid_on_dv->fac_c_loc_dev), \
                                &(grid_on_dv->fun_loc_dev), \
                                &(grid_on_dv->T0v_loc_dev), \
                                &(grid_on_dv->tau_loc_dev), \
                                &grid_on_dv->dr_host, &grid_on_dv->dt_host, &grid_on_dv->dp_host, \
                                &(grid_on_dv->fac_f_loc_dev), \
                                &(grid_on_dv->T0r_loc_dev  ), \
                                &(grid_on_dv->T0t_loc_dev  ), \
                                &(grid_on_dv->T0p_loc_dev  ), \
                                &(grid_on_dv->is_changed_dev) \
                                };
        if (grid_on_dv->stencil_order_host == 3)
            print_CUDA_error_if_any(cudaLaunchCooperativeKernel((void*) cuda_do_sweep_level_kernel_3rd, grid_on_dv->grid_sweep_host, grid_on_dv->threads_sweep_host, kernelArgs), 30001);
        else if (grid_on_dv->stencil_order_host == 1)
            print_CUDA_error_if_any(cudaLaunchCooperativeKernel((void*) cuda_do_sweep_level_kernel_1st, grid_on_dv->grid_sweep_host, grid_on_dv->threads_sweep_host, kernelArgs), 30002);

    }


    // update boundary
    cuda_calculate_boundary_nodes(grid_on_dv);

}


__global__ void cuda_calc_T_plus_tau_kernel(\
    int loc_I, int loc_J, int loc_K, \
    CUSTOMREAL* T0v_loc, \
    CUSTOMREAL* tau_loc, \
    CUSTOMREAL* T_loc) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    // skip if the assigned ijk is greater than the grid size
    if (i >= loc_I || j >= loc_J || k >= loc_K)
        return;

    int ijk = I2V_cuda(i,j,k, loc_I, loc_J);
    T_loc[ijk] = T0v_loc[ijk] * tau_loc[ijk];

}


void cuda_calc_T_plus_tau(Grid_on_device* grid_on_dv) {
    dim3 blocks, threads;
    get_thread_block_for_3d_loop(grid_on_dv->loc_I_host, grid_on_dv->loc_J_host, grid_on_dv->loc_K_host, &threads, &blocks);

    cuda_calc_T_plus_tau_kernel<<<blocks, threads>>>(grid_on_dv->loc_I_host, grid_on_dv->loc_J_host, grid_on_dv->loc_K_host, \
                                                     grid_on_dv->T0v_loc_dev, \
                                                     grid_on_dv->tau_loc_dev, \
                                                     grid_on_dv->T_loc_dev);
}


void cuda_copy_T_loc_tau_loc_to_host(CUSTOMREAL* h_T_loc, CUSTOMREAL* h_tau_loc, Grid_on_device* grid) {
    print_CUDA_error_if_any(copy_device_to_host_cv(h_T_loc,   grid->T_loc_dev,   grid->loc_I_host*grid->loc_J_host*grid->loc_K_host),20001);
    print_CUDA_error_if_any(copy_device_to_host_cv(h_tau_loc, grid->tau_loc_dev, grid->loc_I_host*grid->loc_J_host*grid->loc_K_host),20002);
}


void initialize_sweep_params(Grid_on_device* grid_on_dv) {

    // set threads and grid for 3d loop
    get_thread_block_for_3d_loop(grid_on_dv->loc_I_host, grid_on_dv->loc_J_host, grid_on_dv->loc_K_host, &(grid_on_dv->threads_3d_full_incr_host), &(grid_on_dv->grid_3d_full_incr_host));

    // for L1 calculation
    int n_blocks_x, n_blocks_y;
    int n_theads_per_block = CUDA_L1_BLOCK_SIZE;
    int n_total_pt = grid_on_dv->loc_I_host*grid_on_dv->loc_J_host*grid_on_dv->loc_K_host;
    get_block_xy(ceil(n_total_pt/n_theads_per_block+0.5), &n_blocks_x, &n_blocks_y);
    grid_on_dv->grid_L1_host     = dim3(n_blocks_x, n_blocks_y, 1);
    grid_on_dv->threads_L1_host  = dim3(n_theads_per_block, 1, 1);
    grid_on_dv->n_blocks_L1_host = grid_on_dv->grid_L1_host.x * grid_on_dv->grid_L1_host.y;

    // check the numBlockPerSm and set the block size accordingly
    int numBlocksPerSm = 0;
    int block_size = CUDA_SWEEPING_BLOCK_SIZE;

    int device;
    cudaGetDevice(&device);

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, device);
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm, cuda_do_sweep_level_kernel_3rd, CUDA_SWEEPING_BLOCK_SIZE, 0);

//    printf("numBlocksPerSm = %d\n", numBlocksPerSm);
//    printf("multiprocessor_count = %d\n", deviceProp.multiProcessorCount);
//    printf("max number of blocks for cooperative kernel launch: %d\n", deviceProp.multiProcessorCount*numBlocksPerSm);
//    printf("max number of points on each level: %d\n", deviceProp.multiProcessorCount*numBlocksPerSm*deviceProp.maxThreadsPerBlock);
//    printf("current number of points on each level: %d\n", grid_on_dv->n_max_points_level_host);

    int max_cooperative_blocks = deviceProp.multiProcessorCount*numBlocksPerSm;

    grid_on_dv->threads_sweep_host = dim3(block_size, 1, 1);
    grid_on_dv->grid_sweep_host = dim3(max_cooperative_blocks, 1, 1);


    // don't use unified kernel if the device has insufficient memory
    // abort if max_cooperative_blocks*threads is smaller than the max number of points on each level
    if (max_cooperative_blocks*block_size < grid_on_dv->n_max_points_level_host || FORCE_UNUSE_UNIFIED_KERNEL) {
        //#TODO: for the case that the block size for cooperative launch is not enogh,
        // we need to divide the large level int to smaller levels.

//        printf("max_cooperative_blocks*block_size < grid_on_dv->n_max_points_level_host\n");
//        printf("max_cooperative_blocks*block_size = %d\n", max_cooperative_blocks*block_size);
//        printf("grid_on_dv->n_max_points_level_host = %d\n", grid_on_dv->n_max_points_level_host);
        //printf("abort\n");
        //exit(1);
        // allocate level_streams
        // required only for kernel launching for each level
        //
        //int n_levels = grid_on_dv->ed_level_host - grid_on_dv->st_level_host + 1;
        grid_on_dv->level_streams = (cudaStream_t*)malloc(CUDA_MAX_NUM_STREAMS*sizeof(cudaStream_t));
        // spawn streams
        for (int i = 0; i < CUDA_MAX_NUM_STREAMS; i++) {
            cudaStreamCreate(&(grid_on_dv->level_streams[i]));
        }
        grid_on_dv->use_unified_kernel_host = false;
    } else {
        printf("use unified kernel for all levels\n");
        grid_on_dv->use_unified_kernel_host = true;
    }

}


void finalize_sweep_params(Grid_on_device* grid_on_dv){
    if (grid_on_dv->use_unified_kernel_host == false)
        free(grid_on_dv->level_streams);
}


void cuda_run_iterate_forward(Grid_on_device* grid_on_dv, \
                              CUSTOMREAL tolerance, \
                              int max_iter, \
                              int stencil_order_host) {
    //
    //printf("cuda_run_iterate_forward\n");

    grid_on_dv->stencil_order_host = stencil_order_host;

    initialize_sweep_params(grid_on_dv);

    const int nswp = 8;
    int iter_count = 0;

    // calculate initial L1
    cuda_calculate_L1(grid_on_dv);

    // start iteration
    while (true) {
        // store tau to tau_old
        cuda_tau2old_tau(grid_on_dv);

        // do sweeping for all direction
        for (int iswp = nswp-1; iswp > -1; iswp--) {
            // do sweeping
            cuda_do_sweep_level(iswp, grid_on_dv);
        }

        // copy the values of communication nodes to the ghost nodes of adjacent subdomains
        cuda_send_recev_boundary_data(grid_on_dv);

        // calculate the obj for this subdomain and mpi reduce
        cuda_calculate_L1(grid_on_dv);

        // check convergence
        if (grid_on_dv->L1_host < tolerance) {
            goto iter_end;
        } else if (iter_count >= max_iter) {
            goto iter_end;
        } else {
            if (grid_on_dv->myrank_host== 0)
                printf("iter_count, L1 , tolerance : %d, %5.16f, %5.16f\n", iter_count, grid_on_dv->L1_host, tolerance);
            iter_count++;
        }
    }

iter_end:
    // after convergence
    printf("GPU iteration end, iter_count, L1 = %d, %5.16f\n", iter_count, grid_on_dv->L1_host);

    // calculate T
    cuda_calc_T_plus_tau(grid_on_dv);

    finalize_sweep_params(grid_on_dv);

}