#ifndef VECTORIZED_SWEEP_H
#define VECTORIZED_SWEEP_H


#ifdef USE_AVX
    #include <immintrin.h>



#include<vector>

#define NSIMD 4
#define COEF 1.0
//#define CUSTOMREAL double

//int NP, NT, NR;
//#define I2V(i,j,k) ((k)*NP*NT+(j)*NP+(i))



inline void fake_stencil_3rd_pre_simd(__m256i& v_iip, __m256i& v_jjt, __m256i& v_kkr, __m256d& v_c__,  __m256d& v_p__,     __m256d& v_m__,     __m256d& v__p_,    __m256d& v__m_,    __m256d& v___p,    __m256d& v___m, \
                                                                                                       __m256d& v_pp____,  __m256d& v_mm____,  __m256d& v___pp__, __m256d& v___mm__, __m256d& v_____pp, __m256d& v_____mm, \
                                                                                                       __m256d& v_pp1,     __m256d& v_pp2,     __m256d& v_pt1,    __m256d& v_pt2,    __m256d& v_pr1,    __m256d& v_pr2, \
                                                                                                       CUSTOMREAL& DP, CUSTOMREAL& DT, CUSTOMREAL& DR, \
                                                                                                       int& NP, int& NT, int& NR){

    __m256d v_wp1, v_wp2, v_wt1, v_wt2, v_wr1, v_wr2;
    __m256d v_eps = _mm256_set1_pd(1e-6);
    __m256d v_DR_inv = _mm256_set1_pd(1.0/DR);
    __m256d v_DT_inv = _mm256_set1_pd(1.0/DT);
    __m256d v_DP_inv = _mm256_set1_pd(1.0/DP);

    // initialize
    v_wp1 = _mm256_setzero_pd();
    v_wp2 = _mm256_setzero_pd();
    v_wt1 = _mm256_setzero_pd();
    v_wt2 = _mm256_setzero_pd();
    v_wr1 = _mm256_setzero_pd();
    v_wr2 = _mm256_setzero_pd();

    // debug
    // set one
    v_pp1 = _mm256_set1_pd(1.0);
    v_pp2 = _mm256_set1_pd(1.0);
    v_pt1 = _mm256_set1_pd(1.0);
    v_pt2 = _mm256_set1_pd(1.0);
    v_pr1 = _mm256_set1_pd(1.0);
    v_pr2 = _mm256_set1_pd(1.0);

    // check if the stencil is in the domain or partially out of the domain
    __m256i mask_i_eq_1         = _mm256_cmpeq_epi32(v_iip, _mm256_set1_epi32(1));           // if iip == 1
    __m256i mask_i_eq_N_minus_1 = _mm256_cmpeq_epi32(v_iip, _mm256_set1_epi32(NP-1)); // if iip == N-1
    __m256i mask_i_else         = _mm256_andnot_si256(_mm256_or_si256(mask_i_eq_1, mask_i_eq_N_minus_1), _mm256_set1_epi32(-1));
    __m256i mask_j_eq_1         = _mm256_cmpeq_epi32(v_jjt, _mm256_set1_epi32(1));           // if jjt == 1
    __m256i mask_j_eq_N_minus_1 = _mm256_cmpeq_epi32(v_jjt, _mm256_set1_epi32(NT-1)); // if jjt == N-1
    __m256i mask_j_else         = _mm256_andnot_si256(_mm256_or_si256(mask_j_eq_1, mask_j_eq_N_minus_1), _mm256_set1_epi32(-1));
    __m256i mask_k_eq_1         = _mm256_cmpeq_epi32(v_kkr, _mm256_set1_epi32(1));           // if kkr == 1
    __m256i mask_k_eq_N_minus_1 = _mm256_cmpeq_epi32(v_kkr, _mm256_set1_epi32(NR-1)); // if kkr == N-1
    __m256i mask_k_else         = _mm256_andnot_si256(_mm256_or_si256(mask_k_eq_1, mask_k_eq_N_minus_1), _mm256_set1_epi32(-1));

    // if _i_eq_1 == true
    v_pp1 = _mm256_blendv_pd(v_pp1,_mm256_mul_pd(_mm256_sub_pd(v_c__,v_m__),v_DP_inv), _mm256_castsi256_pd(mask_i_eq_1));
    // v_eps + sqrt(v_c__ - 2.0*v_p__ + v_pp____)
    __m256d tmp1 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v_c__,_mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v_p__),v_pp____))));
    // v_eps + sqrt(v_m__ - 2.0*v_c__ + v_p__)
    __m256d tmp2 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v_m__,_mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v_c__),v_p__))));
    // v_wp2 = 1.0/(1.0 + 2.0 * sqrt(tmp1/tmp2))
    v_wp2 = _mm256_blendv_pd(v_wp2,_mm256_div_pd(_mm256_set1_pd(1.0),_mm256_add_pd(_mm256_set1_pd(1.0),_mm256_mul_pd(_mm256_set1_pd(2.0),_mm256_sqrt_pd(_mm256_div_pd(tmp1,tmp2))))), _mm256_castsi256_pd(mask_i_eq_1));
    // v_pp2 = (1.0 - v_wp2) * (v_p__ - v_m__) / 2.0 / D
    // + v_wp2 * (-3.0* v_c__ + 4.0 * v_p__ - v_pp____) / 2.0 / D
    v_pp2 = _mm256_blendv_pd(v_pp2,_mm256_add_pd(\
                          _mm256_mul_pd(_mm256_sub_pd(_mm256_set1_pd(1.0),v_wp2),_mm256_mul_pd(_mm256_sub_pd(v_p__,v_m__),_mm256_set1_pd(0.5*DP))),\
                          _mm256_mul_pd(v_wp2,_mm256_mul_pd(_mm256_sub_pd(_mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(4.0),v_p__),_mm256_mul_pd(_mm256_set1_pd(3.0),v_c__)),v_pp____),_mm256_set1_pd(0.5*DP)))), _mm256_castsi256_pd(mask_i_eq_1));

    // if_i_eq_N_minus_1 == true
    // v_eps + sqrt(v_mm____ - 2.0*v_m__ + v_c__)
    tmp1 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v_mm____,_mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v_m__),v_c__))));
    // v_eps + sqrt(v_m__ - 2.0*v_c__ + v_p__)
    tmp2 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v_m__,   _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v_c__),v_p__))));
    // v_wp1 = 1.0/(1.0 + 2.0 * sqrt(tmp1/tmp2))
    v_wp1 = _mm256_blendv_pd(v_wp1,_mm256_div_pd(_mm256_set1_pd(1.0),_mm256_add_pd(_mm256_set1_pd(1.0),_mm256_mul_pd(_mm256_set1_pd(2.0),_mm256_sqrt_pd(_mm256_div_pd(tmp1,tmp2))))), _mm256_castsi256_pd(mask_i_eq_N_minus_1));
    // v_pp1 = (1.0 - v_wp1) * (v_p__ - v_m__) / 2.0 / D
    // + v_wp1 * (v_mm____ - 4.0 * v_m__ + 3.0 * v_c__) / 2.0 / D
    v_pp1 = _mm256_blendv_pd(v_pp1,_mm256_add_pd(\
                          _mm256_mul_pd(_mm256_sub_pd(_mm256_set1_pd(1.0),v_wp1),_mm256_mul_pd(_mm256_sub_pd(v_p__,v_m__),_mm256_set1_pd(0.5*DP))),\
                          _mm256_mul_pd(v_wp1,_mm256_mul_pd(_mm256_sub_pd(v_mm____,_mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(4.0),v_m__),_mm256_mul_pd(_mm256_set1_pd(3.0),v_c__))),_mm256_set1_pd(0.5*DP)))), _mm256_castsi256_pd(mask_i_eq_N_minus_1));
    // v_pp2 = (v_p__ - v_c__) * v_D_inv
    v_pp2 = _mm256_blendv_pd(v_pp2,_mm256_mul_pd(_mm256_sub_pd(v_p__,v_c__),v_DP_inv), _mm256_castsi256_pd(mask_i_eq_N_minus_1));

    // else
    // v_eps + sqrt(v_mm____ - 2.0*v_m__ + v_c__)
    tmp1 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v_mm____,_mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v_m__),v_c__))));
    // v_eps + sqrt(v_m__ - 2.0*v_c__ + v_p__)
    tmp2 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v_m__,   _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v_c__),v_p__))));
    // v_wp1 = 1.0/(1.0 + 2.0 * sqrt(tmp1/tmp2))
    v_wp1 = _mm256_blendv_pd(v_wp1,_mm256_div_pd(_mm256_set1_pd(1.0),_mm256_add_pd(_mm256_set1_pd(1.0),_mm256_mul_pd(_mm256_set1_pd(2.0),_mm256_sqrt_pd(_mm256_div_pd(tmp1,tmp2))))), _mm256_castsi256_pd(mask_i_else));
    // v_pp1 = (1.0 - v_wp1) * (v_p__ - v_m__) / 2.0 / D
    // + v_wp1 * (v_mm____ - 4.0 * v_m__ + 3.0 * v_c__) / 2.0 / D
    v_pp1 = _mm256_blendv_pd(v_pp1,_mm256_add_pd(\
                          _mm256_mul_pd(_mm256_sub_pd(_mm256_set1_pd(1.0),v_wp1),_mm256_mul_pd(_mm256_sub_pd(v_p__,v_m__),_mm256_set1_pd(0.5*DP))),\
                          _mm256_mul_pd(v_wp1,_mm256_mul_pd(_mm256_sub_pd(v_mm____,_mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(4.0),v_m__),_mm256_mul_pd(_mm256_set1_pd(3.0),v_c__))),_mm256_set1_pd(0.5*DP)))), _mm256_castsi256_pd(mask_i_else));
    // v_esp + sqrt(v_c__ - 2.0*v_p__ + v_pp____)
    tmp1 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v_c__,   _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v_p__),v_pp____))));
    // v_eps + sqrt(v_m__ - 2.0*v_c__ + v_p__)
    tmp2 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v_m__,   _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v_c__),v_p__))));
    // v_wp2 = 1.0/(1.0 + 2.0 * sqrt(tmp1/tmp2))
    v_wp2 = _mm256_blendv_pd(v_wp2,_mm256_div_pd(_mm256_set1_pd(1.0),_mm256_add_pd(_mm256_set1_pd(1.0),_mm256_mul_pd(_mm256_set1_pd(2.0),_mm256_sqrt_pd(_mm256_div_pd(tmp1,tmp2))))), _mm256_castsi256_pd(mask_i_else));
    // v_pp2 = (1.0 - v_wp2) * (v_p__ - v_m__) / 2.0 / D
    // + v_wp2 * (-3.0 * v_c__ + 4.0 * v_p__ - v_pp____) / 2.0 / D
    v_pp2 = _mm256_blendv_pd(v_pp2,_mm256_add_pd(\
                          _mm256_mul_pd(_mm256_sub_pd(_mm256_set1_pd(1.0),v_wp2),_mm256_mul_pd(_mm256_sub_pd(v_p__,v_m__),_mm256_set1_pd(0.5*DP))),\
                          _mm256_mul_pd(v_wp2,_mm256_mul_pd(_mm256_sub_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(-3.0),v_c__),_mm256_mul_pd(_mm256_set1_pd(4.0),v_p__)),v_pp____),_mm256_set1_pd(0.5*DP)))), _mm256_castsi256_pd(mask_i_else));

    // if _j_eq_1 == true
    // v_pt1 = (v_c__ - v__m_) * v_D_inv
    v_pt1 = _mm256_blendv_pd(v_pt1,_mm256_mul_pd(_mm256_sub_pd(v_c__,v__m_),v_DT_inv), _mm256_castsi256_pd(mask_j_eq_1));
    // v_esp + sqrt(v_c__ - 2.0*v__p_ + v___pp__)
    tmp1 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v_c__,   _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v__p_),v___pp__))));
    // v_eps + sqrt(v__m_ - 2.0*v_c__ + v__p_)
    tmp2 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v__m_,   _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v_c__),v__p_))));
    // v_wt2 = 1.0/(1.0 + 2.0 * sqrt(tmp1/tmp2))
    v_wt2 = _mm256_blendv_pd(v_wp1,_mm256_div_pd(_mm256_set1_pd(1.0),_mm256_add_pd(_mm256_set1_pd(1.0),_mm256_mul_pd(_mm256_set1_pd(2.0),_mm256_sqrt_pd(_mm256_div_pd(tmp1,tmp2))))), _mm256_castsi256_pd(mask_j_eq_1));
    // v_pt2 = (1.0 - v_wt2) * (v__p_ - v__m_) / 2.0 / D
    // + v_wt2 * (-3.0 * v_c__ + 4.0 * v__p_ - v___mm__) / 2.0 / D
    v_pt2 = _mm256_blendv_pd(v_pt2,_mm256_add_pd(\
                          _mm256_mul_pd(_mm256_sub_pd(_mm256_set1_pd(1.0),v_wt2),_mm256_mul_pd(_mm256_sub_pd(v__p_,v__m_),_mm256_set1_pd(0.5*DT))),\
                          _mm256_mul_pd(v_wt2,_mm256_mul_pd(_mm256_sub_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(-3.0),v_c__),_mm256_mul_pd(_mm256_set1_pd(4.0),v__p_)),v___mm__),_mm256_set1_pd(0.5*DT)))), _mm256_castsi256_pd(mask_j_eq_1));
    // if _j_eq_N_minus_1 == true
    // v_esp + sqrt(v___mm__ - 2.0*v__m_ + v_c__)
    tmp1 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v___mm__,_mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v__m_),v_c__))));
    // v_eps + sqrt(v__m_ - 2.0*v_c__ + v__m_)
    tmp2 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v__m_,   _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v_c__),v__m_))));
    // v_wt1 = 1.0/(1.0 + 2.0 * sqrt(tmp1/tmp2))
    v_wt1 = _mm256_blendv_pd(v_wt2,_mm256_div_pd(_mm256_set1_pd(1.0),_mm256_add_pd(_mm256_set1_pd(1.0),_mm256_mul_pd(_mm256_set1_pd(2.0),_mm256_sqrt_pd(_mm256_div_pd(tmp1,tmp2))))), _mm256_castsi256_pd(mask_j_eq_N_minus_1));
    // v_pt1 = (1.0 - v_wt1) * (v__p_ - v__m_) / 2.0 / D
    // + v_wt1 * (v___mm__ - 4.0 * v__m_ + 3.0 * v_c__) / 2.0 / D
    v_pt1 = _mm256_blendv_pd(v_pt1,_mm256_add_pd(\
                          _mm256_mul_pd(_mm256_sub_pd(_mm256_set1_pd(1.0),v_wt1),_mm256_mul_pd(_mm256_sub_pd(v__p_,v__m_),_mm256_set1_pd(0.5*DT))),\
                          _mm256_mul_pd(v_wt1,_mm256_mul_pd(_mm256_add_pd(_mm256_sub_pd(v___mm__,_mm256_mul_pd(_mm256_set1_pd(4.0),v__m_)),_mm256_mul_pd(_mm256_set1_pd(3.0),v_c__)),_mm256_set1_pd(0.5*DT)))), _mm256_castsi256_pd(mask_j_eq_N_minus_1));
    // v_pt2 = (v__p_ - v_c__) * v_D_inv
    v_pt2 = _mm256_blendv_pd(v_pt2,_mm256_mul_pd(_mm256_sub_pd(v__p_,v_c__),v_DT_inv), _mm256_castsi256_pd(mask_j_eq_N_minus_1));
    // else
    // v_esp + sqrt(v___mm__ - 2.0*v__m_ + v_c__)
    tmp1 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v___mm__,_mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v__m_),v_c__))));
    // v_eps + sqrt(v__m_ - 2.0*v_c__ + v__m_)
    tmp2 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v__m_,   _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v_c__),v__m_))));
    // v_wt1 = 1.0/(1.0 + 2.0 * sqrt(tmp1/tmp2))
    v_wt1 = _mm256_blendv_pd(v_wt2,_mm256_div_pd(_mm256_set1_pd(1.0),_mm256_add_pd(_mm256_set1_pd(1.0),_mm256_mul_pd(_mm256_set1_pd(2.0),_mm256_sqrt_pd(_mm256_div_pd(tmp1,tmp2))))), _mm256_castsi256_pd(mask_j_else));
    // v_pt1 = (1.0 - v_wt1) * (v__p_ - v__m_) / 2.0 / D
    // + v_wt1 * (v___mm__ - 4.0 * v__m_ + 3.0 * v_c__) / 2.0 / D
    v_pt1 = _mm256_blendv_pd(v_pt1,_mm256_add_pd(\
                          _mm256_mul_pd(_mm256_sub_pd(_mm256_set1_pd(1.0),v_wt1),_mm256_mul_pd(_mm256_sub_pd(v__p_,v__m_),_mm256_set1_pd(0.5*DT))),\
                          _mm256_mul_pd(v_wt1,_mm256_mul_pd(_mm256_add_pd(_mm256_sub_pd(v___mm__,_mm256_mul_pd(_mm256_set1_pd(4.0),v__m_)),_mm256_mul_pd(_mm256_set1_pd(3.0),v_c__)),_mm256_set1_pd(0.5*DT)))), _mm256_castsi256_pd(mask_j_else));
    tmp1 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v_c__,   _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v__p_),v___pp__))));
    // v_eps + sqrt(v__m_ - 2.0*v_c__ + v__p_)
    tmp2 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v__m_,   _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v_c__),v__p_))));
    // v_wt2 = 1.0/(1.0 + 2.0 * sqrt(tmp1/tmp2))
    v_wt2 = _mm256_blendv_pd(v_wp1,_mm256_div_pd(_mm256_set1_pd(1.0),_mm256_add_pd(_mm256_set1_pd(1.0),_mm256_mul_pd(_mm256_set1_pd(2.0),_mm256_sqrt_pd(_mm256_div_pd(tmp1,tmp2))))), _mm256_castsi256_pd(mask_j_else));
    // v_pt2 = (1.0 - v_wt2) * (v__p_ - v__m_) / 2.0 / D
    // + v_wt2 * (-3.0 * v_c__ + 4.0 * v__p_ - v___mm__) / 2.0 / D
    v_pt2 = _mm256_blendv_pd(v_pt2,_mm256_add_pd(\
                          _mm256_mul_pd(_mm256_sub_pd(_mm256_set1_pd(1.0),v_wt2),_mm256_mul_pd(_mm256_sub_pd(v__p_,v__m_),_mm256_set1_pd(0.5*DT))),\
                          _mm256_mul_pd(v_wt2,_mm256_mul_pd(_mm256_sub_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(-3.0),v_c__),_mm256_mul_pd(_mm256_set1_pd(4.0),v__p_)),v___mm__),_mm256_set1_pd(0.5*DT)))), _mm256_castsi256_pd(mask_j_else));

    // if _k_eq_1 == true
    // v_pr1 = (v_c__ - v___m) * v_D_inv
    v_pr1 = _mm256_blendv_pd(v_pr1,_mm256_mul_pd(_mm256_sub_pd(v_c__,v___m),v_DR_inv), _mm256_castsi256_pd(mask_k_eq_1));
    // v_eps + sqrt(v_c__ - 2.0*v___p + v_____pp)
    tmp1 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v_c__,   _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v___p),v_____pp))));
    // v_eps + sqrt(v___m - 2.0*v_c__ + v___p)
    tmp2 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v___m,   _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v_c__),v___p))));
    // v_wr2 = 1.0/(1.0 + 2.0 * sqrt(tmp1/tmp2))
    v_wr2 = _mm256_blendv_pd(v_wr2,_mm256_div_pd(_mm256_set1_pd(1.0),_mm256_add_pd(_mm256_set1_pd(1.0),_mm256_mul_pd(_mm256_set1_pd(2.0),_mm256_sqrt_pd(_mm256_div_pd(tmp1,tmp2))))), _mm256_castsi256_pd(mask_k_eq_1));
    // v_pr2 = (1.0 - v_wr2) * (v___p - v___m) / 2.0 / D
    // + v_wr2 * (-3.0 * v_c__ + 4.0 * v___p - v_____pp) / 2.0 / D
    v_pr2 = _mm256_blendv_pd(v_pr2,_mm256_add_pd(\
                          _mm256_mul_pd(_mm256_sub_pd(_mm256_set1_pd(1.0),v_wr2),_mm256_mul_pd(_mm256_sub_pd(v___p,v___m),_mm256_set1_pd(0.5*DR))),\
                          _mm256_mul_pd(v_wr2,_mm256_mul_pd(_mm256_sub_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(-3.0),v_c__),_mm256_mul_pd(_mm256_set1_pd(4.0),v___p)),v_____pp),_mm256_set1_pd(0.5*DR)))), _mm256_castsi256_pd(mask_k_eq_1));
    // if _k_eq_N_minus_1 == true
    // eps + sqrt(v_c__ - 2.0*v___m + v_____mm)
    tmp1 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v_c__,   _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v___m),v_____mm))));
    // eps + sqrt(v___m - 2.0*v_c__ + v___p)
    tmp2 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v___m,   _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v_c__),v___p))));
    // v_wr1 = 1.0/(1.0 + 2.0 * sqrt(tmp1/tmp2))
    v_wr1 = _mm256_blendv_pd(v_wr1,_mm256_div_pd(_mm256_set1_pd(1.0),_mm256_add_pd(_mm256_set1_pd(1.0),_mm256_mul_pd(_mm256_set1_pd(2.0),_mm256_sqrt_pd(_mm256_div_pd(tmp1,tmp2))))), _mm256_castsi256_pd(mask_k_eq_N_minus_1));
    // v_pr1 = (1.0 - v_wr1) * (v___p - v___m) / 2.0 / D
    // + v_wr1 * (v_____mm - 4.0 * v___m + 3.0 * v_c__) / 2.0 / D
    v_pr1 = _mm256_blendv_pd(v_pr1,_mm256_add_pd(\
                          _mm256_mul_pd(_mm256_sub_pd(_mm256_set1_pd(1.0),v_wr1),_mm256_mul_pd(_mm256_sub_pd(v___p,v___m),_mm256_set1_pd(0.5*DR))),\
                          _mm256_mul_pd(v_wr1,_mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(v_____mm,_mm256_mul_pd(_mm256_set1_pd(-4.0),v___m)),_mm256_mul_pd(_mm256_set1_pd(3.0),v_c__)),_mm256_set1_pd(0.5*DR)))), _mm256_castsi256_pd(mask_k_eq_N_minus_1));
    // v_pr2 = (v___p - v_c__) * v_D_inv
    v_pr2 = _mm256_blendv_pd(v_pr2,_mm256_mul_pd(_mm256_sub_pd(v___p,v_c__),v_DR_inv), _mm256_castsi256_pd(mask_k_eq_N_minus_1));
    // else
    // eps + sqrt(v_c__ - 2.0*v___m + v_____mm)
    tmp1 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v_c__,   _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v___m),v_____mm))));
    // eps + sqrt(v___m - 2.0*v_c__ + v___p)
    tmp2 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v___m,   _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v_c__),v___p))));
    // v_wr1 = 1.0/(1.0 + 2.0 * sqrt(tmp1/tmp2))
    v_wr1 = _mm256_blendv_pd(v_wr1,_mm256_div_pd(_mm256_set1_pd(1.0),_mm256_add_pd(_mm256_set1_pd(1.0),_mm256_mul_pd(_mm256_set1_pd(2.0),_mm256_sqrt_pd(_mm256_div_pd(tmp1,tmp2))))), _mm256_castsi256_pd(mask_k_else));
    // v_pr1 = (1.0 - v_wr1) * (v___p - v___m) / 2.0 / D
    // + v_wr1 * (v_____mm - 4.0 * v___m + 3.0 * v_c__) / 2.0 / D
    v_pr1 = _mm256_blendv_pd(v_pr1,_mm256_add_pd(\
                          _mm256_mul_pd(_mm256_sub_pd(_mm256_set1_pd(1.0),v_wr1),_mm256_mul_pd(_mm256_sub_pd(v___p,v___m),_mm256_set1_pd(0.5*DR))),\
                          _mm256_mul_pd(v_wr1,_mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(v_____mm,_mm256_mul_pd(_mm256_set1_pd(-4.0),v___m)),_mm256_mul_pd(_mm256_set1_pd(3.0),v_c__)),_mm256_set1_pd(0.5*DR)))), _mm256_castsi256_pd(mask_k_else));
    tmp1 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v_c__,   _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v___p),v_____pp))));
    // v_eps + sqrt(v___m - 2.0*v_c__ + v___p)
    tmp2 = _mm256_add_pd(v_eps,_mm256_sqrt_pd(_mm256_sub_pd(v___m,   _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0),v_c__),v___p))));
    // v_wr2 = 1.0/(1.0 + 2.0 * sqrt(tmp1/tmp2))
    v_wr2 = _mm256_blendv_pd(v_wr2,_mm256_div_pd(_mm256_set1_pd(1.0),_mm256_add_pd(_mm256_set1_pd(1.0),_mm256_mul_pd(_mm256_set1_pd(2.0),_mm256_sqrt_pd(_mm256_div_pd(tmp1,tmp2))))), _mm256_castsi256_pd(mask_k_else));
    // v_pr2 = (1.0 - v_wr2) * (v___p - v___m) / 2.0 / D
    // + v_wr2 * (-3.0 * v_c__ + 4.0 * v___p - v_____pp) / 2.0 / D
    v_pr2 = _mm256_blendv_pd(v_pr2,_mm256_add_pd(\
                          _mm256_mul_pd(_mm256_sub_pd(_mm256_set1_pd(1.0),v_wr2),_mm256_mul_pd(_mm256_sub_pd(v___p,v___m),_mm256_set1_pd(0.5*DR))),\
                          _mm256_mul_pd(v_wr2,_mm256_mul_pd(_mm256_sub_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(-3.0),v_c__),_mm256_mul_pd(_mm256_set1_pd(4.0),v___p)),v_____pp),_mm256_set1_pd(0.5*DR)))), _mm256_castsi256_pd(mask_k_else));

}


// tau fac_a fac_b fac_c fac_f T0v T0p T0t T0r fun
inline void fake_stencil_3rd_apre_simd(__m256d& v_tau,__m256d& v_fac_a,__m256d& v_fac_b, __m256d& v_fac_c, __m256d& v_fac_f, \
                                       __m256d& v_T0v, __m256d& v_T0p, __m256d& v_T0t, __m256d& v_T0r, __m256d& v_fun, \
                                       __m256d& v_pp1, __m256d& v_pp2, __m256d& v_pt1, __m256d& v_pt2, __m256d& v_pr1, __m256d& v_pr2, \
                                       CUSTOMREAL& DP, CUSTOMREAL& DT, CUSTOMREAL& DR){

    // sigr = COEF * sqrt(v_fac_a)*v_T0v;
    __m256d sigr = _mm256_mul_pd(_mm256_mul_pd(_mm256_set1_pd(COEF),_mm256_sqrt_pd(v_fac_a)),v_T0v);
    // sigt = COEF * sqrt(v_fac_b)*v_T0v;
    __m256d sigt = _mm256_mul_pd(_mm256_mul_pd(_mm256_set1_pd(COEF),_mm256_sqrt_pd(v_fac_b)),v_T0v);
    // sigp = COEF * sqrt(v_fac_c)*v_T0v;
    __m256d sigp = _mm256_mul_pd(_mm256_mul_pd(_mm256_set1_pd(COEF),_mm256_sqrt_pd(v_fac_c)),v_T0v);

    // coe = 1.0 / (sigr/D + sigt/D + sigz/D);
    __m256d DP_inv = _mm256_set1_pd(1.0/DP);
    __m256d DT_inv = _mm256_set1_pd(1.0/DT);
    __m256d DR_inv = _mm256_set1_pd(1.0/DR);
    __m256d coe = _mm256_div_pd(_mm256_set1_pd(1.0),_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(sigr,DR_inv),_mm256_mul_pd(sigt,DT_inv)),_mm256_mul_pd(sigp,DP_inv)));

    // Htau  = sqrt(v_fac_a * sqrt(v_T0r * v_tau + v_T0v * (v_pr1 + v_pr2)*0.5));
    __m256d Htau = _mm256_sqrt_pd(_mm256_mul_pd(v_fac_a,_mm256_sqrt_pd(_mm256_add_pd(_mm256_mul_pd(v_T0r,v_tau),_mm256_mul_pd(v_T0v,_mm256_mul_pd(_mm256_add_pd(v_pr1,v_pr2),_mm256_set1_pd(0.5)))))));

    // Htau += sqrt(v_fac_b * sqrt(v_T0r * v_tau + v_T0v * (v_pt1 + v_pt2)*0.5))
    Htau = _mm256_add_pd(Htau,_mm256_sqrt_pd(_mm256_mul_pd(v_fac_b,_mm256_sqrt_pd(_mm256_add_pd(_mm256_mul_pd(v_T0r,v_tau),_mm256_mul_pd(v_T0v,_mm256_mul_pd(_mm256_add_pd(v_pt1,v_pt2),_mm256_set1_pd(0.5))))))));
    // Htau += sqrt(v_fac_c * sqrt(v_T0p * v_tau + v_T0v * (v_pp1 + v_pp2)*0.5))
    Htau = _mm256_add_pd(Htau,_mm256_sqrt_pd(_mm256_mul_pd(v_fac_c,_mm256_sqrt_pd(_mm256_add_pd(_mm256_mul_pd(v_T0p,v_tau),_mm256_mul_pd(v_T0v,_mm256_mul_pd(_mm256_add_pd(v_pp1,v_pp2),_mm256_set1_pd(0.5))))))));

    // tmp1 = ( v_T0t * v_tau + v_T0v * (v_pt1 + v_pt2)*0.5)
    __m256d tmp1 = _mm256_add_pd(_mm256_mul_pd(v_T0t,v_tau),_mm256_mul_pd(v_T0v,_mm256_mul_pd(_mm256_add_pd(v_pt1,v_pt2),_mm256_set1_pd(0.5))));
    // tmp2 = ( v_T0p * v_tau + v_T0 * (v_pp1 + v_pp2)*0.5))
    __m256d tmp2 = _mm256_add_pd(_mm256_mul_pd(v_T0p,v_tau),_mm256_mul_pd(v_T0v,_mm256_mul_pd(_mm256_add_pd(v_pp1,v_pp2),_mm256_set1_pd(0.5))));
    // Htau -= 2.0 * v_fac_f * tmp1 * tmp2;
    Htau = _mm256_sub_pd(Htau,_mm256_mul_pd(_mm256_set1_pd(2.0),_mm256_mul_pd(v_fac_f,_mm256_mul_pd(tmp1,tmp2))));

    // tmp = (sigr*(v_pr1 - v_pr2) + sigt*(v_pt1 - v_pt2) + sigz*(v_pp1 - v_pp2))*0.5
    __m256d tmp = _mm256_mul_pd(_mm256_set1_pd(0.5),_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(sigr,_mm256_sub_pd(v_pr1,v_pr2)),_mm256_mul_pd(sigt,_mm256_sub_pd(v_pt1,v_pt2))),_mm256_mul_pd(sigp,_mm256_sub_pd(v_pp1,v_pp2))));

    /// v_tau += coe * ((v_fun - Htau) + tmp);
    v_tau = _mm256_add_pd(v_tau,_mm256_mul_pd(coe,_mm256_add_pd(_mm256_sub_pd(v_fun,Htau),tmp)));

}

// tau fac_a fac_b fac_c fac_f fun T0v T0r T0t T0p

inline void load_stencil_data(CUSTOMREAL* tau, CUSTOMREAL* fac_a, CUSTOMREAL* fac_b, CUSTOMREAL* fac_c, CUSTOMREAL* fac_f, \
                        CUSTOMREAL* T0v, CUSTOMREAL* T0r, CUSTOMREAL* T0t, CUSTOMREAL* T0p, CUSTOMREAL* fun, \
                        int& iip, int& jjt, int& kkr, int& i_, \
                        CUSTOMREAL* dump_tau, CUSTOMREAL* dump_fac_a, CUSTOMREAL* dump_fac_b, CUSTOMREAL* dump_fac_c, CUSTOMREAL* dump_fac_f, \
                        CUSTOMREAL* dump_T0v, CUSTOMREAL* dump_T0r, CUSTOMREAL* dump_T0t, CUSTOMREAL* dump_T0p, CUSTOMREAL* dump_fun, \
                        CUSTOMREAL* dump_p__, CUSTOMREAL* dump_m__, CUSTOMREAL* dump__p_, CUSTOMREAL* dump__m_, CUSTOMREAL* dump___p, CUSTOMREAL* dump___m, \
                        CUSTOMREAL* dump_pp____, CUSTOMREAL* dump_mm____, CUSTOMREAL* dump___pp__, CUSTOMREAL* dump___mm__, CUSTOMREAL* dump_____pp, CUSTOMREAL* dump_____mm, \
                        int& NP, int& NT, int& NR) {

    if(iip != 0 && jjt != 0 && kkr != 0 && iip != NP-1 && jjt != NT-1 && kkr != NR-1){
        dump_tau[i_]   = tau[I2V(iip,jjt,kkr)];
        dump_fac_a[i_] = fac_a[I2V(iip,jjt,kkr)];
        dump_fac_b[i_] = fac_b[I2V(iip,jjt,kkr)];
        dump_fac_c[i_] = fac_c[I2V(iip,jjt,kkr)];
        dump_fac_f[i_] = fac_f[I2V(iip,jjt,kkr)];
        dump_T0v[i_]   = T0v[I2V(iip,jjt,kkr)];
        dump_T0r[i_]   = T0r[I2V(iip,jjt,kkr)];
        dump_T0t[i_]   = T0t[I2V(iip,jjt,kkr)];
        dump_T0p[i_]   = T0p[I2V(iip,jjt,kkr)];
        dump_fun[i_]   = fun[I2V(iip,jjt,kkr)];
    }

    if (iip == 0) {}
    else if (iip == 1){
        dump_p__[i_]    = tau[I2V(iip+1,jjt,kkr)];
        dump_m__[i_]    = tau[I2V(iip-1,jjt,kkr)];
        dump_pp____[i_] = tau[I2V(iip+2,jjt,kkr)];
        dump_mm____[i_] = 0.0; // dummy
    } else if (iip < NP-2) {
        dump_p__[i_]    = tau[I2V(iip+1,jjt,kkr)];
        dump_m__[i_]    = tau[I2V(iip-1,jjt,kkr)];
        dump_pp____[i_] = tau[I2V(iip+2,jjt,kkr)];
        dump_mm____[i_] = tau[I2V(iip-2,jjt,kkr)];
    } else if (iip == NP-2){
        dump_p__[i_]    = tau[I2V(iip+1,jjt,kkr)];
        dump_m__[i_]    = tau[I2V(iip-1,jjt,kkr)];
        dump_pp____[i_] = 0.0; // dummy
        dump_mm____[i_] = tau[I2V(iip-2,jjt,kkr)];
    }

    if (jjt == 0) {}
    else if (jjt == 1){
        dump__p_[i_]    = tau[I2V(iip,jjt+1,kkr)];
        dump__m_[i_]    = tau[I2V(iip,jjt-1,kkr)];
        dump___pp__[i_] = tau[I2V(iip,jjt+2,kkr)];
        dump___mm__[i_] = 0.0; // dummy
    } else if (jjt < NT-2) {
        dump__p_[i_]    = tau[I2V(iip,jjt+1,kkr)];
        dump__m_[i_]    = tau[I2V(iip,jjt-1,kkr)];
        dump___pp__[i_] = tau[I2V(iip,jjt+2,kkr)];
        dump___mm__[i_] = tau[I2V(iip,jjt-2,kkr)];
    } else if (jjt == NT-2){
        dump__p_[i_]    = tau[I2V(iip,jjt+1,kkr)];
        dump__m_[i_]    = tau[I2V(iip,jjt-1,kkr)];
        dump___pp__[i_] = 0.0; // dummy
        dump___mm__[i_] = tau[I2V(iip,jjt-2,kkr)];
    }

    if (kkr == 0){}
    else if (kkr == 1){
        dump___p[i_]    = tau[I2V(iip,jjt,kkr+1)];
        dump___m[i_]    = tau[I2V(iip,jjt,kkr-1)];
        dump_____pp[i_] = tau[I2V(iip,jjt,kkr+2)];
        dump_____mm[i_] = 0.0; // dummy
    } else if (kkr < NR-2) {
        dump___p[i_]    = tau[I2V(iip,jjt,kkr+1)];
        dump___m[i_]    = tau[I2V(iip,jjt,kkr-1)];
        dump_____pp[i_] = tau[I2V(iip,jjt,kkr+2)];
        dump_____mm[i_] = tau[I2V(iip,jjt,kkr-2)];
    } else if (kkr == NR-2){
        dump___p[i_]    = tau[I2V(iip,jjt,kkr+1)];
        dump___m[i_]    = tau[I2V(iip,jjt,kkr-1)];
        dump_____pp[i_] = 0.0; // dummy
        dump_____mm[i_] = tau[I2V(iip,jjt,kkr-2)];
    }
}


#endif // USE_AVX
#endif //
