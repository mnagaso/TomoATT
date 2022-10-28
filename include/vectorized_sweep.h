#ifndef VECTORIZED_SWEEP_H
#define VECTORIZED_SWEEP_H


#ifdef USE_AVX
#include <immintrin.h>


#include<vector>
#include "config.h"


__mTd COEF    = _mmT_set1_pd(1.0);
__mTd v_1     = _mmT_set1_pd(1.0);
__mTd v_0     = _mmT_set1_pd(0.0);
__mTd v_half  = _mmT_set1_pd(0.5);
__mTd v_2     = _mmT_set1_pd(2.0);
__mTd v_m2    = _mmT_set1_pd(-2.0);
__mTd v_4     = _mmT_set1_pd(4.0);
__mTd v_m3    = _mmT_set1_pd(-3.0);
__mTd coe_max = _mmT_set1_pd(1e19);
__mTd v_eps   = _mmT_set1_pd(1e-12);

//#define CUSTOMREAL double

//int NP, NT, NR;
//#define I2V(i,j,k) ((k)*NP*NT+(j)*NP+(i))

// square of __mTd
inline __mTd my_square(__mTd const& a){
    return _mmT_mul_pd(a, a);
}

// pp = (a - b) * Dinv
inline __mTd calc_1d_stencil(__mTd const& a, __mTd const& b, __mTd const& Dinv){
    return _mmT_mul_pd(_mmT_sub_pd(a,b),Dinv);
}

/*
    ww  = _1_CR/(_1_CR+_2_CR*my_square((eps + my_square(      a \
                                                       -_2_CR*b \
                                                       +      c) ) \

                                     / (eps + my_square(      d \
                                                       -_2_CR*a \
                                                       +      b) )) );

    pp  = sign  * (_1_CR - ww) * (         b \                  <-- sign = +1 or -1
                                 -         d) * 0.5 * Dinv \
                         + ww  * ( -3.0 * a \
                                 +  4.0 * b \
                                 -  1.0 * c ) * 0.5 * Dinv );
*/
inline __mTd calc_3d_stencil(__mTd const& a, __mTd const& b, __mTd const&c, __mTd const& d, __mTd const& Dinv_half, int& sign){

    // v_eps + sqrt(a - 2.0*b + c)
    __mTd tmp1 = _mmT_add_pd(v_eps,my_square(_mmT_add_pd(a,_mmT_add_pd(_mmT_mul_pd(v_m2,b),c))));
    // v_eps + sqrt(d - 2.0*a + b)
    __mTd tmp2 = _mmT_add_pd(v_eps,my_square(_mmT_add_pd(d,_mmT_add_pd(_mmT_mul_pd(v_m2,a),b))));
    // ww = 1.0/(1.0 + 2.0 * square(tmp1/tmp2))
    __mTd ww = _mmT_div_pd(v_1,_mmT_add_pd(v_1,_mmT_mul_pd(v_2,my_square(_mmT_div_pd(tmp1,tmp2)))));
    // pp = sign* ((1.0 - ww) * (d - b) / 2.0 / D
    //                  + ww  * (-3.0* a + 4.0 * b - c) * Dinv_half)
    return _mmT_mul_pd(_mmT_set1_pd(sign), \
                _mmT_add_pd(\
                            _mmT_mul_pd(_mmT_sub_pd(v_1,ww),_mmT_mul_pd(_mmT_sub_pd(b,d),Dinv_half)),\
                            _mmT_mul_pd(ww,_mmT_mul_pd(_mmT_sub_pd(_mmT_add_pd(_mmT_mul_pd(v_4,b),_mmT_mul_pd(v_m3,a)),c),Dinv_half))\
                )\
           );


}


inline void fake_stencil_3rd_pre_simd(__mTd& v_iip, __mTd& v_jjt, __mTd& v_kkr, __mTd& v_c__,  __mTd& v_p__,     __mTd& v_m__,     __mTd& v__p_,    __mTd& v__m_,    __mTd& v___p,    __mTd& v___m, \
                                                                                                       __mTd& v_pp____,  __mTd& v_mm____,  __mTd& v___pp__, __mTd& v___mm__, __mTd& v_____pp, __mTd& v_____mm, \
                                                                                                       __mTd& v_pp1,     __mTd& v_pp2,     __mTd& v_pt1,    __mTd& v_pt2,    __mTd& v_pr1,    __mTd& v_pr2, \
                                                                                                       __mTd& v_DP_inv,  __mTd& v_DT_inv,  __mTd& v_DR_inv, \
                                                                                                       __mTd& v_DP_inv_half, __mTd& v_DT_inv_half, __mTd& v_DR_inv_half, \
                                                                                                       int& NP, int& NT, int& NR){

    int PLUS  = 1;
    int MINUS = -1;

    // TODO: this is only for AVX512
#ifdef USE_AVX512
    __mmask8 mask_i_eq_1         = _mm512_cmp_pd_mask(v_iip, v_1,_CMP_EQ_OQ);    // if iip == 1
    __mmask8 mask_j_eq_1         = _mm512_cmp_pd_mask(v_jjt, v_1,_CMP_EQ_OQ);    // if jjt == 1
    __mmask8 mask_k_eq_1         = _mm512_cmp_pd_mask(v_kkr, v_1,_CMP_EQ_OQ);    // if kkr == 1
    __mmask8 mask_i_eq_N_minus_2 = _mm512_cmp_pd_mask(v_iip, _mm512_set1_pd(NP-2),_CMP_EQ_OQ); // if iip == N-2
    __mmask8 mask_j_eq_N_minus_2 = _mm512_cmp_pd_mask(v_jjt, _mm512_set1_pd(NT-2),_CMP_EQ_OQ); // if jjt == N-2
    __mmask8 mask_k_eq_N_minus_2 = _mm512_cmp_pd_mask(v_kkr, _mm512_set1_pd(NR-2),_CMP_EQ_OQ); // if kkr == N-2

    // 1 < iip < N-2
    __mmask8 mask_i_else = _kand_mask8(
                               _mm512_cmp_pd_mask(v_iip, v_1                , _CMP_GT_OQ),
                               _mm512_cmp_pd_mask(_mm512_set1_pd(NP-2),  v_iip, _CMP_GT_OQ)
                         );
    // 1 < jjt < N-2
    __mmask8 mask_j_else = _kand_mask8(
                               _mm512_cmp_pd_mask(v_jjt, v_1                , _CMP_GT_OQ),
                               _mm512_cmp_pd_mask(_mm512_set1_pd(NT-2),  v_jjt, _CMP_GT_OQ)
                         );
    // 1 < kkr < N-2
    __mmask8 mask_k_else = _kand_mask8(
                               _mm512_cmp_pd_mask(v_kkr, v_1                , _CMP_GT_OQ),
                               _mm512_cmp_pd_mask(_mm512_set1_pd(NR-2),  v_kkr, _CMP_GT_OQ)
                         );

    // if _i_eq_1 == true
    v_pp1 = _mm512_mask_blend_pd(mask_i_eq_1, calc_1d_stencil(v_c__, v_m__, v_DP_inv)                            , v_pp1 );
    v_pp2 = _mm512_mask_blend_pd(mask_i_eq_1, calc_3d_stencil(v_c__, v_p__, v_pp____, v_m__, v_DP_inv_half, PLUS), v_pp2 );
    // if_i_eq_N_minus_2 == true
    v_pp1 = _mm512_mask_blend_pd(mask_i_eq_N_minus_2, calc_3d_stencil(v_c__, v_m__, v_mm____, v_p__, v_DP_inv_half, MINUS), v_pp1);
    v_pp2 = _mm512_mask_blend_pd(mask_i_eq_N_minus_2, calc_1d_stencil(v_p__, v_c__, v_DP_inv)                             , v_pp2);
    // else
    v_pp1 = _mm512_mask_blend_pd(mask_i_else, calc_3d_stencil(v_c__, v_m__, v_mm____, v_p__, v_DP_inv_half, MINUS), v_pp1);
    v_pp2 = _mm512_mask_blend_pd(mask_i_else, calc_3d_stencil(v_c__, v_p__, v_pp____, v_m__, v_DP_inv_half, PLUS) , v_pp2);

    // if _j_eq_1 == true
    v_pt1 = _mm512_mask_blend_pd(mask_j_eq_1, calc_1d_stencil(v_c__, v__m_, v_DT_inv)                            , v_pt1);
    v_pt2 = _mm512_mask_blend_pd(mask_j_eq_1, calc_3d_stencil(v_c__, v__p_, v___pp__, v__m_, v_DT_inv_half, PLUS), v_pt2);
    // if _j_eq_N_minus_2 == true
    v_pt1 = _mm512_mask_blend_pd(mask_j_eq_N_minus_2, calc_3d_stencil(v_c__, v__m_, v___mm__, v__p_, v_DT_inv_half, MINUS), v_pt1);
    v_pt2 = _mm512_mask_blend_pd(mask_j_eq_N_minus_2, calc_1d_stencil(v__p_, v_c__, v_DT_inv)                             , v_pt2);
    // else
    v_pt1 = _mm512_mask_blend_pd(mask_j_else, calc_3d_stencil(v_c__, v__m_, v___mm__, v__p_, v_DT_inv_half, MINUS), v_pt1);
    v_pt2 = _mm512_mask_blend_pd(mask_j_else, calc_3d_stencil(v_c__, v__p_, v___pp__, v__m_, v_DT_inv_half, PLUS) , v_pt2);

    // if _k_eq_1 == true
    v_pr1 = _mm512_mask_blend_pd(mask_k_eq_1, calc_1d_stencil(v_c__, v___m, v_DR_inv)                            , v_pr1 );
    v_pr2 = _mm512_mask_blend_pd(mask_k_eq_1, calc_3d_stencil(v_c__, v___p, v_____pp, v___m, v_DR_inv_half, PLUS), v_pr2 );
    // if _k_eq_N_minus_2 == true
    v_pr1 = _mm512_mask_blend_pd(mask_k_eq_N_minus_2, calc_3d_stencil(v_c__, v___m, v_____mm, v___p, v_DR_inv_half, MINUS), v_pr1);
    v_pr2 = _mm512_mask_blend_pd(mask_k_eq_N_minus_2, calc_1d_stencil(v___p, v_c__, v_DR_inv)                             , v_pr2);
    // else
    v_pr1 = _mm512_mask_blend_pd(mask_k_else, calc_3d_stencil(v_c__, v___m, v_____mm, v___p, v_DR_inv_half, MINUS), v_pr1);
    v_pr2 = _mm512_mask_blend_pd(mask_k_else, calc_3d_stencil(v_c__, v___p, v_____pp, v___m, v_DR_inv_half, PLUS) , v_pr2);

#endif

#ifdef USE_AVX
    __m256d mask_i_eq_1         = _mm256_cmp_pd(v_iip, v_1,_CMP_EQ_OQ);    // if iip == 1
    __m256d mask_j_eq_1         = _mm256_cmp_pd(v_jjt, v_1,_CMP_EQ_OQ);    // if jjt == 1
    __m256d mask_k_eq_1         = _mm256_cmp_pd(v_kkr, v_1,_CMP_EQ_OQ);    // if kkr == 1
    __m256d mask_i_eq_N_minus_2 = _mm256_cmp_pd(v_iip, _mm256_set1_pd(NP-2),_CMP_EQ_OQ); // if iip == N-2
    __m256d mask_j_eq_N_minus_2 = _mm256_cmp_pd(v_jjt, _mm256_set1_pd(NT-2),_CMP_EQ_OQ); // if jjt == N-2
    __m256d mask_k_eq_N_minus_2 = _mm256_cmp_pd(v_kkr, _mm256_set1_pd(NR-2),_CMP_EQ_OQ); // if kkr == N-2

    // 1 < iip < N-2
    __m256d mask_i_else = _mm256_and_pd(
                                        _mm256_cmp_pd(v_iip, v_1                , _CMP_GT_OQ),
                                        _mm256_cmp_pd(_mm256_set1_pd(NP-2),v_iip, _CMP_GT_OQ)
                         );
    // 1 < jjt < N-2
    __m256d mask_j_else = _mm256_and_pd(
                                        _mm256_cmp_pd(v_jjt, v_1                , _CMP_GT_OQ),
                                        _mm256_cmp_pd(_mm256_set1_pd(NT-2),v_jjt, _CMP_GT_OQ)
                         );
    // 1 < kkr < N-2
    __m256d mask_k_else = _mm256_and_pd(
                                        _mm256_cmp_pd(v_kkr, v_1                , _CMP_GT_OQ),
                                        _mm256_cmp_pd(_mm256_set1_pd(NR-2),v_kkr, _CMP_GT_OQ)
                         );

    // if _i_eq_1 == true
    v_pp1 = _mm256_blendv_pd(v_pp1, calc_1d_stencil(v_c__, v_m__, v_DP_inv),                             mask_i_eq_1);
    v_pp2 = _mm256_blendv_pd(v_pp2, calc_3d_stencil(v_c__, v_p__, v_pp____, v_m__, v_DP_inv_half, PLUS), mask_i_eq_1);
    // if_i_eq_N_minus_2 == true
    v_pp1 = _mm256_blendv_pd(v_pp1, calc_3d_stencil(v_c__, v_m__, v_mm____, v_p__, v_DP_inv_half, MINUS), mask_i_eq_N_minus_2);
    v_pp2 = _mm256_blendv_pd(v_pp2, calc_1d_stencil(v_p__, v_c__, v_DP_inv),                              mask_i_eq_N_minus_2);
    // else
    v_pp1 = _mm256_blendv_pd(v_pp1, calc_3d_stencil(v_c__, v_m__, v_mm____, v_p__, v_DP_inv_half, MINUS), mask_i_else);
    v_pp2 = _mm256_blendv_pd(v_pp2, calc_3d_stencil(v_c__, v_p__, v_pp____, v_m__, v_DP_inv_half, PLUS),  mask_i_else);

    // if _j_eq_1 == true
    v_pt1 = _mm256_blendv_pd(v_pt1, calc_1d_stencil(v_c__, v__m_, v_DT_inv),                             mask_j_eq_1);
    v_pt2 = _mm256_blendv_pd(v_pt2, calc_3d_stencil(v_c__, v__p_, v___pp__, v__m_, v_DT_inv_half, PLUS), mask_j_eq_1);
    // if _j_eq_N_minus_2 == true
    v_pt1 = _mm256_blendv_pd(v_pt1, calc_3d_stencil(v_c__, v__m_, v___mm__, v__p_, v_DT_inv_half, MINUS), mask_j_eq_N_minus_2);
    v_pt2 = _mm256_blendv_pd(v_pt2, calc_1d_stencil(v__p_, v_c__, v_DT_inv),                              mask_j_eq_N_minus_2);
    // else
    v_pt1 = _mm256_blendv_pd(v_pt1, calc_3d_stencil(v_c__, v__m_, v___mm__, v__p_, v_DT_inv_half, MINUS), mask_j_else);
    v_pt2 = _mm256_blendv_pd(v_pt2, calc_3d_stencil(v_c__, v__p_, v___pp__, v__m_, v_DT_inv_half, PLUS),  mask_j_else);

    // if _k_eq_1 == true
    v_pr1 = _mm256_blendv_pd(v_pr1, calc_1d_stencil(v_c__, v___m, v_DR_inv),                             mask_k_eq_1);
    v_pr2 = _mm256_blendv_pd(v_pr2, calc_3d_stencil(v_c__, v___p, v_____pp, v___m, v_DR_inv_half, PLUS), mask_k_eq_1);
    // if _k_eq_N_minus_2 == true
    v_pr1 = _mm256_blendv_pd(v_pr1, calc_3d_stencil(v_c__, v___m, v_____mm, v___p, v_DR_inv_half, MINUS), mask_k_eq_N_minus_2);
    v_pr2 = _mm256_blendv_pd(v_pr2, calc_1d_stencil(v___p, v_c__, v_DR_inv),                              mask_k_eq_N_minus_2);
    // else
    v_pr1 = _mm256_blendv_pd(v_pr1, calc_3d_stencil(v_c__, v___m, v_____mm, v___p, v_DR_inv_half, MINUS), mask_k_else);
    v_pr2 = _mm256_blendv_pd(v_pr2, calc_3d_stencil(v_c__, v___p, v_____pp, v___m, v_DR_inv_half, PLUS),  mask_k_else);
#endif

}


// tau fac_a fac_b fac_c fac_f T0v T0p T0t T0r fun
inline void fake_stencil_3rd_apre_simd(__mTd& v_tau,  __mTd& v_fac_a,__mTd& v_fac_b, __mTd& v_fac_c, __mTd& v_fac_f, \
                                       __mTd& v_T0v,  __mTd& v_T0p,  __mTd& v_T0t,   __mTd& v_T0r,   __mTd& v_fun, __mTd& v_change,\
                                       __mTd& v_pp1,  __mTd& v_pp2,  __mTd& v_pt1,   __mTd& v_pt2,   __mTd& v_pr1, __mTd& v_pr2, \
                                       __mTd& DP_inv, __mTd& DT_inv, __mTd& DR_inv){

    // sigr = COEF * sqrt(v_fac_a)*v_T0v;
    __mTd sigr = _mmT_mul_pd(_mmT_mul_pd(COEF,_mmT_sqrt_pd(v_fac_a)),v_T0v);
    // sigt = COEF * sqrt(v_fac_b)*v_T0v;
    __mTd sigt = _mmT_mul_pd(_mmT_mul_pd(COEF,_mmT_sqrt_pd(v_fac_b)),v_T0v);
    // sigp = COEF * sqrt(v_fac_c)*v_T0v;
    __mTd sigp = _mmT_mul_pd(_mmT_mul_pd(COEF,_mmT_sqrt_pd(v_fac_c)),v_T0v);

    // coe = 1.0 / (sigr/D + sigt/D + sigz/D);
    __mTd coe = _mmT_div_pd(v_1,_mmT_add_pd(_mmT_add_pd(_mmT_mul_pd(sigr,DR_inv),_mmT_mul_pd(sigt,DT_inv)),_mmT_mul_pd(sigp,DP_inv)));
    // coe becomes inf as sig* goes to 0
    // if coe > 1e19, set coe = 1e19
    coe = _mmT_min_pd(coe,coe_max);

    // Htau  = v_fac_a * square(v_T0r * v_tau + v_T0v * (v_pr1 + v_pr2)*0.5);
    __mTd Htau = _mmT_mul_pd(v_fac_a,my_square(_mmT_add_pd(_mmT_mul_pd(v_T0r,v_tau),_mmT_mul_pd(v_T0v,_mmT_mul_pd(_mmT_add_pd(v_pr1,v_pr2),v_half)))));

    // Htau += v_fac_b * square(v_T0t * v_tau + v_T0v * (v_pt1 + v_pt2)*0.5))
    Htau = _mmT_add_pd(Htau,_mmT_mul_pd(v_fac_b,my_square(_mmT_add_pd(_mmT_mul_pd(v_T0t,v_tau),_mmT_mul_pd(v_T0v,_mmT_mul_pd(_mmT_add_pd(v_pt1,v_pt2), v_half))))));
    // Htau += v_fac_c * square(v_T0p * v_tau + v_T0v * (v_pp1 + v_pp2)*0.5))
    Htau = _mmT_add_pd(Htau,_mmT_mul_pd(v_fac_c,my_square(_mmT_add_pd(_mmT_mul_pd(v_T0p,v_tau),_mmT_mul_pd(v_T0v,_mmT_mul_pd(_mmT_add_pd(v_pp1,v_pp2), v_half))))));

    // tmp1 = ( v_T0t * v_tau + v_T0v * (v_pt1 + v_pt2)*0.5)
    __mTd tmp1 = _mmT_add_pd(_mmT_mul_pd(v_T0t,v_tau),_mmT_mul_pd(v_T0v,_mmT_mul_pd(_mmT_add_pd(v_pt1,v_pt2), v_half)));
    // tmp2 = ( v_T0p * v_tau + v_T0 * (v_pp1 + v_pp2)*0.5))
    __mTd tmp2 = _mmT_add_pd(_mmT_mul_pd(v_T0p,v_tau),_mmT_mul_pd(v_T0v,_mmT_mul_pd(_mmT_add_pd(v_pp1,v_pp2), v_half)));
    // tmp3 = -2.0 * v_fac_f * tmp1 * tmp2;
    __mTd tmp3 = _mmT_mul_pd(v_m2,_mmT_mul_pd(v_fac_f,_mmT_mul_pd(tmp1,tmp2)));

    // Htau = sqrt(Htau + tmp3); // # becamse nan as Htau - tmp3 goes to 0
    Htau = _mmT_sqrt_pd(_mmT_add_pd(Htau,tmp3));

    // tmp = (sigr*(v_pr1 - v_pr2) + sigt*(v_pt1 - v_pt2) + sigz*(v_pp1 - v_pp2))*0.5
    __mTd tmp = _mmT_mul_pd(v_half,_mmT_add_pd(_mmT_add_pd(_mmT_mul_pd(sigr,_mmT_sub_pd(v_pr2,v_pr1)),_mmT_mul_pd(sigt,_mmT_sub_pd(v_pt2,v_pt1))),_mmT_mul_pd(sigp,_mmT_sub_pd(v_pp2,v_pp1))));

    // v_tau += coe * ((v_fun - Htau) + tmp) if mask is true
    v_tau = _mmT_add_pd(v_tau,_mmT_mul_pd(coe,_mmT_add_pd(_mmT_sub_pd(v_fun,Htau),tmp)));

#ifdef USE_AVX512
    // mask if v_change != 1.0
    __mmask8 mask = _mm512_cmp_pd_mask(v_change,v_1,_CMP_NEQ_OQ);
    // set 1 if mask is true
    v_tau = _mm512_mask_blend_pd(mask,v_tau,v_1);
#endif

#ifdef USE_AVX

#endif

}


inline __mTd load_mem_gen_to_mTd(CUSTOMREAL* a, \
                         int* iip, int* jjt, int* kkr){
        const CUSTOMREAL dump_[NSIMD] = {a[I2V(iip[0],jjt[0],kkr[0])], \
                                         a[I2V(iip[1],jjt[1],kkr[1])], \
                                         a[I2V(iip[2],jjt[2],kkr[2])], \
                                         a[I2V(iip[3],jjt[3],kkr[3])]};
        return  _mmT_loadu_pd(dump_);
}

inline __mTd load_mem_bool_to_mTd(bool* a, \
                         int* iip, int* jjt, int* kkr){
        const CUSTOMREAL dump_[NSIMD] = {(CUSTOMREAL) a[I2V(iip[0],jjt[0],kkr[0])], \
                                         (CUSTOMREAL) a[I2V(iip[1],jjt[1],kkr[1])], \
                                         (CUSTOMREAL) a[I2V(iip[2],jjt[2],kkr[2])], \
                                         (CUSTOMREAL) a[I2V(iip[3],jjt[3],kkr[3])]};
        return _mmT_loadu_pd(dump_);
}


/*
inline void load_stencil_tau(CUSTOMREAL* tau, \
                        int& iip, int& jjt, int& kkr, int& i_, \
                        CUSTOMREAL* dump_tau, \
                        CUSTOMREAL* dump_p__, CUSTOMREAL* dump_m__, CUSTOMREAL* dump__p_, CUSTOMREAL* dump__m_, CUSTOMREAL* dump___p, CUSTOMREAL* dump___m, \
                        CUSTOMREAL* dump_pp____, CUSTOMREAL* dump_mm____, CUSTOMREAL* dump___pp__, CUSTOMREAL* dump___mm__, CUSTOMREAL* dump_____pp, CUSTOMREAL* dump_____mm, \
                        int& NP, int& NT, int& NR) {

    dump_tau[i_]   = tau[I2V(iip,jjt,kkr)];
     if (iip <= 0) {}
    else if (iip == 1){
        dump_p__[i_]    = tau[I2V(iip+1,jjt,kkr)];
        dump_m__[i_]    = tau[I2V(iip-1,jjt,kkr)];
        dump_pp____[i_] = tau[I2V(iip+2,jjt,kkr)];
        dump_mm____[i_] = 1.0; // dummy
    } else if (iip < NP-2) {
        dump_p__[i_]    = tau[I2V(iip+1,jjt,kkr)];
        dump_m__[i_]    = tau[I2V(iip-1,jjt,kkr)];
        dump_pp____[i_] = tau[I2V(iip+2,jjt,kkr)];
        dump_mm____[i_] = tau[I2V(iip-2,jjt,kkr)];
    } else if (iip == NP-2){
        dump_p__[i_]    = tau[I2V(iip+1,jjt,kkr)];
        dump_m__[i_]    = tau[I2V(iip-1,jjt,kkr)];
        dump_pp____[i_] = 1.0; // dummy
        dump_mm____[i_] = tau[I2V(iip-2,jjt,kkr)];
    }

    if (jjt <= 0) {}
    else if (jjt == 1){
        dump__p_[i_]    = tau[I2V(iip,jjt+1,kkr)];
        dump__m_[i_]    = tau[I2V(iip,jjt-1,kkr)];
        dump___pp__[i_] = tau[I2V(iip,jjt+2,kkr)];
        dump___mm__[i_] = 1.0; // dummy
    } else if (jjt < NT-2) {
        dump__p_[i_]    = tau[I2V(iip,jjt+1,kkr)];
        dump__m_[i_]    = tau[I2V(iip,jjt-1,kkr)];
        dump___pp__[i_] = tau[I2V(iip,jjt+2,kkr)];
        dump___mm__[i_] = tau[I2V(iip,jjt-2,kkr)];
    } else if (jjt == NT-2){
        dump__p_[i_]    = tau[I2V(iip,jjt+1,kkr)];
        dump__m_[i_]    = tau[I2V(iip,jjt-1,kkr)];
        dump___pp__[i_] = 1.0; // dummy
        dump___mm__[i_] = tau[I2V(iip,jjt-2,kkr)];
    }

    if (kkr <= 0){}
    else if (kkr == 1){
        dump___p[i_]    = tau[I2V(iip,jjt,kkr+1)];
        dump___m[i_]    = tau[I2V(iip,jjt,kkr-1)];
        dump_____pp[i_] = tau[I2V(iip,jjt,kkr+2)];
        dump_____mm[i_] = 1.0; // dummy
    } else if (kkr < NR-2) {
        dump___p[i_]    = tau[I2V(iip,jjt,kkr+1)];
        dump___m[i_]    = tau[I2V(iip,jjt,kkr-1)];
        dump_____pp[i_] = tau[I2V(iip,jjt,kkr+2)];
        dump_____mm[i_] = tau[I2V(iip,jjt,kkr-2)];
    } else if (kkr == NR-2){
        dump___p[i_]    = tau[I2V(iip,jjt,kkr+1)];
        dump___m[i_]    = tau[I2V(iip,jjt,kkr-1)];
        dump_____pp[i_] = 1.0; // dummy
        dump_____mm[i_] = tau[I2V(iip,jjt,kkr-2)];
    }
}

inline void load_mem_gen(CUSTOMREAL* a, CUSTOMREAL* dump, \
                         int& iip, int& jjt, int& kkr, int& i_){
    dump[i_] = a[I2V(iip,jjt,kkr)];
}

inline void load_mem_bool(bool* a, CUSTOMREAL* dump, \
                         int& iip, int& jjt, int& kkr, int& i_){
    dump[i_] = (CUSTOMREAL)a[I2V(iip,jjt,kkr)];
}
*/

/*
// tau fac_a fac_b fac_c fac_f fun T0v T0r T0t T0p
inline void load_stencil_data(CUSTOMREAL* tau, CUSTOMREAL* fac_a, CUSTOMREAL* fac_b, CUSTOMREAL* fac_c, CUSTOMREAL* fac_f, \
                        CUSTOMREAL* T0v, CUSTOMREAL* T0r, CUSTOMREAL* T0t, CUSTOMREAL* T0p, CUSTOMREAL* fun, bool* is_changed, \
                        int& iip, int& jjt, int& kkr, int& i_, \
                        CUSTOMREAL* dump_tau, CUSTOMREAL* dump_fac_a, CUSTOMREAL* dump_fac_b, CUSTOMREAL* dump_fac_c, CUSTOMREAL* dump_fac_f, \
                        CUSTOMREAL* dump_T0v, CUSTOMREAL* dump_T0r, CUSTOMREAL* dump_T0t, CUSTOMREAL* dump_T0p, CUSTOMREAL* dump_fun, CUSTOMREAL* dump_change, \
                        CUSTOMREAL* dump_p__, CUSTOMREAL* dump_m__, CUSTOMREAL* dump__p_, CUSTOMREAL* dump__m_, CUSTOMREAL* dump___p, CUSTOMREAL* dump___m, \
                        CUSTOMREAL* dump_pp____, CUSTOMREAL* dump_mm____, CUSTOMREAL* dump___pp__, CUSTOMREAL* dump___mm__, CUSTOMREAL* dump_____pp, CUSTOMREAL* dump_____mm, \
                        int& NP, int& NT, int& NR) {

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
    // convert bool to CUSTOMREAL
    dump_change[i_] = (CUSTOMREAL)is_changed[I2V(iip,jjt,kkr)];

    if (iip <= 0) {}
    else if (iip == 1){
        dump_p__[i_]    = tau[I2V(iip+1,jjt,kkr)];
        dump_m__[i_]    = tau[I2V(iip-1,jjt,kkr)];
        dump_pp____[i_] = tau[I2V(iip+2,jjt,kkr)];
        dump_mm____[i_] = 1.0; // dummy
    } else if (iip < NP-2) {
        dump_p__[i_]    = tau[I2V(iip+1,jjt,kkr)];
        dump_m__[i_]    = tau[I2V(iip-1,jjt,kkr)];
        dump_pp____[i_] = tau[I2V(iip+2,jjt,kkr)];
        dump_mm____[i_] = tau[I2V(iip-2,jjt,kkr)];
    } else if (iip == NP-2){
        dump_p__[i_]    = tau[I2V(iip+1,jjt,kkr)];
        dump_m__[i_]    = tau[I2V(iip-1,jjt,kkr)];
        dump_pp____[i_] = 1.0; // dummy
        dump_mm____[i_] = tau[I2V(iip-2,jjt,kkr)];
    }

    if (jjt <= 0) {}
    else if (jjt == 1){
        dump__p_[i_]    = tau[I2V(iip,jjt+1,kkr)];
        dump__m_[i_]    = tau[I2V(iip,jjt-1,kkr)];
        dump___pp__[i_] = tau[I2V(iip,jjt+2,kkr)];
        dump___mm__[i_] = 1.0; // dummy
    } else if (jjt < NT-2) {
        dump__p_[i_]    = tau[I2V(iip,jjt+1,kkr)];
        dump__m_[i_]    = tau[I2V(iip,jjt-1,kkr)];
        dump___pp__[i_] = tau[I2V(iip,jjt+2,kkr)];
        dump___mm__[i_] = tau[I2V(iip,jjt-2,kkr)];
    } else if (jjt == NT-2){
        dump__p_[i_]    = tau[I2V(iip,jjt+1,kkr)];
        dump__m_[i_]    = tau[I2V(iip,jjt-1,kkr)];
        dump___pp__[i_] = 1.0; // dummy
        dump___mm__[i_] = tau[I2V(iip,jjt-2,kkr)];
    }

    if (kkr <= 0){}
    else if (kkr == 1){
        dump___p[i_]    = tau[I2V(iip,jjt,kkr+1)];
        dump___m[i_]    = tau[I2V(iip,jjt,kkr-1)];
        dump_____pp[i_] = tau[I2V(iip,jjt,kkr+2)];
        dump_____mm[i_] = 1.0; // dummy
    } else if (kkr < NR-2) {
        dump___p[i_]    = tau[I2V(iip,jjt,kkr+1)];
        dump___m[i_]    = tau[I2V(iip,jjt,kkr-1)];
        dump_____pp[i_] = tau[I2V(iip,jjt,kkr+2)];
        dump_____mm[i_] = tau[I2V(iip,jjt,kkr-2)];
    } else if (kkr == NR-2){
        dump___p[i_]    = tau[I2V(iip,jjt,kkr+1)];
        dump___m[i_]    = tau[I2V(iip,jjt,kkr-1)];
        dump_____pp[i_] = 1.0; // dummy
        dump_____mm[i_] = tau[I2V(iip,jjt,kkr-2)];
    }
}
*/

#endif // USE_AVX or USE_AVX512



#endif //