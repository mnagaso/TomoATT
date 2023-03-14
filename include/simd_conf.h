#ifndef SIMD_CONF_H
#define SIMD_CONF_H

#ifdef USE_SIMD

// flags for self controling the SIMD usage
#define USE_AVX     __AVX__
#define USE_AVX512  __AVX512F__
#define USE_ARM_SVE __ARM_FEATURE_SVE

// forcely make USE_AVX and USE_AVX512 to be false if USE_ARM_SVE is true
#if USE_ARM_SVE
#undef USE_AVX
#define USE_AVX     0
#undef USE_AVX512
#define USE_AVX512  0
#endif

#if USE_AVX || USE_AVX512
#include <immintrin.h>
#elif USE_ARM_SVE
#include <arm_sve.h>
#endif // __


#ifdef SINGLE_PRECISION

#if USE_AVX && !USE_AVX512
const int ALIGN = 32;
const int NSIMD = 8;
#define __mT __m256
#define _mmT_set1_pT _mm256_set1_ps
#define _mmT_loadu_pT _mm256_loadu_ps
#define _mmT_mul_pT _mm256_mul_ps
#define _mmT_sub_pT _mm256_sub_ps
#define _mmT_add_pT _mm256_add_ps
#define _mmT_div_pT _mm256_div_ps
#define _mmT_min_pT _mm256_min_ps
#define _mmT_sqrt_pT _mm256_sqrt_ps
#define _mmT_store_pT _mm256_store_ps
#define _mmT_fmadd_pT _mm256_fmadd_ps

#define _mm256_cmp_pT _mm256_cmp_ps

#elif USE_AVX512
const int ALIGN = 64;
const int NSIMD = 16;
#define __mT __m512
#define _mmT_set1_pT _mm512_set1_ps
#define _mmT_loadu_pT _mm512_loadu_ps
#define _mmT_mul_pT _mm512_mul_ps
#define _mmT_sub_pT _mm512_sub_ps
#define _mmT_add_pT _mm512_add_ps
#define _mmT_div_pT _mm512_div_ps
#define _mmT_min_pT _mm512_min_ps
#define _mmT_max_pT _mm512_max_ps
#define _mmT_sqrt_pT _mm512_sqrt_ps
#define _mmT_store_pT _mm512_store_ps
#define _mmT_fmadd_pT _mm512_fmadd_ps

#define __mmaskT __mmask16
#define _mm512_cmp_pT_mask _mm512_cmp_ps_mask
#define _mm512_mask_blend_pT _mm512_mask_blend_ps
#define _kand_maskT _kand_mask16


#elif USE_ARM_SVE // NOT TESTED YET
const int ALIGN = 8*svcntd(); // use svcntw() for float
const int NSIMD = svcntd(); // Vector Length
#define __mT svfloat32_t
//#define _mmT_set1_pT svdup_f64
//#define _mmT_loadu_pT svld1_f64_f64_f64
//#define _mmT_mul_pT svmul_f64_m
//#define _mmT_sub_pT svsub_f64_m
//#define _mmT_add_pT svadd_f64_m
//#define _mmT_div_pT svdiv_f64_m
//#define _mmT_min_pT svmin_f64_m
//#define _mmT_sqrt_pT svsqrt_f64_m
//#define _mmT_store_pT svst1_f64

#endif // __ARM_FEATURE_SVE


#else // DOUBLE_PRECISION

#if USE_AVX && !USE_AVX512
const int ALIGN = 32;
const int NSIMD = 4;
#define __mT __m256d
#define _mmT_set1_pT _mm256_set1_pd
#define _mmT_loadu_pT _mm256_loadu_pd
#define _mmT_mul_pT _mm256_mul_pd
#define _mmT_sub_pT _mm256_sub_pd
#define _mmT_add_pT _mm256_add_pd
#define _mmT_div_pT _mm256_div_pd
#define _mmT_min_pT _mm256_min_pd
#define _mmT_sqrt_pT _mm256_sqrt_pd
#define _mmT_store_pT _mm256_store_pd
#define _mmT_fmadd_pT _mm256_fmadd_pd

#define _mm256_cmp_pT _mm256_cmp_pd

#elif USE_AVX512
const int ALIGN = 64;
const int NSIMD = 8;
#define __mT __m512d
#define _mmT_set1_pT _mm512_set1_pd
#define _mmT_loadu_pT _mm512_loadu_pd
#define _mmT_mul_pT _mm512_mul_pd
#define _mmT_sub_pT _mm512_sub_pd
#define _mmT_add_pT _mm512_add_pd
#define _mmT_div_pT _mm512_div_pd
#define _mmT_min_pT _mm512_min_pd
#define _mmT_max_pT _mm512_max_pd
#define _mmT_sqrt_pT _mm512_sqrt_pd
#define _mmT_store_pT _mm512_store_pd
#define _mmT_fmadd_pT _mm512_fmadd_pd

#define __mmaskT __mmask8
#define _mm512_cmp_pT_mask _mm512_cmp_pd_mask
#define _mm512_mask_blend_pT _mm512_mask_blend_pd
#define _kand_maskT _kand_mask8


#elif USE_ARM_SVE
const int ALIGN = 8*svcntd(); // use svcntw() for float
const int NSIMD = svcntd(); // Vector Length
#define __mT svfloat64_t
//#define _mmT_set1_pT svdup_f64
//#define _mmT_loadu_pT svld1_f64_f64_f64
//#define _mmT_mul_pT svmul_f64_m
//#define _mmT_sub_pT svsub_f64_m
//#define _mmT_add_pT svadd_f64_m
//#define _mmT_div_pT svdiv_f64_m
//#define _mmT_min_pT svmin_f64_m
//#define _mmT_sqrt_pT svsqrt_f64_m
//#define _mmT_store_pT svst1_f64

#endif // __ARM_FEATURE_SVE

#endif // DOUBLE_PRECISION

inline void print_simd_type() {
    std::cout << "SIMD type: ";
#if USE_AVX && !USE_AVX512
    std::cout << "AVX" << std::endl;
#elif USE_AVX512
    std::cout << "AVX512" << std::endl;
#elif USE_ARM_SVE
    std::cout << "ARM SVE" << std::endl;
#endif // __ARM_FEATURE_SVE
}

#else // USE_CUDA but not AVX dummy
const int ALIGN = 32;
const int NSIMD = 4;
#endif // USE_SIMD



#endif // SIMD_CONF_H