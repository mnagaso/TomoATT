#ifndef SIMD_CONF_H
#define SIMD_CONF_H

#ifdef USE_SIMD


#if defined __AVX__ || defined __AVX512F__
#include <immintrin.h>
#elif defined __ARM_FEATURE_SVE
#include <arm_sve.h>
#endif // __


#ifdef __AVX__
const int ALIGN = 32;
const int NSIMD = 4;
#define __mTd __m256d
#define _mmT_set1_pd _mm256_set1_pd
#define _mmT_loadu_pd _mm256_loadu_pd
#define _mmT_mul_pd _mm256_mul_pd
#define _mmT_sub_pd _mm256_sub_pd
#define _mmT_add_pd _mm256_add_pd
#define _mmT_div_pd _mm256_div_pd
#define _mmT_min_pd _mm256_min_pd
#define _mmT_sqrt_pd _mm256_sqrt_pd
#define _mmT_store_pd _mm256_store_pd
#define _mmT_fmadd_pd _mm256_fmadd_pd

#elif defined __AVX512F__
const int ALIGN = 64;
const int NSIMD = 8;
#define __mTd __m512d
#define _mmT_set1_pd _mm512_set1_pd
#define _mmT_loadu_pd _mm512_loadu_pd
#define _mmT_mul_pd _mm512_mul_pd
#define _mmT_sub_pd _mm512_sub_pd
#define _mmT_add_pd _mm512_add_pd
#define _mmT_div_pd _mm512_div_pd
#define _mmT_min_pd _mm512_min_pd
#define _mmT_sqrt_pd _mm512_sqrt_pd
#define _mmT_store_pd _mm512_store_pd
#define _mmT_fmadd_pd _mm512_fmadd_pd

#elif defined __ARM_FEATURE_SVE
const int ALIGN = 8*svcntd(); // use svcntw() for float
const int NSIMD = svcntd(); // Vector Length
#define __mTd svfloat64_t
//#define _mmT_set1_pd svdup_f64
//#define _mmT_loadu_pd svld1_f64_f64_f64
//#define _mmT_mul_pd svmul_f64_m
//#define _mmT_sub_pd svsub_f64_m
//#define _mmT_add_pd svadd_f64_m
//#define _mmT_div_pd svdiv_f64_m
//#define _mmT_min_pd svmin_f64_m
//#define _mmT_sqrt_pd svsqrt_f64_m
//#define _mmT_store_pd svst1_f64

#endif // __ARM_FEATURE_SVE


inline void print_simd_type() {
    std::cout << "SIMD type: ";
#ifdef __AVX__
    std::cout << "AVX" << std::endl;
#elif defined __AVX512F__
    std::cout << "AVX512" << std::endl;
#elif defined __ARM_FEATURE_SVE
    std::cout << "ARM SVE" << std::endl;
#endif // __ARM_FEATURE_SVE
}


#endif // USE_SIMD

#endif // SIMD_CONF_H