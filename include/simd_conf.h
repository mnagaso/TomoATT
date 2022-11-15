#ifndef SIMD_CONF_H
#define SIMD_CONF_H

#ifdef USE_SIMD

// flags for self controling the SIMD usage
#define USE_AVX     __AVX__
#define USE_AVX512  __AVX512F__
#define USE_ARM_SVE __ARM_FEATURE_SVE

// forcely make USE_AVX and USE_AVX512 to be false if USE_ARM_SVE is true
#if USE_ARM_SVE
#define USE_AVX     0
#define USE_AVX512  0
#endif

#if USE_AVX || USE_AVX512
#include <immintrin.h>
#elif USE_ARM_SVE
#include <arm_sve.h>
#endif // __


#if USE_AVX && !USE_AVX512
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

#elif USE_AVX512
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

#elif USE_ARM_SVE
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