/**
 * @file common.h
 * @brief Common code for pre-computation.
 *
 * Defines T datatype which is used throughout as the coefficient datatype.
 */
#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>
#include <stdio.h>
#include <math.h>

#if !defined(_16BIT_COEFFICIENTS) && !defined(_32BIT_COEFFICIENTS)
#define _16BIT_COEFFICIENTS
# pragma message("default to 16-bit coefficients")
#endif

#if defined(_16BIT_COEFFICIENTS)
typedef uint16_t T;
typedef uint32_t T2;
#elif defined(_32BIT_COEFFICIENTS)
typedef uint32_t T;
typedef uint64_t T2;
#endif

/**
 * @brief Bitreverse an array of length n inplace
 *
 * @param src array
 * @param n length of array
 */
void bitreverse(T *src, size_t n);

/**
 * @brief Computes logarithm with arbitrary base
 *
 * @param x input
 * @param n base
 * @return unsigned int log_n(x)
 */
unsigned int log_base(unsigned int x, unsigned int n);

#endif
