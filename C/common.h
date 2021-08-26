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

// TODO: document this function
void bitreverse(T *src, size_t n);
unsigned int log_base(unsigned int x, unsigned int n);

#endif