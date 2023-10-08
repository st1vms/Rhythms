#ifndef MX_UTILS_H

#ifdef __cplusplus
extern "C"
{
#endif

#include "global_types.h"
#include <stddef.h>

extern unsigned char getMatrixDeterminant(size_t rows, size_t cols, matrix_t matrix[rows][cols], det_t *out);

extern unsigned char matrixMultiply(size_t r1, size_t c1, matrix_t m1[r1][c1], size_t r2, size_t c2, matrix_t m2[r2][c2], matrix_t out[r1][c2]);

extern unsigned char vectorMultiply(size_t s1, matrix_t v1[s1], size_t s2, matrix_t v2[s2], det_t *out);

#ifdef __cplusplus
}
#endif

#define MX_UTILS_H
#endif // MX_UTILS_H
