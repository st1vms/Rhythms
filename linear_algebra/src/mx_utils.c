#include <stdio.h>
#include <stdlib.h>
#include "mx_utils.h"

static unsigned char createSubmatrix(size_t r1, size_t c1, matrix_t matrix[r1][c1], size_t r2, size_t c2, matrix_t out[r2][c2], unsigned long int row, unsigned long int col)
{
    unsigned long int sub_i = 0, sub_j = 0;
    for (unsigned long int i = 0; i < r1; ++i)
    {
        if (i != row)
        {
            sub_j = 0;
            for (unsigned long int j = 0; j < c1; ++j)
            {
                if (j != col)
                {
                    out[sub_i][sub_j] = matrix[i][j];
                    ++sub_j;
                }
            }
            ++sub_i;
        }
    }

    return 1;
}

unsigned char getMatrixDeterminant(size_t rows, size_t cols, matrix_t matrix[rows][cols], det_t *out)
{

    if (!matrix || cols == 0 || rows == 0 || rows != cols)
    {
        return 0;
    }

    if (rows == 2 && cols == 2)
    {
        // 2x2 matrix -> (a*d)-(b*c)
        *out = (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
        return 1;
    }
    else if (rows == 3 && cols == 3)
    {
        // 3x3 Sarrus rule -> (a*e*i)+(b*f*g)+(d*h*c)-(b*d*i)-(f*h*a)-(c*e*g)
        *out = (matrix[0][0] * matrix[1][1] * matrix[2][2]) +
               (matrix[0][1] * matrix[1][2] * matrix[2][0]) +
               (matrix[1][0] * matrix[2][1] * matrix[0][2]) -
               (matrix[1][0] * matrix[0][1] * matrix[2][2]) -
               (matrix[1][2] * matrix[2][1] * matrix[0][0]) -
               (matrix[0][2] * matrix[1][1] * matrix[2][0]);
        return 1;
    }
    else
    {
        // Use Laplace's expansion for larger matrices
        det_t determinant = 0, subDet = 0;
        for (unsigned long int j = 0; j < cols; ++j)
        {
            matrix_t tmp[rows - 1][cols - 1];
            if (!createSubmatrix(rows, cols, matrix, rows - 1, cols - 1, tmp, 0, j))
            {
                return 0;
            }

            if (!getMatrixDeterminant(rows - 1, cols - 1, tmp, &subDet))
            {
                return 0;
            }

            determinant += (j % 2ULL == 0) ? matrix[0][j] * subDet : -(matrix[0][j]) * subDet;
        }
        *out = determinant;
        return 1;
    }
}

unsigned char matrixMultiply(size_t r1, size_t c1, matrix_t m1[r1][c1], size_t r2, size_t c2, matrix_t m2[r2][c2], matrix_t out[r1][c2])
{
    if (c1 != r2 || !out)
    {
        return 0;
    }

    for (unsigned long int i = 0; i < r1; i++)
    {
        for (unsigned long int j = 0; j < c2; j++)
        {
            matrix_t sum = 0;
            for (unsigned long int k = 0; k < c1; k++)
            {
                sum += m1[i][k]*m2[k][j];
            }
            out[i][j] = sum;
        }
    }
    return 1;
}

unsigned char vectorMultiply(size_t s1, matrix_t v1[s1], size_t s2, matrix_t v2[s2], det_t *out)
{
    if (!s1 || s1 != s2)
    {
        return 0;
    }

    for (size_t i = 0; i < s1; i++)
    {
        *out += (v1[i] * v2[i]);
    }
    return 1;
};
