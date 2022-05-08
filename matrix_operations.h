#ifndef MATRIX_OPERATIONS_H_
#define MATRIX_OPERATIONS_H_

#include <stdbool.h>

#define N_STATES (4)

void matrix_vector_multiply(float A[N_STATES][N_STATES], float x[N_STATES], float output[N_STATES]);

void matrix_matrix_multiply(float A[N_STATES][N_STATES], float B[N_STATES][N_STATES], float output[N_STATES][N_STATES]);

float dot_product(float K[N_STATES], float x[N_STATES]);

void vector_sum(float x1[N_STATES], float x2[N_STATES], float output[N_STATES]);

void matrix_sum(float A[N_STATES][N_STATES], float B[N_STATES][N_STATES], float output[N_STATES][N_STATES]);

void vector_diff(float x1[N_STATES], float x2[N_STATES], float output[N_STATES]);

void matrix_diff(float A[N_STATES][N_STATES], float B[N_STATES][N_STATES], float output[N_STATES][N_STATES]);

void vector_scale(float x[N_STATES], float a, float output[N_STATES]);

void matrix_scale(float A[N_STATES][N_STATES], float a, float output[N_STATES][N_STATES]);

void matrix_transpose(float A[N_STATES][N_STATES], float output[N_STATES][N_STATES]);

void matrix_assign(float input[N_STATES][N_STATES], float output[N_STATES][N_STATES]);

bool matrix_inverse(float A[N_STATES][N_STATES], float inverse[N_STATES][N_STATES]);

#endif // MATRIX_OPERATIONS_H_
