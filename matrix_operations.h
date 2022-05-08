#ifndef MATRIX_OPERATIONS_H_
#define MATRIX_OPERATIONS_H_

#include <stdbool.h>
#include "project_specific.h"

void matrix_vector_multiply(const float A[N_STATES][N_STATES], const float x[N_STATES], float output[N_STATES]);

void matrix_matrix_multiply(const float A[N_STATES][N_STATES], const float B[N_STATES][N_STATES], float output[N_STATES][N_STATES]);

float dot_product(const float K[N_STATES], const float x[N_STATES]);

void vector_sum(const float x1[N_STATES], const float x2[N_STATES], float output[N_STATES]);

void matrix_sum(const float A[N_STATES][N_STATES], const float B[N_STATES][N_STATES], float output[N_STATES][N_STATES]);

void vector_diff(const float x1[N_STATES], const float x2[N_STATES], float output[N_STATES]);

void matrix_diff(const float A[N_STATES][N_STATES], const float B[N_STATES][N_STATES], float output[N_STATES][N_STATES]);

void vector_scale(const float x[N_STATES], const float a, float output[N_STATES]);

void matrix_scale(const float A[N_STATES][N_STATES], const float a, float output[N_STATES][N_STATES]);

void matrix_transpose(const float A[N_STATES][N_STATES], float output[N_STATES][N_STATES]);

void matrix_assign(const float input[N_STATES][N_STATES], float output[N_STATES][N_STATES]);

void vector_initialize(float vector[N_STATES], const float init_val);

void matrix_initialize(float matrix[N_STATES][N_STATES], const float init_val);

void identity(float matrix[N_STATES][N_STATES]);

bool matrix_inverse(const float A[N_STATES][N_STATES], float inverse[N_STATES][N_STATES]);

#endif // MATRIX_OPERATIONS_H_
