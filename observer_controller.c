#include <math.h>
#include "observer_controller.h"
#include "matrix_operations.h"

#define PI (3.1415f)
#define TWO_PI (2.0f * PI)
#define TEN_DEG_TO_RAD (0.174533f)

float A[4][4] = {	{0,    1.0000,         0,         0},
					{0,         0,         0,         0},
					{0,         0,         0,    1.0000},
					{0,         0,    1.9620,         0}};

float B[4][4] = {	{0.0f, 0.0f, 0.0f, 0.0f}, 
					{1.0f, 0.0f, 0.0f, 0.0f},
					{0.0f, 0.0f, 0.0f, 0.0f}, 
					{0.2f, 0.0f, 0.0f, 0.0f}};

float C[4][4] = {	{1.0f, 0.0f, 0.0f, 0.0f},
					{0.0f, 0.0f, 0.0f, 0.0f},
					{0.0f, 0.0f, 1.0f, 0.0f},
					{0.0f, 0.0f, 0.0f, 0.0f}};

float I[4][4] = {	{1.0f, 0.0f, 0.0f, 0.0f},
					{0.0f, 1.0f, 0.0f, 0.0f},
					{0.0f, 0.0f, 1.0f, 0.0f},
					{0.0f, 0.0f, 0.0f, 1.0f}};

float A_minus_BK[4][4] = {	{0,         1.0000,         0,         0},
							{1.0000,    5.2796,  -97.4833,  -73.0332},
							{0,         0,              0,    1.0000},
							{0.2000,    1.0559,  -17.5347,  -14.6066}};

float K[4] = {-1.0000, -5.2796, 97.4833, 73.0332};

static float x_hat[N_STATES] = {0.0f};
float L[N_STATES][N_STATES] = {{0.0f}};
float F[N_STATES][N_STATES] = {{0.0f}};
float F_transpose[N_STATES][N_STATES] = {{0.0f}};
float Q_matrix[N_STATES][N_STATES] = {{0.0f}};
float R_matrix[N_STATES][N_STATES] = {{0.0f}};
float B_transpose[N_STATES][N_STATES] = {{0.0f}};
float C_transpose[N_STATES][N_STATES] = {{0.0f}};
float P[N_STATES][N_STATES] = {{0.0f}};

float Q = 1000.0f;
float R = 1e-6f;

void observer_init(float timestep)
{
	for (int i = 0; i < N_STATES; i++)
	{
		x_hat[i] = 0.0f;
	}
	
	matrix_scale(A_minus_BK, timestep, A_minus_BK);
	
	matrix_scale(A, timestep, F);
	matrix_sum(F, I, F);
	
	matrix_transpose(F, F_transpose);
	matrix_transpose(B, B_transpose);
	matrix_transpose(C, C_transpose);
	
	float B_B_transpose[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(B, B_transpose, B_B_transpose);
	
	matrix_scale(B_B_transpose, Q, Q_matrix);
	matrix_scale(I, R, R_matrix);
}

void observer_step(float measurement[N_STATES], float timestep, bool enable, float x_hat_output[N_STATES])
{
	// P = (F * P * F') + Q
	float F_P[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(F, P, F_P);
	
	float F_P_F_transpose[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(F_P, F_transpose, F_P_F_transpose);
	
	matrix_sum(F_P_F_transpose, Q_matrix, P);
	
	// S = (C * P * C') + R
	float C_P[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(C, P, C_P);
	
	float C_P_C_transpose[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(C_P, C_transpose, C_P_C_transpose);
	
	float S[N_STATES][N_STATES] = {{0.0f}};
	matrix_sum(C_P_C_transpose, R_matrix, S);
	
	// L = P * C' * S_inverse
	float S_inverse[N_STATES][N_STATES] = {{0.0f}};
	matrix_inverse(S, S_inverse);
	
	float P_C_transpose[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(P, C_transpose, P_C_transpose);
	
	matrix_matrix_multiply(P_C_transpose, S_inverse, L);
	
	// P_next = (I - (L * C)) * P;
	float L_C[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(L, C, L_C);
	
	float I_minus_L_C[N_STATES][N_STATES] = {{0.0f}};
	matrix_diff(I, L_C, I_minus_L_C);
	
	float P_next[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(I_minus_L_C, P, P_next);
	
	matrix_assign(P_next, P);
	
	// x_hat_dot = (A-B*K)*x_hat + L*(y - C*x_hat)
	
	// error = y - C * x_hat_prev
	float C_times_x_hat[N_STATES] = {0.0f};
	matrix_vector_multiply(C, x_hat, C_times_x_hat);

	float error[N_STATES] = {0.0f};
	vector_diff(measurement, C_times_x_hat, error);
	
	// correction = L * error
	float correction[N_STATES] = {0.0f};
	matrix_vector_multiply(L, error, correction);
	
	// prediction = (A-B*K) * x_hat
	float prediction[N_STATES] = {0.0f};
	matrix_vector_multiply(A_minus_BK, x_hat, prediction);

	float prediction_plus_correction[N_STATES] = {0.0f};
	vector_sum(prediction, correction, prediction_plus_correction);

	float angle_wrapped = (x_hat[2] > PI) ? (x_hat[2] - TWO_PI) : x_hat[2];
	
	bool linearity = (fabs(angle_wrapped) < TEN_DEG_TO_RAD);
	
	for (int i = 0; i < N_STATES; i++)
	{
		x_hat[i] += (enable && linearity) ? prediction_plus_correction[i] : 0.0f;
		x_hat_output[i] = x_hat[i];
	}
}

float control_output(float x_hat[N_STATES])
{	
	return ( -1.0f * dot_product(K, x_hat) );
}
