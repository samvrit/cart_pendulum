#include <math.h>
#include "observer_controller.h"
#include "matrix_operations.h"

#define PI (3.1415f)
#define TWO_PI (2.0f * PI)
#define DEG_TO_RAD(x)	((x) * (PI / 180.0f))

#define LPF_A_FROM_TIME_CONSTANT(Fs, Tau)   ( 1.0f / ( 1.0f + ((Tau) * (Fs)) ) )

const float A[4][4] = {	{0,    1.0000,         0,         0},
						{0,         0,         0,         0},
						{0,         0,         0,    1.0000},
						{0,         0,    1.9620,         0}};

const float B[4][4] = {	{0.0f, 0.0f, 0.0f, 0.0f}, 
						{1.0f, 0.0f, 0.0f, 0.0f},
						{0.0f, 0.0f, 0.0f, 0.0f}, 
						{0.2f, 0.0f, 0.0f, 0.0f}};

const float C[4][4] = {	{1.0f, 0.0f, 0.0f, 0.0f},
						{0.0f, 0.0f, 0.0f, 0.0f},
						{0.0f, 0.0f, 1.0f, 0.0f},
						{0.0f, 0.0f, 0.0f, 0.0f}};

const float I[4][4] = {	{1.0f, 0.0f, 0.0f, 0.0f},
						{0.0f, 1.0f, 0.0f, 0.0f},
						{0.0f, 0.0f, 1.0f, 0.0f},
						{0.0f, 0.0f, 0.0f, 1.0f}};

const float K[4][4] = {	{-1.0000, -5.2796, 97.4833, 73.0332},
						{ 0.0000,  0.0000,  0.0000,  0.0000},
						{ 0.0000,  0.0000,  0.0000,  0.0000},
						{ 0.0000,  0.0000,  0.0000,  0.0000}};

static float x_hat[N_STATES] = {0.0f};
float L[N_STATES][N_STATES] = {{0.0f}};
float F[N_STATES][N_STATES] = {{0.0f}};
float F_transpose[N_STATES][N_STATES] = {{0.0f}};
float Q_matrix[N_STATES][N_STATES] = {{0.0f}};
float R_matrix[N_STATES][N_STATES] = {{0.0f}};
float B_transpose[N_STATES][N_STATES] = {{0.0f}};
float C_transpose[N_STATES][N_STATES] = {{0.0f}};
float P[N_STATES][N_STATES] = {{0.0f}};
float A_minus_BK[N_STATES][N_STATES] = {{0.0f}};

const float Q = 1000.0f;
const float R = 1.0f;

void observer_init(const float timestep)
{
	vector_initialize(x_hat, 0.0f);
	
	matrix_initialize(L, 0.0f);
	matrix_initialize(F, 0.0f);
	matrix_initialize(F_transpose, 0.0f);
	matrix_initialize(Q_matrix, 0.0f);
	matrix_initialize(R_matrix, 0.0f);
	matrix_initialize(B_transpose, 0.0f);
	matrix_initialize(C_transpose, 0.0f);
	matrix_initialize(P, 0.0f);
	matrix_initialize(A_minus_BK, 0.0f);
	
	// Discretize state transition matrix, ie., F = I + dt*A
	matrix_scale(A, timestep, F);
	matrix_sum((const float (*)[N_STATES])F, I, F);
	
	matrix_transpose((const float (*)[N_STATES])F, F_transpose);
	matrix_transpose(B, B_transpose);
	matrix_transpose(C, C_transpose);
	
	float B_B_transpose[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(B, (const float (*)[N_STATES])B_transpose, B_B_transpose);
	
	matrix_scale((const float (*)[N_STATES])B_B_transpose, Q*Q, Q_matrix);
	matrix_scale(I, R, R_matrix);
	
	float B_K[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(B, K, B_K);
	
	matrix_diff(A, (const float (*)[N_STATES])B_K, A_minus_BK);
	matrix_scale((const float (*)[N_STATES])A_minus_BK, timestep, A_minus_BK);	
}

void observer_step(const float measurement[N_STATES], const float timestep, const bool enable, float x_hat_output[N_STATES])
{
	// P = (F * P * F') + Q
	float F_P[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply((const float (*)[N_STATES])F, (const float (*)[N_STATES])P, F_P);
	
	float F_P_F_transpose[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply((const float (*)[N_STATES])F_P, (const float (*)[N_STATES])F_transpose, F_P_F_transpose);
	
	matrix_sum((const float (*)[N_STATES])F_P_F_transpose, (const float (*)[N_STATES])Q_matrix, P);
	
	// S = (C * P * C') + R
	float C_P[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(C, (const float (*)[N_STATES])P, C_P);
	
	float C_P_C_transpose[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply((const float (*)[N_STATES])C_P, (const float (*)[N_STATES])C_transpose, C_P_C_transpose);
	
	float S[N_STATES][N_STATES] = {{0.0f}};
	matrix_sum((const float (*)[N_STATES])C_P_C_transpose, (const float (*)[N_STATES])R_matrix, S);
	
	// L = P * C' * S_inverse
	float S_inverse[N_STATES][N_STATES] = {{0.0f}};
	matrix_inverse((const float (*)[N_STATES])S, S_inverse);
	
	float P_C_transpose[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply((const float (*)[N_STATES])P, (const float (*)[N_STATES])C_transpose, P_C_transpose);
	
	matrix_matrix_multiply((const float (*)[N_STATES])P_C_transpose, (const float (*)[N_STATES])S_inverse, L);
	
	// P_next = (I - (L * C)) * P;
	float L_C[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply((const float (*)[N_STATES])L, C, L_C);
	
	float I_minus_L_C[N_STATES][N_STATES] = {{0.0f}};
	matrix_diff(I, (const float (*)[N_STATES])L_C, I_minus_L_C);
	
	float P_next[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply((const float (*)[N_STATES])I_minus_L_C, (const float (*)[N_STATES])P, P_next);
	
	matrix_assign((const float (*)[N_STATES])P_next, P);
	
	// x_hat_dot = (A-B*K)*x_hat + L*(y - C*x_hat)
	
	// error = y - C * x_hat_prev
	float C_times_x_hat[N_STATES] = {0.0f};
	matrix_vector_multiply(C, (const float *)x_hat, C_times_x_hat);

	float error[N_STATES] = {0.0f};
	vector_diff((const float *)measurement, (const float *)C_times_x_hat, error);
	
	// correction = L * error
	float correction[N_STATES] = {0.0f};
	matrix_vector_multiply((const float (*)[N_STATES])L, (const float *)error, correction);
	
	// prediction = (A-B*K) * x_hat
	float prediction[N_STATES] = {0.0f};
	matrix_vector_multiply((const float (*)[N_STATES])A_minus_BK, (const float *)x_hat, prediction);

	float prediction_plus_correction[N_STATES] = {0.0f};
	vector_sum((const float *)prediction, (const float *)correction, prediction_plus_correction);
	
	for (int i = 0; i < N_STATES; i++)
	{
		x_hat[i] += prediction_plus_correction[i];
		x_hat[i] = enable ? x_hat[i] : 0.0f;
		x_hat_output[i] = x_hat[i];
	}
}

float control_output(const float x_hat[N_STATES], const float timestep)
{	
	static float control_output_lpf = 0.0f;
	
	const bool linearity = (fabs(x_hat[2]) < DEG_TO_RAD(20.0f));
	
	const float control_output = -1.0f * dot_product(K[0], x_hat);
	
	control_output_lpf += (control_output - control_output_lpf) * LPF_A_FROM_TIME_CONSTANT(1.0f / timestep, 0.02f);
	
	return linearity ? control_output_lpf : 0.0f;
}
