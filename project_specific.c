#include <math.h>
#include <stdbool.h>
#include "project_specific.h"

#define PI (3.1415f)
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

const float K[4][4] = {	{-1.0000, -5.2796, 97.4833, 73.0332},
						{ 0.0000,  0.0000,  0.0000,  0.0000},
						{ 0.0000,  0.0000,  0.0000,  0.0000},
						{ 0.0000,  0.0000,  0.0000,  0.0000}};

const float Q = 1000.0f;
const float R = 1.0f;

float control_output_process(const float computed_output, const float x_hat[N_STATES], const float timestep)
{
	static float output_lpf = 0.0f;
	
	output_lpf += (computed_output - output_lpf) * LPF_A_FROM_TIME_CONSTANT(1.0f / timestep, 0.1f);
	
	const bool linearity = (fabs(x_hat[2]) < DEG_TO_RAD(10.0f));
	
	return linearity ? output_lpf : 0.0f;	
}