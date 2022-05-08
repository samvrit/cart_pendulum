#ifndef OBSERVER_CONTROLLER_H_
#define OBSERVER_CONTROLLER_H_

#include <stdbool.h>
#include "project_specific.h"

extern const float A[N_STATES][N_STATES];
extern const float B[N_STATES][N_STATES];
extern const float C[N_STATES][N_STATES];
extern const float I[N_STATES][N_STATES];
extern const float K[N_STATES][N_STATES];

extern const float Q;
extern const float R;

void observer_init(const float timestep);

void observer_step(const float measurement[N_STATES], const float timestep, const bool enable, float x_hat_output[N_STATES]);

float control_output(const float x_hat[N_STATES], const float timestep);

#endif // OBSERVER_CONTROLLER_H_
