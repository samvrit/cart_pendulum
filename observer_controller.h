#ifndef OBSERVER_CONTROLLER_H_
#define OBSERVER_CONTROLLER_H_

#include <stdbool.h>

#define N_STATES (4)

void observer_init(const float timestep);

void observer_step(const float measurement[N_STATES], const float timestep, const bool enable, float x_hat_output[N_STATES]);

float control_output(const float x_hat[N_STATES]);

#endif // OBSERVER_CONTROLLER_H_
