#ifndef SOLVER_H
#define SOLVER_H

void solve(Grid*, double, double);
void viscosity_tensor_evaluation(Grid*);
void apply_viscosity(Grid*,double);

#endif
