#ifndef GEOM_H
#define GEOM_H

__host__ __device__ double get_dv_dr(int, double, double);
__host__ __device__ double get_volume(int, double, double);
void linspace(double*, double, double, int);
void logspace(double*, double, double, int);

#endif
