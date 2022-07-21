#ifndef GEOM_H
#define GEOM_H

double get_dv_dr(int, double, double);
double get_volume(int, double, double);
void linspace(double*, double, double, int);
void logspace(double*, double, double, int);
void nonuspace(double*, double, double, int);
void nonuspace_mix(double*, double, double, int);
void nonuspace_half(double*, double, double, int);

#endif
