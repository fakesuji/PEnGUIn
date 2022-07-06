#ifndef INIT_H
#define INIT_H

void init(Grid*);

__host__ __device__ double get_r(double,double,double);
__host__ __device__ double get_p(double,double,double);
__host__ __device__ double get_u(double,double,double);
__host__ __device__ double get_v(double,double,double);
__host__ __device__ double get_w(double,double,double);

__host__ __device__ double get_cs2(double,double,double);
__host__ __device__ double get_nu(double,double,double);
__host__ __device__ double get_h(double,double,double);

__host__ __device__ Cell init_C(double, double, double, double, double, double);
__host__ __device__ Cell init_C(double, double, double);
__host__ __device__ Dust init_CD(double, double, double);

__global__ void init(Grid, Cell*);
#endif
