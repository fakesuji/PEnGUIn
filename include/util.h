#ifndef UTIL_H
#define UTIL_H

void wait_f_r();
string int_to_string (int);
string frame_num (int);
string path_to_cwd();

//######################## GPU CONTROL ############################

void P2P_all_enable(int,int);
void P2P_all_disable(int,int);
void syncdevices(int,int);
void syncdestreams(int,Grid*);

//######################## FILE CONTROL ############################

void open_output_file(ofstream&,string);
void append_output_file(ofstream&,string);
void close_output_file(ofstream&);
void open_binary_file(ofstream&,string);
void open_binary_file(ifstream&,string);

//######################## REDUCTION ############################

__device__ void shift_round_reduc_max(int,int,double*);
__device__ void round_reduc_max(int,double*);
__device__ void round_reduc_sum(int,double*);
__device__ void bin_reduc_max(int,double*);
__device__ void bin_reduc_sum(int,double*);

//######################## MISC ############################

__device__ double min3(double,double,double);
__device__ double max3(double,double,double);
__device__ int periodic_ind(int,int);

#endif
