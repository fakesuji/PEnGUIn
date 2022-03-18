#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "structs.h"

using namespace std;

void wait_f_r()
{
  char c;
  cout << " To go on press return! " << endl;
  cin.get(c);
  return;
}

string int_to_string (int num)
{
  stringstream A;
  A << num;
  return (A.str());
}

string frame_num (int num)
{
  string B;
  string C = int_to_string(num);

  if(C.length()==5) B=C;
  else if(C.length()==4) B="0"+C;
  else if(C.length()==3) B="00"+C;
  else if(C.length()==2) B="000"+C;
  else if(C.length()==1) B="0000"+C;
  return (B);
}

string path_to_cwd()
{
  char* a_cwd = getcwd(NULL, 0);
  string s_cwd(a_cwd);
  free(a_cwd);
  return s_cwd;
}

//######################## GPU CONTROL ############################

void P2P_all_enable(int start, int nDev)
{
	for (int i=start; i<nDev+start; i++)
		for (int j=i+1; j<nDev+start; j++)
		{
			cudaSetDevice(i);
			cudaDeviceEnablePeerAccess(j, 0);
			cudaSetDevice(j);
			cudaDeviceEnablePeerAccess(i, 0);
    		}
	return;
}

void P2P_all_disable(int start, int nDev)
{
	for (int i=start; i<nDev+start; i++)
		for (int j=i+1; j<nDev+start; j++)
		{
			cudaSetDevice(i);
			cudaDeviceDisablePeerAccess(j);
			cudaSetDevice(j);
			cudaDeviceDisablePeerAccess(i);
    		}
	return;
}


void syncdevices(int start, int nDev)
{
	for (int i=start; i<nDev+start; i++)
	{
		cudaSetDevice(i);
		cudaDeviceSynchronize();
	}
	return;
}

void syncstreams(int ndev, Grid* dev)
{
	for (int n=0; n<ndev; n++)
		cudaStreamSynchronize(dev[n].stream);
}

//######################## FILE CONTROL ############################

void open_output_file(ofstream &wfile, string sfname)
{
	wfile.open(sfname.c_str(), ios::out);
	wfile.precision(16);

	return;
}

void append_output_file(ofstream &wfile, string sfname)
{
	wfile.open(sfname.c_str(), ios::app);
	wfile.precision(16);
	return;
}

void close_output_file(ofstream& wfile)
{
	wfile.close();
	return;
}

void close_output_file(ifstream& wfile)
{
	wfile.close();
	return;
}

void open_binary_file(ofstream &wfile, string sfname)
{
	wfile.open(sfname.c_str(), ios::out | ios::binary);
	return;
}

void open_binary_file(ifstream &wfile, string sfname)
{
	wfile.open(sfname.c_str(), ios::in | ios::binary);
	return;
}

//######################## REDUCTION ############################

__device__ int two_round(int x)
{
  int y = 1;
  while (y<x)
    y *= 2;
  return y;
}

__device__ int log2(int x)
{
  int y = 0;
  while (x>1)
  {
    x /= 2;
    y++;
  }
  return y;
}

__device__ void shift_round_reduc_max(int nmin, int len, double *list)
{

  int halfPoint, thread2;
  int len2 = two_round(len); // Total number of threads, rounded up to the next power of two
  int n = threadIdx.x;

  while(len2 > 1)
  {
    halfPoint = len2/2;
    // only the first half of the threads will be active.
 
    if (n < halfPoint)
    {
      thread2 = n + halfPoint;
 
      // Skipping the fictious threads blockDim.x ... blockDim_2-1
      if (thread2 < len)
      {
        if (list[n+nmin]<list[thread2+nmin]) list[n+nmin] = list[thread2+nmin];
      }
    }
    __syncthreads();
 
    // Reducing the binary tree size by two:
    len2 = halfPoint;
  }
  return;
}

__device__ void round_reduc_max(int len, double *list)
{
  int halfPoint, thread2;
  int len2 = two_round(len);
  int n = threadIdx.x + blockDim.x*threadIdx.y;

  while(len2 > 1)
  {
    halfPoint = len2/2;
 
    if (n < halfPoint)
    {
      thread2 = n + halfPoint;
      if (thread2 < len)
      {
        if (list[n]<list[thread2]) list[n] = list[thread2];
      }
    }
    __syncthreads();
 
    len2 = halfPoint;
  }
  return;
}

__device__ void round_reduc_sum(int len, double *list)
{
  int halfPoint, thread2;
  int len2 = two_round(len);
  int n = threadIdx.x + blockDim.x*threadIdx.y;

  while(len2 > 1)
  {
    halfPoint = len2/2;
 
    if (n < halfPoint)
    {
      thread2 = n + halfPoint;
      if (thread2 < len)
      {
        list[n] += list[thread2];
      }
    }
    __syncthreads();
 
    len2 = halfPoint;
  }
  return;
}

__device__ void bin_reduc_max(int len, double *list)
{
  int halfPoint, thread2;
  int n = threadIdx.x + blockDim.x*threadIdx.y;

  while(len > 1)
  {
    halfPoint = len/2;
 
    if (n < halfPoint)
    {
      thread2 = n + halfPoint;
      if (list[n]<list[thread2]) list[n] = list[thread2];
    }
    __syncthreads();
 
    len = halfPoint;
  }
  return;
}

__device__ void bin_reduc_sum(int len, double *list)
{
  int halfPoint, thread2;
  int n = threadIdx.x + blockDim.x*threadIdx.y;

  while(len > 1)
  {
    halfPoint = len/2;
 
    if (n < halfPoint)
    {
      thread2 = n + halfPoint;
      list[n] += list[thread2];
    }
    __syncthreads();
 
    len = halfPoint;
  }
  return;
}

//######################## MISC ############################

__device__ double min3(double A, double B, double C)
{
  return fmin(A, fmin(B, C));
}

__device__ double max3(double A, double B, double C)
{
  return fmax(A, fmax(B, C));
}

__device__ int periodic_ind(int j, int jmax)
{
  if      (j>=jmax) return j%jmax;
  else if (j< 0)    return (j%jmax) + jmax;
  else              return j;
}

