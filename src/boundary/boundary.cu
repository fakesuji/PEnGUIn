#include <unistd.h>
#include <iostream>

#include "parameters.h"
#include "structs.h"
#include "init.h"

void syncallstreams(Grid* dev)
{
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
}

__global__ void bound_x_left(Grid G, Cell* X)
{
	int n = threadIdx.x;
	int j = blockIdx.x;
	int k = blockIdx.y;

	Cell* C = &X[G.xarr*(j+G.yarr*k)];

	if (bound_lft == 2)
	{
		C[n].copy(C[2*xpad-1-n]);
		C[n].u *= -1.0;
	}
	else if (bound_lft == 1)
	{
		C[n].copy(C[xpad]);
	}
	else if (bound_lft == 0)
	{
		C[n].copy(init_C(G.get_xc(n),G.get_yc(j),G.get_zc(k)));
	}

	return;
}

__global__ void bound_x_rght(Grid G, Cell* X)
{
	int nx = G.xres+xpad;
	int n = threadIdx.x;
	int j = blockIdx.x;
	int k = blockIdx.y;

	Cell* C = &X[G.xarr*(j+G.yarr*k)];

	if (bound_rgh == 2)
	{
		C[nx+n].copy(C[nx-1-n]);
		C[nx+n].u *= -1.0;
	}
	else if (bound_rgh == 1)
	{
		C[nx+n].copy(C[nx-1]);
	}
	else if (bound_rgh == 0)
	{
		C[nx+n].copy(init_C(G.get_xc(nx+n),G.get_yc(j),G.get_zc(k)));
	}

	return;
}

__global__ void bound_x_left(Grid G, Dust* X)
{
	int n = threadIdx.x;
	int j = blockIdx.x;
	int k = blockIdx.y;

	Dust* D = &X[G.xarr*(j+G.yarr*k)];

	if (bound_lft == 2)
	{
		D[n].copy(D[2*xpad-1-n]);
		D[n].u *= -1.0;
	}
	else if (bound_lft == 1)
	{
		D[n].copy(D[xpad]);
	}
	else if (bound_lft == 0)
	{
		D[n].copy(init_CD(G.get_xc(n),G.get_yc(j),G.get_zc(k)));
	}

	return;
}

__global__ void bound_x_rght(Grid G, Dust* X)
{
	int nx = G.xres+xpad;
	int n = threadIdx.x;
	int j = blockIdx.x;
	int k = blockIdx.y;

	Dust* D = &X[G.xarr*(j+G.yarr*k)];

	if (bound_rgh == 2)
	{
		D[nx+n].copy(D[nx-1-n]);
		D[nx+n].u *= -1.0;
	}
	else if (bound_rgh == 1)
	{
		D[nx+n].copy(D[nx-1]);
	}
	else if (bound_rgh == 0)
	{
		D[nx+n].copy(init_CD(G.get_xc(nx+n),G.get_yc(j),G.get_zc(k)));
	}

	return;
}

__device__ Cell left_state()
{
	Cell Q;

	Q.r = 8.0;
	Q.p = 116.5;
	Q.u = 7.14470958122;
	Q.v =-4.125;
	Q.w = 0.0;

	return Q;
}

__device__ Cell right_state()
{
	Cell Q;

	Q.r = 1.4;
	Q.p = 1.0;
	Q.u = 0.0;
	Q.v = 0.0;
	Q.w = 0.0;

	return Q;
}

__global__ void boundy(Grid G, Cell* C, double t)
{
	int n = threadIdx.x;
	int i = blockIdx.x;
	int k = blockIdx.y + zpad;
	int ib = blockIdx.z;
	Cell T;
	double x,y;

	if (ib==0)
	{
		if (bound_bak == 3)
		{
			T.copy(C[G.get_ind(i,G.yres+n,k)]);
		}
		else if (bound_bak == 2)
		{
			T.copy(C[G.get_ind(i,2*ypad-1-n,k)]);
			#if init_flag == 10
			x = G.get_xc(i);
			if (x>1.0/6.0) T.v *= -1.0;
			else T.copy(left_state());
			#else
			T.v *= -1.0;
			#endif
			
		}
		else if (bound_bak == 1)
		{
			T.copy(C[G.get_ind(i,ypad,k)]);
		}
		else if (bound_bak == 0)
		{
			x = G.get_xc(i);
			y = G.get_yc(n);
			T = init_C(x,y,G.get_zc(k));
		}
		C[G.get_ind(i,n,k)].copy(T);
	}
	else if (ib==1)
	{
		if (bound_frn == 3)
		{
			T.copy(C[G.get_ind(i,ypad+n,k)]);
		}
		else if (bound_frn == 2)
		{
			T.copy(C[G.get_ind(i,G.yres+ypad-1-n,k)]);
			T.v *= -1.0;
		}
		else if (bound_frn == 1)
		{
			T.copy(C[G.get_ind(i,G.yres+ypad-1,k)]);
		}
		else if (bound_frn == 0)
		{
			x = G.get_xc(i);
			y = G.get_yc(G.yres+ypad+n);
			#if init_flag == 10
			if (x>1.0/6.0 + 0.57735026919*y + t*7.14470958122) T.copy(right_state());
			else                                               T.copy(left_state());
			#else
			T = init_C(x,y,G.get_zc(k));
			#endif
		}
		C[G.get_ind(i,G.yres+ypad+n,k)].copy(T);
	}
	return;
}

__global__ void boundy_dust(Grid G, Dust* C, double t)
{
	int n = threadIdx.x;
	int i = blockIdx.x;
	int k = blockIdx.y + zpad;
	int ib = blockIdx.z;
	Dust T;
	double x,y;

	if (ib==0)
	{
		if (bound_bak == 3)
		{
			T.copy(C[G.get_ind(i,G.yres+n,k)]);
		}
		else if (bound_bak == 2)
		{
			T.copy(C[G.get_ind(i,2*ypad-1-n,k)]);
			T.v *= -1.0;
		}
		else if (bound_bak == 1)
		{
			T.copy(C[G.get_ind(i,ypad,k)]);
		}
		else if (bound_bak == 0)
		{
			x = G.get_xc(i);
			y = G.get_yc(n);
			T = init_CD(x,y,G.get_zc(k));
		}
		C[G.get_ind(i,n,k)].copy(T);
		//printf("%i %e %e, %e\n", i, C[G.get_ind(i,G.yres+n,k)].r,C[G.get_ind(i,n,k)].r, C[G.get_ind(i,G.yres+n,k)].r-C[G.get_ind(i,n,k)].r);
	}
	else if (ib==1)
	{
		if (bound_frn == 3)
		{
			T.copy(C[G.get_ind(i,ypad+n,k)]);
		}
		else if (bound_frn == 2)
		{
			T.copy(C[G.get_ind(i,G.yres+ypad-1-n,k)]);
			T.v *= -1.0;
		}
		else if (bound_frn == 1)
		{
			T.copy(C[G.get_ind(i,G.yres+ypad-1,k)]);
		}
		else if (bound_frn == 0)
		{
			x = G.get_xc(i);
			y = G.get_yc(G.yres+ypad+n);
			T = init_CD(x,y,G.get_zc(k));
		}
		C[G.get_ind(i,G.yres+ypad+n,k)].copy(T);
	}
	return;
}

__global__ void boundz(Grid G, Cell* C)
{
	int n = threadIdx.x;
	int i = blockIdx.x;
	int j = blockIdx.y;
	int ib = blockIdx.z;
	Cell T;

	if (ib==0)
	{
		if (bound_bom == 3)
		{
			T.copy(C[G.get_ind(i,j,G.zres+n)]);
		}
		else if (bound_bom == 2)
		{
			T.copy(C[G.get_ind(i,j,2*zpad-1-n)]);
			T.w *= -1.0;
		}
		else if (bound_bom == 1)
		{
			T.copy(C[G.get_ind(i,j,zpad)]);
		}
		else if (bound_bom == 0)
		{
			T = init_C(G.get_xc(i),G.get_yc(j),G.get_zc(n));
		}
		C[G.get_ind(i,j,n)].copy(T);
	}
	else if (ib==1)
	{
		if (bound_top == 3)
		{
			T.copy(C[G.get_ind(i,j,zpad+n)]);
		}
		else if (bound_top == 2)
		{
			T.copy(C[G.get_ind(i,j,G.zres+zpad-1-n)]);
			T.w *= -1.0;
		}
		else if (bound_top == 1)
		{
			T.copy(C[G.get_ind(i,j,G.zres+zpad-1)]);
		}
		else if (bound_top == 0)
		{
			T = init_C(G.get_xc(i),G.get_yc(j),G.get_zc(G.zres+zpad+n));
		}
		C[G.get_ind(i,j,G.zres+zpad+n)].copy(T);
	}
	return;
}

__global__ void boundz_dust(Grid G, Dust* C)
{
	int n = threadIdx.x;
	int i = blockIdx.x;
	int j = blockIdx.y;
	int ib = blockIdx.z;
	Dust T;

	if (ib==0)
	{
		if (bound_bom == 3)
		{
			T.copy(C[G.get_ind(i,j,G.zres+n)]);
		}
		else if (bound_bom == 2)
		{
			T.copy(C[G.get_ind(i,j,2*zpad-1-n)]);
			T.w *= -1.0;
		}
		else if (bound_bom == 1)
		{
			T.copy(C[G.get_ind(i,j,zpad)]);
		}
		else if (bound_bom == 0)
		{
			T = init_CD(G.get_xc(i),G.get_yc(j),G.get_zc(n));
		}
		C[G.get_ind(i,j,n)].copy(T);
	}
	else if (ib==1)
	{
		if (bound_top == 3)
		{
			T.copy(C[G.get_ind(i,j,zpad+n)]);
		}
		else if (bound_top == 2)
		{
			T.copy(C[G.get_ind(i,j,G.zres+zpad-1-n)]);
			T.w *= -1.0;
		}
		else if (bound_top == 1)
		{
			T.copy(C[G.get_ind(i,j,G.zres+zpad-1)]);
		}
		else if (bound_top == 0)
		{
			T = init_CD(G.get_xc(i),G.get_yc(j),G.get_zc(G.zres+zpad+n));
		}
		C[G.get_ind(i,j,G.zres+zpad+n)].copy(T);
	}
	return;
}

//=========================================================================================

__global__ void buff_left(int x_res, int x_arr, Cell* C, Cell* B)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int k = threadIdx.z + blockIdx.z*blockDim.z;

	B[i + xpad*(j + yarr*k)].copy(C[i+xpad + x_arr*(j + yarr*k)]);
	return;
}

__global__ void buff_rght(int x_res, int x_arr, Cell* C, Cell* B)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int k = threadIdx.z + blockIdx.z*blockDim.z;

	B[i + xpad*(j + yarr*k)].copy(C[i + x_res + x_arr*(j + yarr*k)]);
	return;
}

__global__ void dump_left(int x_res, int x_arr, Cell* C, Cell* B)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int k = threadIdx.z + blockIdx.z*blockDim.z;

	C[i + x_arr*(j + yarr*k)].copy(B[i + xpad*(j + yarr*k)]);
	return;
}

__global__ void dump_rght(int x_res, int x_arr, Cell* C, Cell* B)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int k = threadIdx.z + blockIdx.z*blockDim.z;

	C[i+x_res+xpad + x_arr*(j + yarr*k)].copy(B[i + xpad*(j + yarr*k)]);
	return;
}

//=========================================================================================

__global__ void buff_left(int x_res, int x_arr, Dust* C, Dust* B)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int k = threadIdx.z + blockIdx.z*blockDim.z;

	B[i + xpad*(j + yarr*k)].copy(C[i+xpad + x_arr*(j + yarr*k)]);
	return;
}

__global__ void buff_rght(int x_res, int x_arr, Dust* C, Dust* B)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int k = threadIdx.z + blockIdx.z*blockDim.z;

	B[i + xpad*(j + yarr*k)].copy(C[i + x_res + x_arr*(j + yarr*k)]);
	return;
}

__global__ void dump_left(int x_res, int x_arr, Dust* C, Dust* B)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int k = threadIdx.z + blockIdx.z*blockDim.z;

	C[i + x_arr*(j + yarr*k)].copy(B[i + xpad*(j + yarr*k)]);
	return;
}

__global__ void dump_rght(int x_res, int x_arr, Dust* C, Dust* B)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int k = threadIdx.z + blockIdx.z*blockDim.z;

	C[i+x_res+xpad + x_arr*(j + yarr*k)].copy(B[i + xpad*(j + yarr*k)]);
	return;
}

//=========================================================================================

void boundx(Grid* dev)
{
	int gdimy, gdimz, bdimy, bdimz;
	if (yarr%16==0)
	{
		gdimy = yarr/16;
		bdimy = 16;
	}
	else
	{
		gdimy = yarr;
		bdimy = 1;
	}
	if (zarr%16==0)
	{
		gdimz = zarr/16;
		bdimz = 16;
	}
	else
	{
		gdimz = zarr;
		bdimz = 1;
	}

	dim3 gdim(1,gdimy,gdimz);
	dim3 bdim(xpad,bdimy,bdimz);

	////////////////////////////////

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);
		buff_rght<<< gdim, bdim, 0, dev[n].stream >>> (dev[n].xres, dev[n].xarr, dev[n].C, dev[n].BuffR);
	}
	for (int n=0; n<ndev-1; n++)
	{
		cudaMemcpyAsync( dev[n+1].BuffL, dev[n].BuffR, xpad*yarr*zarr*sizeof(Cell), cudaMemcpyDeviceToDevice, dev[n].stream);
	}
		cudaMemcpyAsync( dev[0].BuffL, dev[ndev-1].BuffR, xpad*yarr*zarr*sizeof(Cell), cudaMemcpyDeviceToDevice, dev[ndev-1].stream);
	syncallstreams(dev);

	for (int n=0; n<ndev-1; n++)
	{
		cudaSetDevice(n+1);
		dump_left<<< gdim, bdim, 0, dev[n+1].stream >>> (dev[n+1].xres, dev[n+1].xarr, dev[n+1].C, dev[n+1].BuffL);
	}
	syncallstreams(dev);

	////////////////////////////////

	if (bound_lft == 3)
	{
		cudaSetDevice(0);
		dump_left<<< gdim, bdim, 0, dev[0].stream >>> (dev[0].xres, dev[0].xarr, dev[0].C, dev[0].BuffL);
	}
	else
	{
		cudaSetDevice(0);
		bound_x_left<<< dim3(dev[0].yarr,dev[0].zarr,1) , dim3(xpad,1,1) >>> (dev[0], dev[0].C);
	}
	syncallstreams(dev);

	////////////////////////////////

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);
		buff_left<<< gdim, bdim, 0, dev[n].stream >>> (dev[n].xres, dev[n].xarr, dev[n].C, dev[n].BuffL);
	}
	for (int n=1; n<ndev; n++)
	{
		cudaMemcpyAsync( dev[n-1].BuffR, dev[n].BuffL, xpad*yarr*zarr*sizeof(Cell), cudaMemcpyDeviceToDevice, dev[n].stream);
	}
		cudaMemcpyAsync( dev[ndev-1].BuffR, dev[0].BuffL, xpad*yarr*zarr*sizeof(Cell), cudaMemcpyDeviceToDevice, dev[0].stream);
	syncallstreams(dev);

	for (int n=1; n<ndev; n++)
	{
		cudaSetDevice(n-1);
		dump_rght<<< gdim, bdim, 0, dev[n-1].stream >>> (dev[n-1].xres, dev[n-1].xarr, dev[n-1].C, dev[n-1].BuffR);
	}
	syncallstreams(dev);

	////////////////////////////////

	if (bound_rgh == 3)
	{
		cudaSetDevice(ndev-1);
		dump_rght<<< gdim, bdim, 0, dev[ndev-1].stream >>> (dev[ndev-1].xres, dev[ndev-1].xarr, dev[ndev-1].C, dev[ndev-1].BuffR);
	}
	else
	{
		cudaSetDevice(ndev-1);
		bound_x_rght<<< dim3(dev[ndev-1].yarr,dev[ndev-1].zarr,1) , dim3(xpad,1,1) >>> (dev[ndev-1], dev[ndev-1].C);
	}
	syncallstreams(dev);

	return;
}

void boundy(Grid* dev, double t=0.0)
{
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		boundy<<< dim3(dev[n].xarr,dev[n].zres,2), dim3(ypad,1,1), 0, dev[n].stream >>>(dev[n], dev[n].C, t);
	}

	syncallstreams(dev);

	return;
}

void boundz(Grid* dev)
{
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		boundz<<< dim3(dev[n].xarr,dev[n].yarr,2), dim3(zpad,1,1), 0, dev[n].stream >>>(dev[n], dev[n].C);
	}

	syncallstreams(dev);

	return;
}

//=========================================================================================

void boundx_dust(Grid* dev)
{
	int gdimy, gdimz, bdimy, bdimz;
	if (yarr%16==0)
	{
		gdimy = yarr/16;
		bdimy = 16;
	}
	else
	{
		gdimy = yarr;
		bdimy = 1;
	}
	if (zarr%16==0)
	{
		gdimz = zarr/16;
		bdimz = 16;
	}
	else
	{
		gdimz = zarr;
		bdimz = 1;
	}

	dim3 gdim(1,gdimy,gdimz);
	dim3 bdim(xpad,bdimy,bdimz);

	//============================================================

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);
		buff_rght<<< gdim, bdim, 0, dev[n].stream >>> (dev[n].xres, dev[n].xarr, dev[n].CD, dev[n].BuffRD);
	}
	for (int n=0; n<ndev-1; n++)
	{
		cudaMemcpyAsync( dev[n+1].BuffLD, dev[n].BuffRD, xpad*yarr*zarr*sizeof(Dust), cudaMemcpyDeviceToDevice, dev[n].stream);
	}
		cudaMemcpyAsync( dev[0].BuffLD, dev[ndev-1].BuffRD, xpad*yarr*zarr*sizeof(Dust), cudaMemcpyDeviceToDevice, dev[ndev-1].stream);
	syncallstreams(dev);

	for (int n=0; n<ndev-1; n++)
	{
		cudaSetDevice(n+1);
		dump_left<<< gdim, bdim, 0, dev[n+1].stream >>> (dev[n+1].xres, dev[n+1].xarr, dev[n+1].CD, dev[n+1].BuffLD);
	}
	syncallstreams(dev);

	////////////////////////////////

	if (bound_lft == 3)
	{
		cudaSetDevice(0);
		dump_left<<< gdim, bdim, 0, dev[0].stream >>> (dev[0].xres, dev[0].xarr, dev[0].CD, dev[0].BuffLD);
	}
	else
	{
		cudaSetDevice(0);
		bound_x_left<<< dim3(dev[0].yarr,dev[0].zarr,1) , dim3(xpad,1,1) >>> (dev[0], dev[0].CD);
	}
	syncallstreams(dev);

	//============================================================

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);
		buff_left<<< gdim, bdim, 0, dev[n].stream >>> (dev[n].xres, dev[n].xarr, dev[n].CD, dev[n].BuffLD);
	}
	for (int n=1; n<ndev; n++)
	{
		cudaMemcpyAsync( dev[n-1].BuffRD, dev[n].BuffLD, xpad*yarr*zarr*sizeof(Dust), cudaMemcpyDeviceToDevice, dev[n].stream);
	}
		cudaMemcpyAsync( dev[ndev-1].BuffRD, dev[0].BuffLD, xpad*yarr*zarr*sizeof(Dust), cudaMemcpyDeviceToDevice, dev[0].stream);
	syncallstreams(dev);

	for (int n=1; n<ndev; n++)
	{
		cudaSetDevice(n-1);
		dump_rght<<< gdim, bdim, 0, dev[n-1].stream >>> (dev[n-1].xres, dev[n-1].xarr, dev[n-1].CD, dev[n-1].BuffRD);
	}
	syncallstreams(dev);

	////////////////////////////////

	if (bound_rgh == 3)
	{
		cudaSetDevice(ndev-1);
		dump_rght<<< gdim, bdim, 0, dev[ndev-1].stream >>> (dev[ndev-1].xres, dev[ndev-1].xarr, dev[ndev-1].CD, dev[ndev-1].BuffRD);
	}
	else
	{
		cudaSetDevice(ndev-1);
		bound_x_rght<<< dim3(dev[ndev-1].yarr,dev[ndev-1].zarr,1) , dim3(xpad,1,1) >>> (dev[ndev-1], dev[ndev-1].CD);
	}
	syncallstreams(dev);

	return;
}

void boundy_dust(Grid* dev, double t=0.0)
{
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		boundy_dust<<< dim3(dev[n].xarr,dev[n].zres,2), dim3(ypad,1,1), 0, dev[n].stream >>>(dev[n], dev[n].CD, t);
	}

	syncallstreams(dev);

	return;
}

void boundz_dust(Grid* dev)
{
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		boundz_dust<<< dim3(dev[n].xarr,dev[n].yarr,2), dim3(zpad,1,1), 0, dev[n].stream >>>(dev[n], dev[n].CD);
	}

	syncallstreams(dev);

	return;
}

//=========================================================================================

void boundx2(Grid* dev)
{
	int gdimy, gdimz, bdimy, bdimz;
	if (yarr%16==0)
	{
		gdimy = yarr/16;
		bdimy = 16;
	}
	else
	{
		gdimy = yarr;
		bdimy = 1;
	}
	if (zarr%16==0)
	{
		gdimz = zarr/16;
		bdimz = 16;
	}
	else
	{
		gdimz = zarr;
		bdimz = 1;
	}

	dim3 gdim(1,gdimy,gdimz);
	dim3 bdim(xpad,bdimy,bdimz);

	////////////////////////////////

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);
		buff_rght<<< gdim, bdim, 0, dev[n].stream >>> (dev[n].xres, dev[n].xarr, dev[n].T, dev[n].BuffR);
	}
	for (int n=0; n<ndev-1; n++)
	{
		cudaMemcpyAsync( dev[n+1].BuffL, dev[n].BuffR, xpad*yarr*zarr*sizeof(Cell), cudaMemcpyDeviceToDevice, dev[n].stream);
	}
		cudaMemcpyAsync( dev[0].BuffL, dev[ndev-1].BuffR, xpad*yarr*zarr*sizeof(Cell), cudaMemcpyDeviceToDevice, dev[ndev-1].stream);
	syncallstreams(dev);

	for (int n=0; n<ndev-1; n++)
	{
		cudaSetDevice(n+1);
		dump_left<<< gdim, bdim, 0, dev[n+1].stream >>> (dev[n+1].xres, dev[n+1].xarr, dev[n+1].T, dev[n+1].BuffL);
	}
	syncallstreams(dev);

	////////////////////////////////

	if (bound_lft == 3)
	{
		cudaSetDevice(0);
		dump_left<<< gdim, bdim, 0, dev[0].stream >>> (dev[0].xres, dev[0].xarr, dev[0].T, dev[0].BuffL);
	}
	else
	{
		cudaSetDevice(0);
		bound_x_left<<< dim3(dev[0].yarr,dev[0].zarr,1) , dim3(xpad,1,1) >>> (dev[0], dev[0].T);
	}
	syncallstreams(dev);

	////////////////////////////////

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);
		buff_left<<< gdim, bdim, 0, dev[n].stream >>> (dev[n].xres, dev[n].xarr, dev[n].T, dev[n].BuffL);
	}
	for (int n=1; n<ndev; n++)
	{
		cudaMemcpyAsync( dev[n-1].BuffR, dev[n].BuffL, xpad*yarr*zarr*sizeof(Cell), cudaMemcpyDeviceToDevice, dev[n].stream);
	}
		cudaMemcpyAsync( dev[ndev-1].BuffR, dev[0].BuffL, xpad*yarr*zarr*sizeof(Cell), cudaMemcpyDeviceToDevice, dev[0].stream);
	syncallstreams(dev);

	for (int n=1; n<ndev; n++)
	{
		cudaSetDevice(n-1);
		dump_rght<<< gdim, bdim, 0, dev[n-1].stream >>> (dev[n-1].xres, dev[n-1].xarr, dev[n-1].T, dev[n-1].BuffR);
	}
	syncallstreams(dev);

	////////////////////////////////

	if (bound_rgh == 3)
	{
		cudaSetDevice(ndev-1);
		dump_rght<<< gdim, bdim, 0, dev[ndev-1].stream >>> (dev[ndev-1].xres, dev[ndev-1].xarr, dev[ndev-1].T, dev[ndev-1].BuffR);
	}
	else
	{
		cudaSetDevice(ndev-1);
		bound_x_rght<<< dim3(dev[ndev-1].yarr,dev[ndev-1].zarr,1) , dim3(xpad,1,1) >>> (dev[ndev-1], dev[ndev-1].T);
	}
	syncallstreams(dev);

	return;
}

void boundy2(Grid* dev, double t=0.0)
{
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		boundy<<< dim3(dev[n].xarr,dev[n].zres,2), dim3(ypad,1,1), 0, dev[n].stream >>>(dev[n], dev[n].T, t);
	}

	syncallstreams(dev);

	return;
}

void boundz2(Grid* dev)
{
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		boundz<<< dim3(dev[n].xarr,dev[n].yarr,2), dim3(zpad,1,1), 0, dev[n].stream >>>(dev[n], dev[n].T);
	}

	syncallstreams(dev);

	return;
}
