__global__ void bound_x_left(Grid G)
{
	int n = threadIdx.x;
	int j = blockIdx.x;
	int k = blockIdx.y;

	Cell* C = &G.C[G.xarr*(j+G.yarr*k)];

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

__global__ void bound_x_rght(Grid G)
{
	int nx = G.xres+xpad;
	int n = threadIdx.x;
	int j = blockIdx.x;
	int k = blockIdx.y;

	Cell* C = &G.C[G.xarr*(j+G.yarr*k)];

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

__global__ void boundy(Grid G)
{
	int n = threadIdx.x;
	int i = blockIdx.x + xpad;
	int k = blockIdx.y;
	int ib = blockIdx.z;
	Cell T;

	if (ib==0)
	{
		if (bound_bak == 3)
		{
			T.copy(G.get_cell(i,G.yres+n,k));
		}
		else if (bound_bak == 2)
		{
			T.copy(G.get_cell(i,2*ypad-1-n,k));
			T.v *= -1.0;
		}
		else if (bound_bak == 1)
		{
			T.copy(G.get_cell(i,ypad,k));
		}
		else if (bound_bak == 0)
		{
			T = init_C(G.get_xc(i),G.get_yc(n),G.get_zc(k));
		}
		G.C[G.get_ind(i,n,k)].copy(T);
	}
	else if (ib==1)
	{
		if (bound_frn == 3)
		{
			T.copy(G.get_cell(i,ypad+n,k));
		}
		else if (bound_frn == 2)
		{
			T.copy(G.get_cell(i,G.yres+ypad-1-n,k));
			T.v *= -1.0;
		}
		else if (bound_frn == 1)
		{
			T.copy(G.get_cell(i,G.yres+ypad-1,k));
		}
		else if (bound_frn == 0)
		{
			T = init_C(G.get_xc(i),G.get_yc(G.yres+ypad+n),G.get_zc(k));
		}
		G.C[G.get_ind(i,G.yres+ypad+n,k)].copy(T);
	}
	return;
}

__global__ void boundz(Grid G)
{
	int n = threadIdx.x;
	int i = blockIdx.x + xpad;
	int j = blockIdx.y + ypad;
	int ib = blockIdx.z;
	Cell T;

	if (ib==0)
	{
		if (bound_bom == 3)
		{
			T = G.get_cell(i,j,G.zres+n);
		}
		else if (bound_bom == 2)
		{
			T = G.get_cell(i,j,2*zpad-1-n);
			T.w *= -1.0;
		}
		else if (bound_bom == 1)
		{
			T = G.get_cell(i,j,zpad);
		}
		else if (bound_bom == 0)
		{
			T = init_C(G.get_xc(i),G.get_yc(j),G.get_zc(n));
		}
		G.write_cell(i,j,n,T);
	}
	else if (ib==1)
	{
		if (bound_top == 3)
		{
			T = G.get_cell(i,j,zpad+n);
		}
		else if (bound_top == 2)
		{
			T = G.get_cell(i,j,G.zres+zpad-1-n);
			T.w *= -1.0;
		}
		else if (bound_top == 1)
		{
			T = G.get_cell(i,j,G.zres+zpad-1);
		}
		else if (bound_top == 0)
		{
			T = init_C(G.get_xc(i),G.get_yc(j),G.get_zc(G.zres+zpad+n));
		}
		G.write_cell(i,j,G.zres+zpad+n,T);
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
		bound_x_rght<<< dim3(dev[ndev-1].yarr,dev[ndev-1].zarr,1) , dim3(xpad,1,1) >>> (dev[ndev-1]);
	}

	if (bound_lft == 3)
	{
		cudaSetDevice(0);
		dump_left<<< gdim, bdim, 0, dev[0].stream >>> (dev[0].xres, dev[0].xarr, dev[0].C, dev[0].BuffL);
	}
	else
	{
		cudaSetDevice(0);
		bound_x_left<<< dim3(dev[0].yarr,dev[0].zarr,1) , dim3(xpad,1,1) >>> (dev[0]);
	}
	syncallstreams(dev);

	return;
}

void boundy(Grid* dev)
{
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		boundy<<< dim3(dev[n].xres,dev[n].zarr,2), dim3(ypad,1,1), 0, dev[n].stream >>>(dev[n]);
	}

	syncallstreams(dev);

	return;
}

void boundz(Grid* dev)
{
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		boundz<<< dim3(dev[n].xres,dev[n].yres,2), dim3(zpad,1,1), 0, dev[n].stream >>>(dev[n]);
	}

	syncallstreams(dev);

	return;
}
