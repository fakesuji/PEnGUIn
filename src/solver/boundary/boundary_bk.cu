__global__ void boundx(double* xa, double* ya, double* za, int nx, int mx, int my, Cell* G)
{
	int n = threadIdx.x;
	int j = blockIdx.x + ypad;
	int k = blockIdx.y + zpad;
	int ib = blockIdx.z;

	Cell* C = &G[mx*(j+my*k)];

	if (ib==0)
	{
		if (bound_lft == 3)
		{
			C[n].copy(C[nx+n]);
		}
		else if (bound_lft == 2)
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
			C[n].copy(init_C(xa[n],xa[n+1],ya[j],ya[j+1],za[k],za[k+1]));
		}
	}
	else if (ib==1)
	{
		if (bound_rgh == 3)
		{
			C[nx+xpad+n].copy(C[xpad+n]);
		}
		else if (bound_rgh == 2)
		{
			C[nx+xpad+n].copy(C[nx+xpad-1-n]);
			C[nx+xpad+n].u *= -1.0;
		}
		else if (bound_rgh == 1)
		{
			C[nx+xpad+n].copy(C[nx+xpad-1]);
		}
		else if (bound_rgh == 0)
		{
			C[nx+xpad+n].copy(init_C(xa[nx+xpad+n],xa[nx+xpad+n+1],ya[j],ya[j+1],za[k],za[k+1]));
		}
	}
	return;
}

__global__ void bound_x_left(double* xa, double* ya, double* za, int nx, int mx, int my, Cell* G)
{
	int n = threadIdx.x;
	int j = blockIdx.x + ypad;
	int k = blockIdx.y + zpad;

	Cell* C = &G[mx*(j+my*k)];

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
		C[n].copy(init_C(xa[n],xa[n+1],ya[j],ya[j+1],za[k],za[k+1]));
	}
	return;
}

__global__ void bound_x_rght(double* xa, double* ya, double* za, int nx, int mx, int my, Cell* G)
{
	int n = threadIdx.x;
	int j = blockIdx.x + ypad;
	int k = blockIdx.y + zpad;

	Cell* C = &G[mx*(j+my*k)];

	if (bound_rgh == 2)
	{
		C[nx+xpad+n].copy(C[nx+xpad-1-n]);
		C[nx+xpad+n].u *= -1.0;
	}
	else if (bound_rgh == 1)
	{
		C[nx+xpad+n].copy(C[nx+xpad-1]);
	}
	else if (bound_rgh == 0)
	{
		C[nx+xpad+n].copy(init_C(xa[nx+xpad+n],xa[nx+xpad+n+1],ya[j],ya[j+1],za[k],za[k+1]));
	}
	//C[nx+xpad+n].r=0.0;
	return;
}

__global__ void boundy(double* xa, double* ya, double* za, int ny, int mx, int my, Cell* C)
{
	int n = threadIdx.x;
	int i = blockIdx.x + xpad;
	int k = blockIdx.y + zpad;
	int ib = blockIdx.z;

	if (ib==0)
	{
		if (bound_bak == 3)
		{
			C[glo_idx(i,n,k,mx,my)].copy(C[glo_idx(i,ny+n,k,mx,my)]);
		}
		else if (bound_bak == 2)
		{
			C[glo_idx(i,n,k,mx,my)].copy(C[glo_idx(i,2*ypad-1-n,k,mx,my)]);
			C[glo_idx(i,n,k,mx,my)].v *= -1.0;
		}
		else if (bound_bak == 1)
		{
			C[glo_idx(i,n,k,mx,my)].copy(C[glo_idx(i,ypad,k,mx,my)]);
		}
		else if (bound_bak == 0)
		{
			C[glo_idx(i,n,k,mx,my)].copy(init_C(xa[i],xa[i+1],ya[n],ya[n+1],za[k],za[k+1]));
		}
	}
	else if (ib==1)
	{
		if (bound_frn == 3)
		{
			C[glo_idx(i,ny+ypad+n,k,mx,my)].copy(C[glo_idx(i,ypad+n,k,mx,my)]);
		}
		else if (bound_frn == 2)
		{
			C[glo_idx(i,ny+ypad+n,k,mx,my)].copy(C[glo_idx(i,ny+ypad-1-n,k,mx,my)]);
			C[glo_idx(i,ny+ypad+n,k,mx,my)].v *= -1.0;
		}
		else if (bound_frn == 1)
		{
			C[glo_idx(i,ny+ypad+n,k,mx,my)].copy(C[glo_idx(i,ny+ypad-1,k,mx,my)]);
		}
		else if (bound_frn == 0)
		{
			C[glo_idx(i,ny+ypad+n,k,mx,my)].copy(init_C(xa[i],xa[i+1],ya[ny+ypad+n],ya[ny+ypad+n+1],za[k],za[k+1]));
		}
	}
	return;
}

__global__ void boundz(Grid G);
{
	int n = threadIdx.x;
	int i = blockIdx.x + xpad;
	int j = blockIdx.y + ypad;
	int ib = blockIdx.z;

	if (ib==0)
	{
		if (bound_bom == 3)
		{
			//G.write_shf(i,j,n,0);
			//G.write_cell(i,j,n,G.get_cell(i,j,G.zres+n));
			C[glo_idx(i,j,n,mx,my)].copy(C[glo_idx(i,j,nz+n,mx,my)]);
		}
		else if (bound_bom == 2)
		{
			C[glo_idx(i,j,n,mx,my)].copy(C[glo_idx(i,j,2*zpad-1-n,mx,my)]);
			C[glo_idx(i,j,n,mx,my)].w *= -1.0;
		}
		else if (bound_bom == 1)
		{
			C[glo_idx(i,j,n,mx,my)].copy(C[glo_idx(i,j,zpad,mx,my)]);
		}
		else if (bound_bom == 0)
		{
			C[glo_idx(i,j,n,mx,my)].copy(init_C(xa[i],xa[i+1],ya[j],ya[j+1],za[n],za[n+1]));
		}
	}
	else if (ib==1)
	{
		if (bound_top == 3)
		{
			C[glo_idx(i,j,nz+zpad+n,mx,my)].copy(C[glo_idx(i,j,zpad+n,mx,my)]);
		}
		else if (bound_top == 2)
		{
			C[glo_idx(i,j,nz+zpad+n,mx,my)].copy(C[glo_idx(i,j,nz+zpad-1-n,mx,my)]);
			C[glo_idx(i,j,nz+zpad+n,mx,my)].w *= -1.0;
		}
		else if (bound_top == 1)
		{
			C[glo_idx(i,j,nz+zpad+n,mx,my)].copy(C[glo_idx(i,j,nz+zpad-1,mx,my)]);
		}
		else if (bound_top == 0)
		{
			C[glo_idx(i,j,nz+zpad+n,mx,my)].copy(init_C(xa[i],xa[i+1],ya[j],ya[j+1],za[nz+zpad+n],za[nz+zpad+n+1]));
		}
	}
	return;
}

//=========================================================================================

__global__ void buff_left(int x_res, int x_arr, Cell* C, Cell* B)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int k = threadIdx.z + blockIdx.z*blockDim.z;

	B[i + xpad*(j + yres*k)].copy(C[i+xpad + x_arr*(j+ypad + yarr*(k+zpad))]);
}

__global__ void buff_rght(int x_res, int x_arr, Cell* C, Cell* B)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int k = threadIdx.z + blockIdx.z*blockDim.z;

	B[i + xpad*(j + yres*k)].copy(C[i + x_res + x_arr*(j+ypad + yarr*(k+zpad))]);
}

__global__ void dump_left(int x_res, int x_arr, Cell* C, Cell* B)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int k = threadIdx.z + blockIdx.z*blockDim.z;

	C[i + x_arr*(j+ypad + yarr*(k+zpad))].copy(B[i + xpad*(j + yres*k)]);
}

__global__ void dump_rght(int x_res, int x_arr, Cell* C, Cell* B)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int k = threadIdx.z + blockIdx.z*blockDim.z;

	C[i+x_res+xpad + x_arr*(j+ypad + yarr*(k+zpad))].copy(B[i + xpad*(j + yres*k)]);
}

void boundx(Grid* dev)
{
	int gdimy, gdimz, bdimy, bdimz;
	if (yres%16==0)
	{
		gdimy = yres/16;
		bdimy = 16;
	}
	else
	{
		gdimy = yres;
		bdimy = 1;
	}
	if (zres%16==0)
	{
		gdimz = zres/16;
		bdimz = 16;
	}
	else
	{
		gdimz = zres;
		bdimz = 1;
	}

	dim3 gdim(1,gdimy,gdimz);
	dim3 bdim(xpad,bdimy,bdimz);

	////////////////////////////////

	for (int n=0; n<ndev-1; n++)
	{
		cudaSetDevice(n);
		buff_rght<<< gdim, bdim, 0, dev[n].stream >>> (dev[n].xres, dev[n].xarr, dev[n].C, dev[n].BuffR);
		cudaMemcpyAsync( dev[n+1].BuffL, dev[n].BuffR, xpad*yres*zres*sizeof(Cell), cudaMemcpyDeviceToDevice, dev[n].stream);
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
	for (int n=0; n<ndev-1; n++)
	{
		cudaSetDevice(n+1);
		dump_left<<< gdim, bdim, 0, dev[n+1].stream >>> (dev[n+1].xres, dev[n+1].xarr, dev[n+1].C, dev[n+1].BuffL);
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);

	////////////////////////////////

	for (int n=1; n<ndev; n++)
	{
		cudaSetDevice(n);
		buff_left<<< gdim, bdim, 0, dev[n].stream >>> (dev[n].xres, dev[n].xarr, dev[n].C, dev[n].BuffL);
		cudaMemcpyAsync( dev[n-1].BuffR, dev[n].BuffL, xpad*yres*zres*sizeof(Cell), cudaMemcpyDeviceToDevice, dev[n].stream);
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
	for (int n=1; n<ndev; n++)
	{
		cudaSetDevice(n-1);
		dump_rght<<< gdim, bdim, 0, dev[n-1].stream >>> (dev[n-1].xres, dev[n-1].xarr, dev[n-1].C, dev[n-1].BuffR);
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);

	////////////////////////////////

	if (bound_rgh == 3)
	{
		cudaSetDevice(0);
		buff_left<<< gdim, bdim, 0, dev[0].stream >>> (dev[0].xres, dev[0].xarr, dev[0].C, dev[0].BuffL);
		cudaMemcpyAsync( dev[ndev-1].BuffR, dev[0].BuffL, xpad*yres*zres*sizeof(Cell), cudaMemcpyDeviceToDevice, dev[0].stream);
		cudaStreamSynchronize(dev[0].stream);
		cudaSetDevice(ndev-1);
		dump_rght<<< gdim, bdim, 0, dev[ndev-1].stream >>> (dev[ndev-1].xres, dev[ndev-1].xarr, dev[ndev-1].C, dev[ndev-1].BuffR);
	}
	else
	{
		cudaSetDevice(ndev-1);
		bound_x_rght<<< dim3(dev[ndev-1].yres,dev[ndev-1].zres,1) , dim3(xpad,1,1) >>>
		            (&dev[ndev-1].xa[dev[ndev-1].xbgn], dev[ndev-1].ya, dev[ndev-1].za, dev[ndev-1].xres, dev[ndev-1].xarr, dev[ndev-1].yarr, dev[ndev-1].C);
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);

	if (bound_lft == 3)
	{
		cudaSetDevice(ndev-1);
		buff_rght<<< gdim, bdim, 0, dev[ndev-1].stream >>> (dev[ndev-1].xres, dev[ndev-1].xarr, dev[ndev-1].C, dev[ndev-1].BuffR);
		cudaMemcpyAsync( dev[0].BuffL, dev[ndev-1].BuffR, xpad*yres*zres*sizeof(Cell), cudaMemcpyDeviceToDevice, dev[ndev-1].stream);
		cudaStreamSynchronize(dev[ndev-1].stream);
		cudaSetDevice(0);
		dump_left<<< gdim, bdim, 0, dev[0].stream >>> (dev[0].xres, dev[0].xarr, dev[0].C, dev[0].BuffL);
	}
	else
	{
		cudaSetDevice(0);
		bound_x_left<<< dim3(dev[0].yres,dev[0].zres,1) , dim3(xpad,1,1) >>>
		            (&dev[0].xa[dev[0].xbgn], dev[0].ya, dev[0].za, dev[0].xres, dev[0].xarr, dev[0].yarr, dev[0].C);
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);


	return;
}

//=========================================================================================

__global__ void buff_shf_left(int x_res, int x_arr, int* S, int* B)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int k = threadIdx.y + blockIdx.y*blockDim.y;

	B[i + xpad*k] = S[i+xpad + x_arr*(k+zpad)];
}

__global__ void buff_shf_rght(int x_res, int x_arr, int* S, int* B)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int k = threadIdx.y + blockIdx.y*blockDim.y;

	B[i + xpad*k] = S[i+x_res + x_arr*(k+zpad)];
}

__global__ void dump_shf_left(int x_res, int x_arr, int* S, int* B)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int k = threadIdx.y + blockIdx.y*blockDim.y;

	S[i + x_arr*(k+zpad)] = B[i + xpad*k];
}

__global__ void dump_shf_rght(int x_res, int x_arr, int* S, int* B)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int k = threadIdx.y + blockIdx.y*blockDim.y;

	S[i+x_res+xpad + x_arr*(k+zpad)] = B[i + xpad*k];
}

void orb_boundx(Grid* dev)
{
	int bdimz, gdimz;
	if (zres%16==0)
	{
		gdimz = zres/16;
		bdimz = 16;
	}
	else
	{
		gdimz = zres;
		bdimz = 1;
	}

	dim3 gdim(1,gdimz,1);
	dim3 bdim(xpad,bdimz,1);

	////////////////////////////////

	for (int n=0; n<ndev-1; n++)
	{
		cudaSetDevice(n);
		buff_shf_rght<<< gdim, bdim, 0, dev[n].stream >>> (dev[n].xres, dev[n].xarr, dev[n].orb_shf, dev[n].shfR);
		cudaMemcpyAsync( dev[n+1].shfL, dev[n].shfR, xpad*zres*sizeof(int), cudaMemcpyDeviceToDevice, dev[n].stream);
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
	for (int n=0; n<ndev-1; n++)
	{
		cudaSetDevice(n+1);
		dump_shf_left<<< gdim, bdim, 0, dev[n+1].stream >>> (dev[n+1].xres, dev[n+1].xarr, dev[n+1].orb_shf, dev[n+1].shfL);
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);

	////////////////////////////////

	for (int n=1; n<ndev; n++)
	{
		cudaSetDevice(n);
		buff_shf_left<<< gdim, bdim, 0, dev[n].stream >>> (dev[n].xres, dev[n].xarr, dev[n].orb_shf, dev[n].shfL);
		cudaMemcpyAsync( dev[n-1].shfR, dev[n].shfL, xpad*zres*sizeof(int), cudaMemcpyDeviceToDevice, dev[n].stream);
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
	for (int n=1; n<ndev; n++)
	{
		cudaSetDevice(n-1);
		dump_shf_rght<<< gdim, bdim, 0, dev[n-1].stream >>> (dev[n-1].xres, dev[n-1].xarr, dev[n-1].orb_shf, dev[n-1].shfR);
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);


	return;
}
