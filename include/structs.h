#ifndef STRUCTS_H
#define STRUCTS_H

struct State
{
	double rl;
	double pl;
	double ul;
	double vl;
	double wl;

	double rr;
	double pr;
	double ur;
	double vr;
	double wr;
};


struct Cell
{
	double r;
	double p;
	double u;
	double v;
	double w;
	__host__ __device__ void copy(Cell B)
	{
		r = B.r;
		p = B.p;
		u = B.u;
		v = B.v;
		w = B.w;
	}
	__host__ __device__ void copy1(Cell B)
	{
		r = B.r;
		p = B.p;
		u = B.v;
		v = B.w;
		w = B.u;
	}
	__host__ __device__ void copy2(Cell B)
	{
		r = B.r;
		p = B.p;
		u = B.w;
		v = B.u;
		w = B.v;
	}
	__host__ __device__ void add(Cell B)
	{
		r += B.r;
		p += B.p;
		u += B.u;
		v += B.v;
		w += B.w;
	}
	__host__ __device__ void multiply(double x)
	{
		r *= x;
		p *= x;
		u *= x;
		v *= x;
		w *= x;
	}
	__host__ __device__ void zero()
	{
		r = 0.0;
		p = 0.0;
		u = 0.0;
		v = 0.0;
		w = 0.0;
	}
};

struct body
{
	double m;

	double x;
	double y;
	double z;

	double vx;
	double vy;
	double vz;

	double rs;
};

struct Grid
{
	cudaStream_t stream;

	int xbgn;
	int xres;
	int xarr;
	
	int ybgn;
	int yres;
	int yarr;
	
	int zbgn;
	int zres;
	int zarr;

	double* xa;
	double* xv;

	double* ya;
	double* yv;

	double* za;
	double* zv;

	double* dt;
	double* Buff;

	double* orb_rot;
	double* orb_res;
	int* orb_shf;
	int* orb_inc;
	int* shfL;
	int* shfR;

	Cell* C;
	Cell* T;
	Cell* F;

	Cell* BuffL;
	Cell* BuffR;

	double* vis_tensor;
	double* fx;
	double* fy;
	double* fz;

	double* Du;
	double* Dv;
	double* Dw;
	double* De;

	double* ext;

	body* planets;

	__host__ __device__ int get_ind(int i, int j, int k)
	{
		return i + xarr*j + xarr*yarr*k;
	}

	////////////////////////////////////////
	__host__ __device__ double get_xa(int i)
	{
		return xa[i+xbgn];
	}
	__host__ __device__ double get_xv(int i)
	{
		return xv[i+xbgn];
	}
	__host__ __device__ double get_xc(int i)
	{
		return 0.5*(xa[i+1+xbgn]+xa[i+xbgn]);
	}
	////////////////////////////////////////
	__host__ __device__ double get_ya(int j)
	{
		return ya[j+ybgn];
	}
	__host__ __device__ double get_yv(int j)
	{
		return yv[j+ybgn];
	}
	__host__ __device__ double get_yc(int j)
	{
		return 0.5*(ya[j+1+ybgn]+ya[j+ybgn]);
	}
	////////////////////////////////////////
	__host__ __device__ double get_za(int k)
	{
		return za[k+zbgn];
	}
	__host__ __device__ double get_zv(int k)
	{
		return zv[k+zbgn];
	}
	__host__ __device__ double get_zc(int k)
	{
		return 0.5*(za[k+1+zbgn]+za[k+zbgn]);
	}
	////////////////////////////////////////
	__host__ __device__ double get_r(int i, int j, int k)
	{
		return C[get_ind(i,j,k)].r;
	}
	__host__ __device__ double get_p(int i, int j, int k)
	{
		return C[get_ind(i,j,k)].p;
	}
	__host__ __device__ double get_u(int i, int j, int k)
	{
		return C[get_ind(i,j,k)].u;
	}
	__host__ __device__ double get_v(int i, int j, int k)
	{
		return C[get_ind(i,j,k)].v;
	}
	__host__ __device__ double get_w(int i, int j, int k)
	{
		return C[get_ind(i,j,k)].w;
	}
	////////////////////////////////////////
	__host__ __device__ Cell get_cell(int i, int j, int k)
	{
		return C[get_ind(i,j,k)];
	}
	__host__ __device__ void write_cell(int i, int j, int k, Cell W)
	{
		C[get_ind(i,j,k)].copy(W);
		return;
	}
	////////////////////////////////////////
	__host__ __device__ void write_r(int i, int j, int k, double r)
	{
		C[get_ind(i,j,k)].r=r;
		return;
	}
	__host__ __device__ void write_p(int i, int j, int k, double p)
	{
		C[get_ind(i,j,k)].p=p;
		return;
	}
	__host__ __device__ void write_u(int i, int j, int k, double u)
	{
		C[get_ind(i,j,k)].u=u;
		return;
	}
	__host__ __device__ void write_v(int i, int j, int k, double v)
	{
		C[get_ind(i,j,k)].v=v;
		return;
	}
	__host__ __device__ void write_w(int i, int j, int k, double w)
	{
		C[get_ind(i,j,k)].w=w;
		return;
	}
	////////////////////////////////////////
	__host__ __device__ double get_rot(int i, int k)
	{
		return orb_rot[i+xarr*k];
	}
	__host__ __device__ void write_rot(int i, int k, double rot)
	{
		orb_rot[i+xarr*k] = rot;
		return;
	}
	////////////////////////////////////////
	__host__ __device__ double get_res(int i, int k)
	{
		return orb_res[i+xarr*k];
	}
	__host__ __device__ void write_res(int i, int k, double res)
	{
		orb_res[i+xarr*k] = res;
		return;
	}
	////////////////////////////////////////
	__host__ __device__ int get_shf(int i, int k)
	{
		return orb_shf[i+xarr*k];
	}
	__host__ __device__ void write_shf(int i, int k, int shf)
	{
		orb_shf[i+xarr*k] = shf;
		return;
	}
	////////////////////////////////////////
	__host__ __device__ int get_inc(int i, int k)
	{
		return orb_inc[i+xarr*k];
	}
	__host__ __device__ void write_inc(int i, int k, int inc)
	{
		orb_inc[i+xarr*k] = inc;
		return;
	}
	////////////////////////////////////////
};

/*
struct GPU_plan
{
  int id;

  int N_ring;

  int istart;
  int iblk;

  int kstart;
  int kblk;

  long int memsize;

  ring *rings;
  ring *h_rings;
  
  ConsVar *D;
  #ifdef dust_flag
  ConsDus *F;
  #endif

  ring *lft;
  ring *rgh;
  ring *udr;
  ring *top;

  ring_geometry *ygrid;
  ring_geometry *h_ygrid;

  PrimVar *sym;
  #ifdef dust_flag
  PrimDus *sym_d;
  #endif

  dim3 sx_grid;
  dim3 sy_grid;
  dim3 sz_grid;

  int sx_block;
  int sy_block;
  int sz_block;

  dim3 t_grid;

  double *dt;
  double *h_dt;
  double *dt_reduc;

  body *planet;

  double *h_output;
  double *d_output;

  cudaEvent_t event;

  cudaStream_t stream;
};
*/
#endif
