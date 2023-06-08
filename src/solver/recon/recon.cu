#include "PLM.cu"
#include "PPM.cu"
#include "PEM.cu"

//=======================================================================================

__device__ void get_2nd_MOC_parameters(int i, int geom, double* r, double* dr, double* dv, double* a, double* par)
{
	double a1, a2, a3, x1, x2, x3, tmp;

	a1 = a[i-1];
	a2 = a[i];
	a3 = a[i+1];
	
	if ((a3-a2)*(a2-a1)<=0.0)
	{
		par[0] = 0.0;
		par[1] = a2;
		par[2] = 0.0;
		return;	
	}

	x1 = dr[i-1];
	x2 = dr[i];
	x3 = dr[i+1];

	//===============================================================================

	tmp = (a3-a1)*x2/(x1+x3+2.0*x2);
  	tmp = copysign( fmin(fmin(fabs(a2-a1),fabs(a3-a2)),fabs(tmp)) , tmp );

	par[0] = tmp;
	par[1] = a2;
	par[2] = tmp;
    
	return;
}

//=======================================================================================

__device__ void get_2nd_VAN_parameters(int i, int geom, double* r, double* dr, double* dv, double* a, double* par)
{
	double a1, a2, a3, x1, x2, x3, tmp;

	a1 = a[i-1];
	a2 = a[i];
	a3 = a[i+1];
	
	if ((a3-a2)*(a2-a1)<=0.0)
	{
		par[0] = 0.0;
		par[1] = a2;
		par[2] = 0.0;
		return;	
	}

	x1 = dr[i-1];
	x2 = dr[i];
	x3 = dr[i+1];

	//===============================================================================

	tmp = (a3-a2)*(a2-a1)/(a3-a1);
	tmp*= x2*(x1+x3+2.0*x2)/(x3+x2)/(x2+x1);

	par[0] = tmp;
	par[1] = a2;
	par[2] = tmp;
    
	return;
}

//=======================================================================================

__device__ void get_2nd_parameters(int i, int geom, double* x, double* dx, double* dv, double* a, double* par)
{
	double a0, a1, a2, sL, sR, c1, x0, x1, x2, x3;

	a0 = a[i-1];
	a1 = a[i];
	a2 = a[i+1];

	x0 = x[i-1];
	x1 = x[i];
	x2 = x[i+1];
	x3 = x[i+2];

	double B02,B12,B22;
	double D02;

	if (geom==1)
	{
		B02 = (x1*x1*x1-x0*x0*x0)/(x1*x1-x0*x0)/1.5;
		B12 = (x2*x2*x2-x1*x1*x1)/(x2*x2-x1*x1)/1.5;
		B22 = (x3*x3*x3-x2*x2*x2)/(x3*x3-x2*x2)/1.5;
	}	
	else
	{
		B02 = (x1*x1-x0*x0)/(x1-x0)/2.0;
		B12 = (x2*x2-x1*x1)/(x2-x1)/2.0;
		B22 = (x3*x3-x2*x2)/(x3-x2)/2.0;
	}

	D02 = B12-B02;

	c1 = (a1-a0)/D02;
	sL = c1*(B12-x1);

	//=====================================================

	D02 = B22-B12;

	c1 = (a2-a1)/D02;
	sR = c1*(x1-B12);

	sL = copysign(fmin(fabs(sL),fabs(a1-a0)),a1-a0);
	sR = copysign(fmin(fabs(sR),fabs(a2-a1)),a2-a1);

	if (sR*sL<=0.0)
	{
		par[0] = 0.0;
		par[1] = a1;
		par[2] = 0.0;
	}
	else
	{
		par[0] = sL;
		par[1] = a1;
		par[2] = sR;
	}

	return;
}

//=======================================================================================

__device__ void get_3rd_parameters(int i, int geom, double* x, double* dx, double* dv, double* a, double* par)
{
	double a0, a1, a2, sL, sR, c1, c2, x0, x1, x2, x3;

	a0 = a[i-1];
	a1 = a[i];
	a2 = a[i+1];

	x0 = x[i-1];
	x1 = x[i];
	x2 = x[i+1];
	x3 = x[i+2];

	double B02,B12,B22,B03,B13,B23;
	double D02,D12,D03,D13;

	if (geom==2)
	{
		B02 = (x1*x1*x1*x1-x0*x0*x0*x0)/(x1*x1*x1-x0*x0*x0)/(4.0/3.0);
		B12 = (x2*x2*x2*x2-x1*x1*x1*x1)/(x2*x2*x2-x1*x1*x1)/(4.0/3.0);
		B22 = (x3*x3*x3*x3-x2*x2*x2*x2)/(x3*x3*x3-x2*x2*x2)/(4.0/3.0);

		B03 = (x1*x1*x1*x1*x1-x0*x0*x0*x0*x0)/(x1*x1*x1-x0*x0*x0)/(5.0/3.0);
		B13 = (x2*x2*x2*x2*x2-x1*x1*x1*x1*x1)/(x2*x2*x2-x1*x1*x1)/(5.0/3.0);
		B23 = (x3*x3*x3*x3*x3-x2*x2*x2*x2*x2)/(x3*x3*x3-x2*x2*x2)/(5.0/3.0);
	}
	else if (geom==1)
	{
		B02 = (x1*x1*x1-x0*x0*x0)/(x1*x1-x0*x0)/1.5;
		B12 = (x2*x2*x2-x1*x1*x1)/(x2*x2-x1*x1)/1.5;
		B22 = (x3*x3*x3-x2*x2*x2)/(x3*x3-x2*x2)/1.5;

		B03 = (x1*x1*x1*x1-x0*x0*x0*x0)/(x1*x1-x0*x0)/2.0;
		B13 = (x2*x2*x2*x2-x1*x1*x1*x1)/(x2*x2-x1*x1)/2.0;
		B23 = (x3*x3*x3*x3-x2*x2*x2*x2)/(x3*x3-x2*x2)/2.0;
	}	
	else
	{
		B02 = (x1*x1-x0*x0)/(x1-x0)/2.0;
		B12 = (x2*x2-x1*x1)/(x2-x1)/2.0;
		B22 = (x3*x3-x2*x2)/(x3-x2)/2.0;

		B03 = (x1*x1*x1-x0*x0*x0)/(x1-x0)/3.0;
		B13 = (x2*x2*x2-x1*x1*x1)/(x2-x1)/3.0;
		B23 = (x3*x3*x3-x2*x2*x2)/(x3-x2)/3.0;
	}

	D02 = B12-B02;
	D12 = B22-B12;
	D03 = B13-B03;
	D13 = B23-B13;

	c1 = ( D03*(a2-a1) - D13*(a1-a0))/(D12*D03-D02*D13);
	c2 = (-D02*(a2-a1) + D12*(a1-a0))/(D12*D03-D02*D13);
	
	sL = c1*(B12-x1) + c2*(B13-x1*x1);
	sR = c1*(x2-B12) + c2*(x2*x2-B13);

	sL = copysign(fmin(fabs(sL),fabs(a1-a0)),a1-a0);
	sR = copysign(fmin(fabs(sR),fabs(a2-a1)),a2-a1);

	if (sR*sL<=0.0)
	{
		par[0] = 0.0;
		par[1] = a1;
		par[2] = 0.0;
	}
	else
	{
		par[0] = sL;
		par[1] = a1;
		par[2] = sR;
	}

	return;
}

//=======================================================================================

__device__ void get_4th_parameters(int i, int geom, double* x, double* dx, double* dv, double* a, double* par)
{
	double a0, a1, a2, a3, a4, sL, sR, c1, c2, c3, x0, x1, x2, x3, x4, x5;

	a0 = a[i-2];
	a1 = a[i-1];
	a2 = a[i];
	a3 = a[i+1];
	a4 = a[i+2];

	x0 = x[i-2];
	x1 = x[i-1];
	x2 = x[i];
	x3 = x[i+1];
	x4 = x[i+2];
	x5 = x[i+3];

	double B02,B12,B22,B32,B42;
	double B03,B13,B23,B33,B43;
	double B04,B14,B24,B34,B44;

	double D02,D12,D22;
	double D03,D13,D23;
	double D04,D14,D24;


	if (geom==2)
	{
		B02 = (x1*x1*x1*x1-x0*x0*x0*x0)/(x1*x1*x1-x0*x0*x0)/(4.0/3.0);
		B12 = (x2*x2*x2*x2-x1*x1*x1*x1)/(x2*x2*x2-x1*x1*x1)/(4.0/3.0);
		B22 = (x3*x3*x3*x3-x2*x2*x2*x2)/(x3*x3*x3-x2*x2*x2)/(4.0/3.0);
		B32 = (x4*x4*x4*x4-x3*x3*x3*x3)/(x4*x4*x4-x3*x3*x3)/(4.0/3.0);
		B42 = (x5*x5*x5*x5-x4*x4*x4*x4)/(x5*x5*x5-x4*x4*x4)/(4.0/3.0);

		B03 = (x1*x1*x1*x1*x1-x0*x0*x0*x0*x0)/(x1*x1*x1-x0*x0*x0)/(5.0/3.0);
		B13 = (x2*x2*x2*x2*x2-x1*x1*x1*x1*x1)/(x2*x2*x2-x1*x1*x1)/(5.0/3.0);
		B23 = (x3*x3*x3*x3*x3-x2*x2*x2*x2*x2)/(x3*x3*x3-x2*x2*x2)/(5.0/3.0);
		B33 = (x4*x4*x4*x4*x4-x3*x3*x3*x3*x3)/(x4*x4*x4-x3*x3*x3)/(5.0/3.0);
		B43 = (x5*x5*x5*x5*x5-x4*x4*x4*x4*x4)/(x5*x5*x5-x4*x4*x4)/(5.0/3.0);

		B04 = (x1*x1*x1*x1*x1*x1-x0*x0*x0*x0*x0*x0)/(x1*x1*x1-x0*x0*x0)/2.0;
		B14 = (x2*x2*x2*x2*x2*x2-x1*x1*x1*x1*x1*x1)/(x2*x2*x2-x1*x1*x1)/2.0;
		B24 = (x3*x3*x3*x3*x3*x3-x2*x2*x2*x2*x2*x2)/(x3*x3*x3-x2*x2*x2)/2.0;
		B34 = (x4*x4*x4*x4*x4*x4-x3*x3*x3*x3*x3*x3)/(x4*x4*x4-x3*x3*x3)/2.0;
		B44 = (x5*x5*x5*x5*x5*x5-x4*x4*x4*x4*x4*x4)/(x5*x5*x5-x4*x4*x4)/2.0;
	}
	else if (geom==1)
	{
		B02 = (x1*x1*x1-x0*x0*x0)/(x1*x1-x0*x0)/1.5;
		B12 = (x2*x2*x2-x1*x1*x1)/(x2*x2-x1*x1)/1.5;
		B22 = (x3*x3*x3-x2*x2*x2)/(x3*x3-x2*x2)/1.5;
		B32 = (x4*x4*x4-x3*x3*x3)/(x4*x4-x3*x3)/1.5;
		B42 = (x5*x5*x5-x4*x4*x4)/(x5*x5-x4*x4)/1.5;

		B03 = (x1*x1*x1*x1-x0*x0*x0*x0)/(x1*x1-x0*x0)/2.0;
		B13 = (x2*x2*x2*x2-x1*x1*x1*x1)/(x2*x2-x1*x1)/2.0;
		B23 = (x3*x3*x3*x3-x2*x2*x2*x2)/(x3*x3-x2*x2)/2.0;
		B33 = (x4*x4*x4*x4-x3*x3*x3*x3)/(x4*x4-x3*x3)/2.0;
		B43 = (x5*x5*x5*x5-x4*x4*x4*x4)/(x5*x5-x4*x4)/2.0;

		B04 = (x1*x1*x1*x1*x1-x0*x0*x0*x0*x0)/(x1*x1-x0*x0)/2.5;
		B14 = (x2*x2*x2*x2*x2-x1*x1*x1*x1*x1)/(x2*x2-x1*x1)/2.5;
		B24 = (x3*x3*x3*x3*x3-x2*x2*x2*x2*x2)/(x3*x3-x2*x2)/2.5;
		B34 = (x4*x4*x4*x4*x4-x3*x3*x3*x3*x3)/(x4*x4-x3*x3)/2.5;
		B44 = (x5*x5*x5*x5*x5-x4*x4*x4*x4*x4)/(x5*x5-x4*x4)/2.5;
	}	
	else
	{
		B02 = (x1*x1-x0*x0)/(x1-x0)/2.0;
		B12 = (x2*x2-x1*x1)/(x2-x1)/2.0;
		B22 = (x3*x3-x2*x2)/(x3-x2)/2.0;
		B32 = (x4*x4-x3*x3)/(x4-x3)/2.0;
		B42 = (x5*x5-x4*x4)/(x5-x4)/2.0;

		B03 = (x1*x1*x1-x0*x0*x0)/(x1-x0)/3.0;
		B13 = (x2*x2*x2-x1*x1*x1)/(x2-x1)/3.0;
		B23 = (x3*x3*x3-x2*x2*x2)/(x3-x2)/3.0;
		B33 = (x4*x4*x4-x3*x3*x3)/(x4-x3)/3.0;
		B43 = (x5*x5*x5-x4*x4*x4)/(x5-x4)/3.0;

		B04 = (x1*x1*x1*x1-x0*x0*x0*x0)/(x1-x0)/4.0;
		B14 = (x2*x2*x2*x2-x1*x1*x1*x1)/(x2-x1)/4.0;
		B24 = (x3*x3*x3*x3-x2*x2*x2*x2)/(x3-x2)/4.0;
		B34 = (x4*x4*x4*x4-x3*x3*x3*x3)/(x4-x3)/4.0;
		B44 = (x5*x5*x5*x5-x4*x4*x4*x4)/(x5-x4)/4.0;
	}

	double tmp;

	D02 = B12-B02;
	D12 = B22-B12;
	D22 = B32-B22;

	D03 = B13-B03;
	D13 = B23-B13;
	D23 = B33-B23;

	D04 = B14-B04;
	D14 = B24-B14;
	D24 = B34-B24;

	//=================================================================================================

	c3 = (D02*(a3-a2) - D22*(a1-a0))/(D23*D02-D03*D22) - ((a2-a1)*D02 - (a1-a0)*D12)/(D13*D02 - D03*D12);
	tmp = (D02*D24 - D22*D04)/(D23*D02-D03*D22) - (D14*D02 - D04*D12)/(D13*D02 - D03*D12);

	c3 /= tmp;

	tmp = ((a2-a1)*D02 - (a1-a0)*D12)/(D13*D02 - D03*D12);
	c2 = tmp - c3*(D14*D02 - D04*D12)/(D13*D02 - D03*D12);

	tmp = (a1-a0)/D02;
	c1 = tmp - c2*D03/D02 - c3*D04/D02;
	
	sL = c1*(B22-x2) + c2*(B23-x2*x2) + c3*(B24-x2*x2*x2);

	//=================================================================================================

	D02 = B22-B12;
	D12 = B32-B22;
	D22 = B42-B32;

	D03 = B23-B13;
	D13 = B33-B23;
	D23 = B43-B33;

	D04 = B24-B14;
	D14 = B34-B24;
	D24 = B44-B34;

	c3 = (D02*(a4-a3) - D22*(a2-a1))/(D23*D02-D03*D22) - ((a3-a2)*D02 - (a2-a1)*D12)/(D13*D02 - D03*D12);
	tmp = (D02*D24 - D22*D04)/(D23*D02-D03*D22) - (D14*D02 - D04*D12)/(D13*D02 - D03*D12);

	c3 /= tmp;

	tmp = ((a3-a2)*D02 - (a2-a1)*D12)/(D13*D02 - D03*D12);
	c2 = tmp - c3*(D14*D02 - D04*D12)/(D13*D02 - D03*D12);

	tmp = (a2-a1)/D02;
	c1 = tmp - c2*D03/D02 - c3*D04/D02;
	
	sR = c1*(x3-B22) + c2*(x3*x3-B23) + c3*(x3*x3*x3-B24);

	//=================================================================================================

	sL = copysign(fmin(fabs(sL),fabs(a2-a1)),a2-a1);
	sR = copysign(fmin(fabs(sR),fabs(a3-a2)),a3-a2);

	if (sR*sL<=0.0)
	{
		par[0] = 0.0;
		par[1] = a2;
		par[2] = 0.0;
	}
	else
	{
		par[0] = sL;
		par[1] = a2;
		par[2] = sR;
	}

	return;
}

//=======================================================================================


__device__ void get_CON_parameters(int i, int geom, double* x, double* dx, double* dv, double* a, double* par)
{
	#if recon_flag==0
	get_2nd_VAN_parameters(i, geom, x, dx, dv, a, par);
	#elif recon_flag==1
	get_2nd_MOC_parameters(i, geom, x, dx, dv, a, par);
	#elif recon_flag==2
	get_2nd_parameters(i, geom, x, dx, dv, a, par);
	#elif recon_flag==3
	get_3rd_parameters(i, geom, x, dx, dv, a, par);
	#elif recon_flag==4
	get_3rd_parameters(i, geom, x, dx, dv, a, par);
	adjust_PPM_parameters(i,par);
	#elif recon_flag==5
	get_4th_parameters(i, geom, x, dx, dv, a, par);
	#elif recon_flag==6
	get_4th_parameters(i, geom, x, dx, dv, a, par);
	adjust_PPM_parameters(i,par);
	#endif
	return;
}

__device__ double get_CON_aveR(int geom, double x, double* par, double lx, double ly)
{
	#if recon_flag==0
	return get_PLM_aveR(geom, x, par);
	#elif recon_flag==1
	return get_PLM_aveR(geom, x, par);
	#elif recon_flag==2
	return get_PEM_aveR(geom, x, par);
	#elif recon_flag==3
	return get_PEM_aveR(geom, x, par);
	#elif recon_flag==4
	return get_PPM_aveR(geom, x, par);
	#elif recon_flag==5
	return get_PEM_aveR(geom, x, par);
	#elif recon_flag==6
	return get_PPM_aveR(geom, x, par);
	#endif
}

__device__ double get_CON_aveL(int geom, double x, double* par, double lx, double ly)
{
	#if recon_flag==0
	return get_PLM_aveL(geom, x, par);
	#elif recon_flag==1
	return get_PLM_aveL(geom, x, par);
	#elif recon_flag==2
	return get_PEM_aveL(geom, x, par);
	#elif recon_flag==3
	return get_PEM_aveL(geom, x, par);
	#elif recon_flag==4
	return get_PPM_aveL(geom, x, par);
	#elif recon_flag==5
	return get_PEM_aveL(geom, x, par);
	#elif recon_flag==6
	return get_PPM_aveL(geom, x, par);
	#endif
}

__device__ void get_PRM_parameters(int i, int geom, double* x, double* dx, double* dv, double* a, double* par)
{
	#if recon_flag==0
	get_2nd_VAN_parameters(i, geom, x, dx, dv, a, par);
	#elif recon_flag==1
	get_2nd_MOC_parameters(i, geom, x, dx, dv, a, par);
	#elif recon_flag==2
	get_2nd_parameters(i, geom, x, dx, dv, a, par);
	#elif recon_flag==3
	get_3rd_parameters(i, geom, x, dx, dv, a, par);
	#elif recon_flag==4
	get_3rd_parameters(i, geom, x, dx, dv, a, par);
	adjust_PPM_parameters(i,par);
	#elif recon_flag==5
	get_4th_parameters(i, geom, x, dx, dv, a, par);
	#elif recon_flag==6
	get_4th_parameters(i, geom, x, dx, dv, a, par);
	adjust_PPM_parameters(i,par);
	#endif
	return;
}

__device__ double get_PRM_aveR(int geom, double x, double* par, double lx, double ly)
{
	#if recon_flag==0
	return get_PLM_aveR(geom, x, par);
	#elif recon_flag==1
	return get_PLM_aveR(geom, x, par);
	#elif recon_flag==2
	return get_PEM_aveR(geom, x, par);
	#elif recon_flag==3
	return get_PEM_aveR(geom, x, par);
	#elif recon_flag==4
	return get_PPM_aveR(geom, x, par);
	#elif recon_flag==5
	return get_PEM_aveR(geom, x, par);
	#elif recon_flag==6
	return get_PPM_aveR(geom, x, par);
	#endif
}

__device__ double get_PRM_aveL(int geom, double x, double* par, double lx, double ly)
{
	#if recon_flag==0
	return get_PLM_aveL(geom, x, par);
	#elif recon_flag==1
	return get_PLM_aveL(geom, x, par);
	#elif recon_flag==2
	return get_PEM_aveL(geom, x, par);
	#elif recon_flag==3
	return get_PEM_aveL(geom, x, par);
	#elif recon_flag==4
	return get_PPM_aveL(geom, x, par);
	#elif recon_flag==5
	return get_PEM_aveL(geom, x, par);
	#elif recon_flag==6
	return get_PPM_aveL(geom, x, par);
	#endif
}

__device__ double get_slopeR(int geom, double x1, double x2, double* par)
{
	#if recon_flag==0
	return get_PEM_slopeR(geom, x1, x2, par);
	#elif recon_flag==1
	return get_PLM_slopeR(geom, x1, x2, par);
	#elif recon_flag==2
	return get_PPM_slopeR(geom, x1, x2, par);
	#endif
}

__device__ double get_slopeL(int geom, double x1, double x2, double* par)
{
	#if recon_flag==0
	return get_PEM_slopeL(geom, x1, x2, par);
	#elif recon_flag==1
	return get_PLM_slopeL(geom, x1, x2, par);
	#elif recon_flag==2
	return get_PPM_slopeL(geom, x1, x2, par);
	#endif
}
