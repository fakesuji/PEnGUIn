__device__ double PLM_div(double x1, double x2, double x3, double a1, double a2, double a3)
{
	double tmp;

	tmp = (a3-a1)*x2/(x1+x3+2.0*x2);
  	tmp = copysign( fmin(fmin(fabs(a2-a1),fabs(a3-a2)),fabs(tmp)) , tmp );

	if ((a3-a2)*(a2-a1)<=0.0)
	{
		return 0.0;
	}
	else
	{
		return 2.0*tmp;
    	}
}

//=======================================================================================

__device__ void get_MOC_parameters(int i, int geom, double* r, double* dr, double* dv, double* a, double* par)
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

__device__ void get_VAN_parameters(int i, int geom, double* r, double* dr, double* dv, double* a, double* par)
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

__device__ double get_PLM_aveR(int geom, double x, double* par)
{
	double sR = par[2];
	double val;

	val = par[1] + sR*x;

	return val;
}

__device__ double get_PLM_aveL(int geom, double x, double* par)
{
	double sL = par[0];
	double val;

	val = par[1] - sL*(1.0-x);

	return val;
}

__device__ double get_PLM_slope(int geom, double x, double* par)
{
	return par[2]+par[0];
}
