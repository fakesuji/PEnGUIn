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

__device__ void get_PLM_parameters(int i, int geom, double* r, double* dr, double* dv, double* a, double* par)
{
	double a1, a2, a3, x1, x2, x3, tmp;

	a1 = a[i-1];
	a2 = a[i];
	a3 = a[i+1];
	
	if ((a3-a2)*(a2-a1)<=0.0)
	{
		par[0] = a2;
		par[1] = a2;
		par[2] = a2;
		par[3] = 1.0;
		return;	
	}

	x1 = dr[i-1];
	x2 = dr[i];
	x3 = dr[i+1];

	//===============================================================================

	//tmp = (a3-a1)*x2/(x1+x3+2.0*x2);
  	//tmp = copysign( fmin(fmin(fabs(a2-a1),fabs(a3-a2)),fabs(tmp)) , tmp );
	tmp = (a3-a2)*(a2-a1)/(a3-a1);
	tmp*= x2*(x1+x3+2.0*x2)/(x3+x2)/(x2+x1);

	if (tmp==a2-a1)
		par[0] = a1;
	else
		par[0] = a2-tmp;

	par[1] = a2;

	if (tmp==a3-a2)
		par[2] = a3;
	else
		par[2] = a2+tmp;
    
	return;
}

//=======================================================================================

__device__ double get_PLM_aveR(int geom, double rL, double r0, double rR, double* par)
{
	double aL = par[0];
	double aR = par[2];
	double x = (r0-rL)/(rR-rL);
	double val;

	val = aL + 0.5*(aR-aL)*(1.0+x);

	return val;
}

__device__ double get_PLM_aveL(int geom, double rL, double r0, double rR, double* par)
{
	double aL = par[0];
	double aR = par[2];
	double x = (r0-rL)/(rR-rL);
	double val;

	val = aL + 0.5*(aR-aL)*x;

	return val;
}

__device__ double get_PLM_ave(int i, int geom, double* r, double* dr, double* dv, double* a, double c, double dt)
{
	double par[4];
	double rL, rR;

	get_PLM_parameters(i,geom,r,dr,dv,a,par);
	rL = r[i];
	rR = r[i+1];

	if (c*dt>0.0) return get_PLM_aveR(geom, rL, rR-c*dt, rR, par);
	else          return get_PLM_aveL(geom, rL, rL-c*dt, rR, par);
}

