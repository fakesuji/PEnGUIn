__device__ double get_g_apprx(double x, double k)
{
	double y = x-1.0;
	double tmp = 1.0;
	double val = tmp;

	tmp *= y*(k+1.0)/2.0;
	val += tmp;
	
	tmp *= y*(k-1.0)/3.0;
	val += tmp;

	tmp *= y*(k-2.0)/4.0;
	val += tmp;

//	val += y*(k+1.0)/2.0;
//	val += y*y*(k+1.0)*(k-1.0)/6.0;
//	val += y*y*y*(k+1.0)*(k-1.0)*(k-2.0)/24.0;

	return val;
}

//=======================================================================================

__device__ void get_PEM_parameters(int i, int geom, double* x, double* dx, double* dv, double* a, double* par)
{
	double a1, a2, a3, sL, sR, c1, c2, c3, x0, x1, x2, x3, d1, d2, d3;

	a1 = a[i-1];
	a2 = a[i];
	a3 = a[i+1];

	x0 = x[i-1];
	x1 = x[i];
	x2 = x[i+1];
	x3 = x[i+2];

	d1 = dx[i-1];
	d2 = dx[i];
	d3 = dx[i+1];

	double d131, d120, d231, d220, d12, d22, deno;

	if (geom==1)
	{
		d12 = (2.0/3.0)*(x1*x1 + x1*x2 + x2*x2)/(x1+x2);
		d22 = (1.0/2.0)*(x1*x1 + x2*x2);

		d131 = (2.0/3.0)*(d2+d3)*(x1*x2 + x2*x3 + x1*x3)/(x1+x2)/(x2+x3);
		d120 = (2.0/3.0)*(d1+d2)*(x0*x1 + x1*x2 + x0*x2)/(x0+x1)/(x1+x2);

		d231 = (1.0/2.0)*(d2+d3)*(x1 + x3);
		d220 = (1.0/2.0)*(d1+d2)*(x0 + x2);
		
		deno = (1.0/3.0)*(d1+d2)*(d2+d3)*(-d1-d2-d3)*(x0*x1*x2 + x1*x2*x3 + x0*x2*x3 + x0*x1*x3)/(x0+x1)/(x1+x2)/(x2+x3);
	}	
	else
	{
		d12 = (1.0/2.0)*(x1+x2);
		d22 = (1.0/3.0)*(x1*x1 + x1*x2 + x2*x2);

		d131 = (1.0/2.0)*(d2+d3);
		d120 = (1.0/2.0)*(d1+d2);

		d231 = (1.0/3.0)*(d2+d3)*(x1 + x2 + x3);
		d220 = (1.0/3.0)*(d1+d2)*(x0 + x1 + x2);

		deno = (1.0/6.0)*(d1+d2)*(d2+d3)*(-d1-d2-d3);
	}

	c3 = ((a2-a1)*d131 - (a3-a2)*d120)/deno;
	c2 = ((a3-a2)*d220 - (a2-a1)*d231)/deno;
	c1 = ((a2-a1)*(d231*d12-d131*d22) + (a3-a2)*(d120*d22-d220*d12))/deno;
	
	sL = -(c1 + x1*c2 + x1*x1*c3);
	sR = (c1 + x2*c2 + x2*x2*c3);

	sL = copysign(fmin(fabs(sL),fabs(a2-a1)),a2-a1);
	sR = copysign(fmin(fabs(sR),fabs(a3-a2)),a3-a2);

	double S;
	
	if (sR*sL<=0.0)
	{
		par[0] = 0.0;
		par[1] = a2;
		par[2] = 0.0;
		return;
	}

	if (sR/sL>1.0)
	{
		if (geom==1) S = (a3*(0.5*(x3+x2))-a1*(0.5*(x2+x1)))/x2;
		else         S = (a3-a2);

		S = 0.5*sL*(sqrt(4.0*S/sL+1.0)-1.0);
		if (sR/S>1.0 && S/sL>1.0) sR = S;
	}
	else if (sL/sR>1.0)
	{
		if (geom==1) S = (a2*(0.5*(x2+x1))-a1*(0.5*(x1+x0)))/x1;
		else         S = (a2-a1);

		S = 0.5*sR*(sqrt(4.0*S/sR+1.0)-1.0);
		if (sL/S>1.0 && S/sR>1.0) sL = S;
	}

	//if      (sR/sL>6.0) sR = 6.0*sL;
	//else if (sL/sR>6.0) sL = 6.0*sR;


	par[0] = sL;
	par[1] = a2;
	par[2] = sR;

	return;
}

//=======================================================================================

__device__ double get_PEM_0(double x, double s0, double aM, double eta, double lx)
{
	double val;
	val = aM + s0*(1.0-__expf(eta*lx));

	return val;
}

__device__ double get_PEM_1(double x, double aM, double s1, double eta, double lx)
{
	double val;
	double y=1.0-x;
	if (eta*y>1e-4) val = aM + s1*lim01(x*(1.0-__expf(eta*lx))/(y*eta));
	else            val = aM + s1*get_g_apprx(x, eta);

	return val;
}

//=======================================================================================

__device__ double get_PEM_aveR(int geom, double x, double* par, double lx, double ly)
{
	double sL = par[0];
	double aM = par[1];
	double sR = par[2];
	double val;

	if      (sR==0.0) val = aM;
	else if (sR/sL>=1.0) 
	{
		if      (x>=1.0) val = aM + sR;
		else if (x<=0.0) val = aM;
		else             val = get_PEM_1(x, aM, sR, sR/sL, lx);
	}
	else
	{
		x = 1.0-x;
		if      (x>=1.0) val = aM;
		else if (x<=0.0) val = aM + sR;
		else             val = get_PEM_0(x, sR, aM, sL/sR, ly);
	}

	return val;
}

__device__ double get_PEM_aveL(int geom, double x, double* par, double lx, double ly)
{
	double sL = par[0];
	double aM = par[1];
	double sR = par[2];
	double val;

	if      (sL==0.0) val = aM;
	else if (sR/sL>=1.0) 
	{
		if      (x>=1.0) val = aM;
		else if (x<=0.0) val = aM - sL;
		else             val = get_PEM_0(x, -sL, aM, sR/sL, lx);
	}
	else
	{
		x = 1.0-x;
		if      (x>=1.0) val = aM - sL;
		else if (x<=0.0) val = aM;
		else             val = get_PEM_1(x, aM, -sL, sL/sR, ly);
	}

	return val;
}

__device__ double get_PEM_slope(int geom, double x, double* par)
{
	double sL = par[0];
	double sR = par[2];
	double eta, val;

	if      (sL==0.0) val = 0.0;
	else if (sR/sL>=1.0) 
	{
		eta = sR/sL;
		if (x>0.0) val = eta*(sR+sL)*__expf((eta-1.0)*__logf(x));
		else val = 0.0;
	}
	else
	{
		x = 1.0-x;
		eta = sL/sR;
		if (x>0.0) val = eta*(sR+sL)*__expf((eta-1.0)*__logf(x));
		else val = 0.0;
	}

	return val;
}
