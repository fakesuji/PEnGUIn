__device__ double get_g_apprx(double x, double k)
{
/*
	double n = 1.0;
	double tmp = -(k+1.0)*y/2.0;
	double sum = 1.0 + tmp;

	do
	{
		tmp *= -(k-n)*y/(n+2.0);
		sum += tmp;
		n += 1.0;
	} while (fabs(tmp/sum)>1.0e-18);

	return sum;
*/
	double y = x-1.0;
	double tmp = 1.0;
	tmp += y*(k+1.0)/2.0;
	tmp += y*y*(k+1.0)*(k-1.0)/6.0;
	tmp += y*y*y*(k+1.0)*(k-1.0)*(k-2.0)/24.0;
	//tmp += y*y*y*y*(k+1.0)*(k-1.0)*(k-2.0)*(k-3.0)/120.0;
	//tmp += y*y*y*y*y*(k+1.0)*(k-1.0)*(k-2.0)*(k-3.0)*(k-4.0)/720.0;

	return tmp;
}

__device__ double div_3rd(int geom, double x0, double x1, double x2, double x3, double a1, double a2, double a3)
{
	double diff, sL, sR, c1, c2, c3, d1, d2, d3;

	d1 = x1-x0;
	d2 = x2-x1;
	d3 = x3-x2;

	double d131, d120, d231, d220, d12, d22, deno;

	if (geom==0)
	{
		d12 = (1.0/2.0)*(x1+x2);
		d22 = (1.0/3.0)*(x1*x1 + x1*x2 + x2*x2);

		d131 = (1.0/2.0)*(d2+d3);
		d120 = (1.0/2.0)*(d1+d2);

		d231 = (1.0/3.0)*(d2+d3)*(x1 + x2 + x3);
		d220 = (1.0/3.0)*(d1+d2)*(x0 + x1 + x2);

		deno = (1.0/6.0)*(d1+d2)*(d2+d3)*(-d1-d2-d3);
	}
	else if (geom==1)
	{
		d12 = (2.0/3.0)*(x1*x1 + x1*x2 + x2*x2)/(x1+x2);
		d22 = (1.0/2.0)*(x1*x1 + x2*x2);

		d131 = (2.0/3.0)*(d2+d3)*(x1*x2 + x2*x3 + x1*x3)/(x1+x2)/(x2+x3);
		d120 = (2.0/3.0)*(d1+d2)*(x0*x1 + x1*x2 + x0*x2)/(x0+x1)/(x1+x2);

		d231 = (1.0/2.0)*(d2+d3)*(x1 + x3);
		d220 = (1.0/2.0)*(d1+d2)*(x0 + x2);
		
		deno = (1.0/3.0)*(d1+d2)*(d2+d3)*(-d1-d2-d3)*(x0*x1*x2 + x1*x2*x3 + x0*x2*x3 + x0*x1*x3)/(x0+x1)/(x1+x2)/(x2+x3);
	}	

	c3 = ((a2-a1)*d131 - (a3-a2)*d120)/deno;
	c2 = ((a3-a2)*d220 - (a2-a1)*d231)/deno;
	c1 = ((a2-a1)*(d231*d12-d131*d22) + (a3-a2)*(d120*d22-d220*d12))/deno;
	
	sL = -(c1 + x1*c2 + x1*x1*c3);
	sR = (c1 + x2*c2 + x2*x2*c3);

	sL = copysign(fmin(fabs(sL),fabs(a2-a1)),a2-a1);
	sR = copysign(fmin(fabs(sR),fabs(a3-a2)),a3-a2);

	if (sR*sL<=0.0)
	{
		diff = 0.0;
	}
	else
	{
		if (sR/sL>2.0) sR = 2.0*sL;
		if (sL/sR>2.0) sL = 2.0*sR;

		diff = (sR+sL)/d2;
    	}

	return diff;
}

//=======================================================================================

__device__ void get_PEM_parameters(int i, int geom, double* x, double* dx, double* dv, double* a, double* par)
{
/*
	double a1, a2, a3, x1, x2, x3, sL, sR, aL, aR;

	x1 = dx[i-1];
	x2 = dx[i];
	x3 = dx[i+1];

	a1 = a[i-1];//*dv[i-1]/x1;
	a2 = a[i];//*dv[i]/x2;
	a3 = a[i+1];//*dv[i+1]/x3;

	//===============================================================================

	sL = (a2-a1)/(x2+x1);
	sR = (a3-a2)/(x3+x2);
    
	aL = x2* ((x3+x2)*sL + x1*sR) / (x1+x2+x3);
	aR = x2* ((x1+x2)*sR + x3*sL) / (x1+x2+x3);

	sL = copysign(fmin(fabs(aL),fabs(a2-a1)),a2-a1);
	sR = copysign(fmin(fabs(aR),fabs(a3-a2)),a3-a2);

	//===============================================================================

	if (sR*sL<=0.0)
	{
		par[0] = 0.0;
		par[1] = a2;
		par[2] = 0.0;
	}
	else
	{
		//if (sR/sL>4.0) sR = 4.0*sL;
		//if (sL/sR>4.0) sL = 4.0*sR;

		par[0] = sL;
		par[1] = a2;
		par[2] = sR;
    	}
	return;
*/

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

	if (geom==0)
	{
		d12 = (1.0/2.0)*(x1+x2);
		d22 = (1.0/3.0)*(x1*x1 + x1*x2 + x2*x2);

		d131 = (1.0/2.0)*(d2+d3);
		d120 = (1.0/2.0)*(d1+d2);

		d231 = (1.0/3.0)*(d2+d3)*(x1 + x2 + x3);
		d220 = (1.0/3.0)*(d1+d2)*(x0 + x1 + x2);

		deno = (1.0/6.0)*(d1+d2)*(d2+d3)*(-d1-d2-d3);
	}
	else if (geom==1)
	{
		d12 = (2.0/3.0)*(x1*x1 + x1*x2 + x2*x2)/(x1+x2);
		d22 = (1.0/2.0)*(x1*x1 + x2*x2);

		d131 = (2.0/3.0)*(d2+d3)*(x1*x2 + x2*x3 + x1*x3)/(x1+x2)/(x2+x3);
		d120 = (2.0/3.0)*(d1+d2)*(x0*x1 + x1*x2 + x0*x2)/(x0+x1)/(x1+x2);

		d231 = (1.0/2.0)*(d2+d3)*(x1 + x3);
		d220 = (1.0/2.0)*(d1+d2)*(x0 + x2);
		
		deno = (1.0/3.0)*(d1+d2)*(d2+d3)*(-d1-d2-d3)*(x0*x1*x2 + x1*x2*x3 + x0*x2*x3 + x0*x1*x3)/(x0+x1)/(x1+x2)/(x2+x3);
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
		S = 0.5*sL*(sqrt(8.0*(a3-a2)/sL+1.0)-1.0);
		if (sR/S>1.0) sR = S;
	}
	else if (sL/sR>1.0)
	{
		S = 0.5*sR*(sqrt(8.0*(a2-a1)/sR+1.0)-1.0);
		if (sL/S>1.0) sL = S;
	}

	par[0] = sL;
	par[1] = a2;
	par[2] = sR;

	return;

/*
	double a1, a2, a3, x1, x2, x3;

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

	x1 = dx[i-1];
	x2 = dx[i];
	x3 = dx[i+1];

	//===============================================================================

	double sL = 0.5*(a2-a1);
	double sR = 0.5*(a3-a2);

	if (sR/sL>6.0) sR = 6.0*sL;
	if (sL/sR>6.0) sL = 6.0*sR;

	par[0] = sL;;
	par[1] = a2;
	par[2] = sR;;

	return;
*/
}

//=======================================================================================

__device__ double lim01(double a)
{
	return fmin(fmax(a,0.0),1.0);
}

__device__ double get_PEM_0(double r0, double rx, double r1, double s0, double aM, double eta)
{
	double dr = r1-r0;
	double x = (rx-r0)/dr;

	double val;
	if (x<1.0e-10*fabs(r1-r0)) val = aM + s0;
	else        val = aM + s0*(1.0-exp(eta*log(x)));

	return val;
}

__device__ double get_PEM_1(double r0, double rx, double r1, double aM, double s1, double eta)
{
	double dr = r1-r0;
	double x = (rx-r0)/dr;
	double y = (r1-rx)/dr;

	double val;
	if          (x<1.0e-10*fabs(r1-r0)) val = aM;
	else if (eta*y>1e-4) val = aM + s1*lim01(x*(1.0-exp(eta*log(x)))/(y*eta));
	else                 val = aM + s1*get_g_apprx(x, eta);

	return val;
}

//=======================================================================================

__device__ double get_PEM_aveR(int geom, double rL, double r0, double rR, double* par)
{
	double sL = par[0];
	double aM = par[1];
	double sR = par[2];
	double val;

	if         (sR==0.0) val = aM;
	else if (sR/sL>=1.0) val = get_PEM_1(rL, r0, rR, aM, sR, sR/sL);
	else                 val = get_PEM_0(rR, r0, rL, sR, aM, sL/sR);

	return val;///get_dv_dr_dev(geom, r0, rR-r0);
}

__device__ double get_PEM_aveL(int geom, double rL, double r0, double rR, double* par)
{
	double sL = par[0];
	double aM = par[1];
	double sR = par[2];
	double val;

	if         (sL==0.0) val = aM;
	else if (sR/sL>=1.0) val = get_PEM_0(rL, r0, rR, -sL, aM, sR/sL);
	else                 val = get_PEM_1(rR, r0, rL, aM, -sL, sL/sR);

	return val;///get_dv_dr_dev(geom, rL, r0-rL);
}

