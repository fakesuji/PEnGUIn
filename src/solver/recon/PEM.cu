
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
	tmp += y*y*y*y*(k+1.0)*(k-1.0)*(k-2.0)*(k-3.0)/120.0;
	tmp += y*y*y*y*y*(k+1.0)*(k-1.0)*(k-2.0)*(k-3.0)*(k-4.0)/720.0;

	return tmp;
}

//=======================================================================================

__device__ void get_PEM_parameters(int i, int geom, double* x, double* dx, double* dv, double* a, double* par)
{
	double a1, a2, a3, x1, x2, x3, sL, sR, aL, aR;

	x1 = dx[i-1];
	x2 = dx[i];
	x3 = dx[i+1];

	a1 = a[i-1]*dv[i-1]/x1;
	a2 = a[i]*dv[i]/x2;
	a3 = a[i+1]*dv[i+1]/x3;

	//===============================================================================

	sL = (a2-a1)/(x2+x1);
	sR = (a3-a2)/(x3+x2);
    
	aL = x2* ((x3+x2)*sL + x1*sR) / (x1+x2+x3);
	aR = x2* ((x1+x2)*sR + x3*sL) / (x1+x2+x3);

	sL = copysign(fmin(fabs(aL),fabs(a2-a1)),a2-a1);
	sR = copysign(fmin(fabs(aR),fabs(a3-a2)),a3-a2);

	aL = a2 - sL;
	aR = a2 + sR;

	//===============================================================================

	if (sR*sL<=0.0)
	{
		par[0] = 0.0;
		par[1] = a2;
		par[2] = 0.0;
		par[3] = 1.0;
	}
	else
	{
		//if (aR/aL>2.0) aR = 2.0*aL;
		//if (aL/aR>2.0) aL = 2.0*aR;
		par[0] = sL;
		par[1] = a2;
		par[2] = sR;
		par[3] = sR/sL;
    	}
	return;
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
	val  = aM + s0*(1.0-exp(eta*log(x)));

	return val;
}

__device__ double get_PEM_1(double r0, double rx, double r1, double aM, double s1, double eta)
{
	double dr = r1-r0;
	double x = (rx-r0)/dr;
	double y = (r1-rx)/dr;

	double val;
	if (eta*y>1e-3) val = aM + s1*lim01(x*(1.0-exp(eta*log(x)))/(y*eta));
	else            val = aM + s1*get_g_apprx(x, eta);

	return val;
}

//=======================================================================================

__device__ double get_PEM_aveR(int geom, double rL, double r0, double rR, double* par)
{
	double aM = par[1];
	double sR = par[2];
	double eta = par[3];
	double val;

	if (eta>=1.0) val = get_PEM_1(rL, r0, rR, aM, sR, eta);
	else          val = get_PEM_0(rR, r0, rL, sR, aM, 1.0/eta);

	return val/get_dv_dr(geom, r0, rR-r0);
}

__device__ double get_PEM_aveL(int geom, double rL, double r0, double rR, double* par)
{
	double aM = par[1];
	double sL = par[0];
	double eta = par[3];
	double val;

	if (eta>=1.0) val = get_PEM_0(rL, r0, rR, -sL, aM, eta);
	else          val = get_PEM_1(rR, r0, rL, aM, -sL, 1.0/eta);

	return val/get_dv_dr(geom, rL, r0-rL);
}

__device__ double get_PEM_ave(int i, int geom, double* r, double* dr, double* dv, double* a, double c, double dt)
{
	double par[4];
	double rL, rR;

	get_PEM_parameters(i,geom,r,dr,dv,a,par);
	rL = r[i];
	rR = r[i+1];

	if (c*dt>0.0) return get_PEM_aveR(geom, rL, rR-c*dt, rR, par);
	else          return get_PEM_aveL(geom, rL, rL-c*dt, rR, par);
}

