//===============================================================================
//use the numbers below for more aggresive flattening
//#define flat_omega1 0.0     // 0.5    //0.75
//#define flat_omega2 3.0     // 10.0   //5.0
//#define flat_epsilon 0.0    // 1.0    //0.33

__device__ void flatten(double *r, double *p, double *u)
{
	extern __shared__ double share[];
	double* flat = &share[0];
	int imax = blockDim.x;
	int i = threadIdx.x;
	int is = i + imax*threadIdx.y;

 	double f = 0.0;
	double dp1f, dp2f, dp1b, dp2b, ddp, Mf, Mb;
	double pb, p0, pf;
	double cs;
	bool shock = false;

	if (i>=2 && i<imax-2)
	{
		pb = p[i-1];
		p0 = p[i];
		pf = p[i+1];
 
		cs = sqrt(gam*p0/r[i]);

		dp1f = pf - p0;
		dp2f = p[i+2] - pb;
  
		dp1b = p0 - pb;
		dp2b = pf - p[i-2];

		Mf = (u[i] - u[i+1])/cs;
		Mb = (u[i-1] - u[i])/cs;

		if ( Mf>0.0 || Mb>0.0 ) shock = true;

		if ( shock )
		{
			if (dp2f!=0.0 && dp2b!=0.0)
			{
				Mf  = 1.0;//fmax(0.0,10.0*Mf);	
				ddp = 10.0*fmax(0.0,fabs(dp1f/dp2f - 0.5) - 0.125);
				ddp = Mf * ddp;
				f = ddp;

				Mb  = 1.0;//fmax(0.0,10.0*Mb);	
				ddp = 10.0*fmax(0.0,fabs(dp1b/dp2b - 0.5) - 0.125);
				ddp = Mf * ddp;
				f = fmax(f, ddp);

				f = fmin(0.5, f * fmax(fabs(dp1f),fabs(dp1b))/p0 );
			}
		}
	}

	flat[is] = f;

  return;
}

__device__ void get_PPM_parameters(int i, int geom, double* r, double* dr, double* dv, double* a, double* par)
{
	extern __shared__ double share[];
	double* tmp = &share[0];
	int imax = blockDim.x;
	int is = i + imax*threadIdx.y;
	double flat = tmp[is];

	double d1,d2,d3,d4,d5,a2,a3,a4;
  	double a21, a32, a43, a54;
 	double daL, daC, daR;
 	double sR, sL;
  	double t1, t2, t3;

	d1 = dr[i-2];
	d2 = dr[i-1];
	d3 = dr[i];
	d4 = dr[i+1];
	d5 = dr[i+2];

	a2 = a[i-1];
	a3 = a[i];
	a4 = a[i+1];

	a21 = a2-a[i-2];
	a32 = a3-a2;
	a43 = a4-a3;
	a54 = a[i+2]-a4;

	/////////////////////////////////////////////////////////////////////////////
	
	daL = (d2/(d1+d2+d3)) * (a32*(d1+d1+d2)/(d3+d2) + a21*(d3+d3+d2)/(d1+d2));
	daC = (d3/(d2+d3+d4)) * (a43*(d2+d2+d3)/(d4+d3) + a32*(d4+d4+d3)/(d2+d3));
	daR = (d4/(d3+d4+d5)) * (a54*(d3+d3+d4)/(d5+d4) + a43*(d5+d5+d4)/(d3+d4));

	if (a32*a21<=0.0) daL = 0.0;
	else daL = copysign( fmin( 2.0*fmin(fabs(a21), fabs(a32)), fabs(daL)) , daL );

	if (a43*a32<=0.0) daC = 0.0;
	else daC = copysign( fmin( 2.0*fmin(fabs(a32), fabs(a43)), fabs(daC)) , daC );

	if (a54*a43<=0.0) daR = 0.0;
	else daR = copysign( fmin( 2.0*fmin(fabs(a43), fabs(a54)), fabs(daR)) , daR );

	/////////////////////////////////////////////////////////////////////////////

	t1 = d2+d3+d4+d5;
	t2 = (d2+d3)/(d3+d3+d4)/t1;
	t3 = (d4+d5)/(d3+d4+d4)/t1;
	sR = a43*(d3/(d3+d4))*(1.0 + 2.0*d4*(t2-t3)) - daR*d3*t2 + daC*d4*t3;

	t1 = d4+d3+d2+d1;
	t2 = (d4+d3)/(d3+d3+d2)/t1;
	t3 = (d2+d1)/(d3+d2+d2)/t1;
	sL = a32*(d3/(d3+d2))*(1.0 + 2.0*d2*(t2-t3)) - daL*d3*t2 + daC*d2*t3;

	/////////////////////////////////////////////////////////////////////////////

	sR *= 1.0 - flat;
	sL *= 1.0 - flat;

	t3 = sR+sL;
	t2 = t3*t3;
	t1 = t3*(-3.0*(sR-sL));

	if (sR*sL<=0.0)
	{
		sR = 0.0;
		sL = 0.0;
	}
	else if (t2 <  t1) sL = 2.0*sR;
	else if (t2 < -t1) sR = 2.0*sL;

	par[0] = sL;
	par[1] = a3;
	par[2] = sR;
    
	return;
}

//=======================================================================================

__device__ double get_PPM_aveR(int geom, double x, double* par)
{
	double sL = par[0];
	double am = par[1];
	double sR = par[2];
	double a6 = -3.0*(sR-sL);
	x = 1.0 - x;

	return am+sR - (x/2.0) * ( sR+sL - (1.0 - x/1.5) * a6);
}

__device__ double get_PPM_aveL(int geom, double x, double* par)
{
	double sL = par[0];
	double am = par[1];
	double sR = par[2];
	double a6 = -3.0*(sR-sL);

	return am-sL + (x/2.0) * ( sR+sL + (1.0 - x/1.5) * a6);
}

__device__ double get_PPM_slopeR(int geom, double x1, double x2, double* par)
{
	double sL = par[0];
	double sR = par[2];
	double a6 = -3.0*(sR-sL);
	x2 = 1.0 - x2;
	x1 = 1.0 - x1;

	return (sR+sL-a6) + a6*(x1+x2)/1.5;
}

__device__ double get_PPM_slopeL(int geom, double x1, double x2, double* par)
{
	double sL = par[0];
	double sR = par[2];
	double a6 = -3.0*(sR-sL);

	return (sR+sL+a6) - a6*(x2+x1)/1.5;
}
