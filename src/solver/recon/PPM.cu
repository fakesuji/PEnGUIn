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
				ddp = 10.0*fmax(0.0,fabs(dp1f/dp2f - 0.5) - 0.25);
				ddp = Mf * ddp;
				f = ddp;

				Mb  = 1.0;//fmax(0.0,10.0*Mb);	
				ddp = 10.0*fmax(0.0,fabs(dp1b/dp2b - 0.5) - 0.25);
				ddp = Mf * ddp;
				f = fmax(f, ddp);

				f = fmin(0.5, f * fmax(fabs(dp1f),fabs(dp1b))/p0 );
			}
		}
	}

	flat[is] = f;

  return;
}

__device__ void adjust_PPM_parameters(int i, double* par)
{
	double sL = par[0];
	double sR = par[2];

	#if recon_flag == 6
	extern __shared__ double share[];
	double* tmp = &share[0];
	int imax = blockDim.x;
	int is = i + imax*threadIdx.y;
	double flat = tmp[is];

	sR *= (1.0 - flat);
	sL *= (1.0 - flat);
	#endif
	
	double t3 = sR+sL;
	double t2 = t3*t3;
	double t1 = t3*3.0*(sL-sR);
	
	if      (t2 <  t1) sL = 2.0*sR;
	else if (t2 < -t1) sR = 2.0*sL;

	par[0] = sL;
	par[2] = sR;

	return;
}

//=======================================================================================

__device__ double get_PPM_aveR(int geom, double x, double* par)
{
	double a6 = 3.0*(par[0]-par[2]);
	double aR = par[1]+par[2];
	double S = par[2]+par[0];

	x = 1.0 - x;

	return aR - (x/2.0) * ( S - (1.0 - x/1.5) * a6);
}

__device__ double get_PPM_aveL(int geom, double x, double* par)
{
	double a6 = 3.0*(par[0]-par[2]);
	double aL = par[1]-par[0];
	double S = par[2]+par[0];

	return aL + (x/2.0) * ( S + (1.0 - x/1.5) * a6);
}
