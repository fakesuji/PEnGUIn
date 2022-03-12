__device__ double get_ramping_fac(double x)
{
  //f  = cpow(csin(f*hpi), 2.0);//1.0/(1.0 + 999.0*pow(1.0-f, 8.0)) - 0.001;//
  double f = sin(hpi*x);
  f *= f; 
  return f;
}


__global__ void killwave(Grid G, double dt)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int j = threadIdx.y + blockIdx.y*blockDim.y + ypad;
	int k = threadIdx.z + blockIdx.z*blockDim.z + zpad;

	double rad, azi, pol;

	double f, tau;
	#ifdef shear_box
	double d_in  = sc_h*kill_width;
	double d_out = sc_h*kill_width;
	#else
	double d_in  = get_h(xmin,0.0,hpi)*kill_width;
	double d_out = get_h(xmax,0.0,hpi)*kill_width;
	#endif

	double inner = xmin+d_in;
	double outer = xmax-d_out;

	Cell C_tmp;
	Cell I_tmp;

	if (i>=xpad && i<G.xarr-xpad)
	if (j>=ypad && j<G.yarr-ypad)
	if (k>=zpad && k<G.zarr-zpad)
	{
		rad = G.get_xc(i);

		if (rad<inner) 
		{
			azi = G.get_yc(j);
			pol = G.get_zc(k);

			f = (inner-rad)/d_in;
			#if geomx == 0
			tau  = 1.0;
			#else
			tau  = 0.5*pow(rad,1.5);
			#endif

			f  = get_ramping_fac(f);
			f *= dt/tau;

			C_tmp = G.get_cell(i,j,k);
			I_tmp = init_C(rad,azi,pol);

			C_tmp.r += f * ( I_tmp.r - C_tmp.r );
			C_tmp.p += f * ( I_tmp.p - C_tmp.p );
			C_tmp.u += f * ( I_tmp.u - C_tmp.u );
			C_tmp.v += f * ( I_tmp.v - C_tmp.v );
			C_tmp.w += f * ( I_tmp.w - C_tmp.w );

			G.write_cell(i,j,k,C_tmp);
		}

		if (rad>outer)
		{
			azi = G.get_yc(j);
			pol = G.get_zc(k);

			f = (rad-outer)/d_out;
			#if geomx == 0
			tau  = 1.0;
			#else
			tau  = 0.5*pow(rad,1.5);
			#endif

			f  = get_ramping_fac(f);
			f *= dt/tau;

			C_tmp = G.get_cell(i,j,k);
			I_tmp = init_C(rad,azi,pol);

			C_tmp.r += f * ( I_tmp.r - C_tmp.r );
			C_tmp.p += f * ( I_tmp.p - C_tmp.p );
			C_tmp.u += f * ( I_tmp.u - C_tmp.u );
			C_tmp.v += f * ( I_tmp.v - C_tmp.v );
			C_tmp.w += f * ( I_tmp.w - C_tmp.w );

			G.write_cell(i,j,k,C_tmp);
		}
	}
	return;
}
