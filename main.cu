#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <chrono>

using namespace std;

#include "parameters.h"
#include "structs.h"
#include "util.h"
#include "output.h"
#include "timestep.h"
#include "geom.h"
#include "init.h"
#include "orbit.h"
#include "planet.h"
#include "solver.h"

bool sanity_check()
{
	bool sane = true;

	if (xres%x_xdiv!=0) sane = false;
	if (xres%y_xdiv!=0) sane = false;
	if (xres%z_xdiv!=0) sane = false;

	#if ndim>1
	if (yres%x_ydiv!=0) sane = false;
	if (yres%y_ydiv!=0) sane = false;
	if (yres%z_ydiv!=0) sane = false;
	#endif

	#if ndim>2
	if (zres%x_zdiv!=0) sane = false;
	if (zres%y_zdiv!=0) sane = false;
	if (zres%z_zdiv!=0) sane = false;
	#endif

	return sane;
}

void cpy_grid_DevicetoHost(Grid* hst, Grid* dev)
{
	
	//for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);
		cudaMemcpyAsync( hst[n].C, dev[n].C, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(Cell), cudaMemcpyDeviceToHost, dev[n].stream );

		#ifdef dust_flag
		cudaMemcpyAsync( hst[n].CD, dev[n].CD, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(Dust), cudaMemcpyDeviceToHost, dev[n].stream );
		#endif

		cudaMemcpyAsync( hst[n].orb_shf, dev[n].orb_shf, dev[n].xarr*dev[n].zarr*sizeof(int), cudaMemcpyDeviceToHost, dev[n].stream );
		cudaMemcpyAsync( hst[n].planets, dev[n].planets, n_planet*sizeof(body), cudaMemcpyDeviceToHost, dev[n].stream );
	}
	for (int n=0; n<ndev; n++) 
	{
		cudaSetDevice(n);
		cudaDeviceSynchronize();
		//cudaStreamSynchronize(dev[n].stream);
	}
	return;
}

void cpy_grid_HosttoDevice(Grid* hst, Grid* dev)
{
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);
		cudaMemcpyAsync( dev[n].xa, hst[n].xa, (xarr+1)*sizeof(double), cudaMemcpyHostToDevice, dev[n].stream );
		cudaMemcpyAsync( dev[n].xv, hst[n].xv, xarr*sizeof(double), cudaMemcpyHostToDevice, dev[n].stream );

		cudaMemcpyAsync( dev[n].ya, hst[n].ya, (yarr+1)*sizeof(double), cudaMemcpyHostToDevice, dev[n].stream );
		cudaMemcpyAsync( dev[n].yv, hst[n].yv, yarr*sizeof(double), cudaMemcpyHostToDevice, dev[n].stream );

		cudaMemcpyAsync( dev[n].za, hst[n].za, (zarr+1)*sizeof(double), cudaMemcpyHostToDevice, dev[n].stream );
		cudaMemcpyAsync( dev[n].zv, hst[n].zv, zarr*sizeof(double), cudaMemcpyHostToDevice, dev[n].stream );

		cudaMemcpyAsync( dev[n].C, hst[n].C, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(Cell), cudaMemcpyHostToDevice, dev[n].stream );

		#ifdef dust_flag
		cudaMemcpyAsync( dev[n].CD, hst[n].CD, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(Dust), cudaMemcpyHostToDevice, dev[n].stream );
		#endif

		cudaMemcpyAsync( dev[n].planets, hst[n].planets, n_planet*sizeof(body), cudaMemcpyHostToDevice, dev[n].stream );
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
	return;
}

void init_grid_dimensions(Grid* G)
{
	for (int n=0; n<ndev; n++)
	{
		G[n].xres = xres/ndev;
		G[n].xarr = G[n].xres + 2*xpad;
		G[n].xbgn = G[n].xres*n;

		G[n].yres = yres;
		G[n].yarr = G[n].yres + 2*ypad;
		G[n].ybgn = 0;

		G[n].zres = zres;
		G[n].zarr = G[n].zres + 2*zpad;
		G[n].zbgn = 0;
	}
	return;
}
/*
void save_check_point(ofstream &check_point, string fname, int &sstep, double cur_time, Grid* hst)
{
	open_binary_file(check_point,fname);
	write_check_point(check_point, cur_time, hst);
	close_output_file(check_point);
	return;
}
*/
int main(int narg, char *args[])
{
	if (!sanity_check())
	{
		printf("incorrect grid setting.\n");
		return 1;
	}

	double sta_time = 0.0;
	double cur_time = sta_time;
	double sav_time = 0.0;
	double dt = end_time-cur_time;
	int nstep = 0;
	int sstep = 0;
	int tmp;
	ofstream check_point;

	string label=create_label();
	string path = "/scratch/fung/";
	string fname;

	Grid* hst;
	Grid* dev;

	cudaMallocHost( (void**)&hst, ndev*sizeof(Grid) );
	cudaMallocHost( (void**)&dev, ndev*sizeof(Grid) );

	init_grid_dimensions(hst);
	init_grid_dimensions(dev);

	///////////////////////////////////////////////////////////
	
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		cudaMallocHost( (void**)&hst[n].xa, (xarr+1)*sizeof(double) );
		cudaMallocHost( (void**)&hst[n].xv, xarr*sizeof(double) );

		cudaMallocHost( (void**)&hst[n].ya, (yarr+1)*sizeof(double) );
		cudaMallocHost( (void**)&hst[n].yv, yarr*sizeof(double) );

		cudaMallocHost( (void**)&hst[n].za, (zarr+1)*sizeof(double) );
		cudaMallocHost( (void**)&hst[n].zv, zarr*sizeof(double) );
	
		cudaMallocHost( (void**)&hst[n].C, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(Cell) );

		#ifdef dust_flag
		cudaMallocHost( (void**)&hst[n].CD, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(Dust) );
		#endif

		cudaMallocHost( (void**)&hst[n].orb_shf, dev[n].xarr*dev[n].zarr*sizeof(int) );

		cudaMallocHost( (void**)&hst[n].planets, n_planet*sizeof(body) );

		cudaMallocHost( (void**)&hst[n].dt, sizeof(double) );


		cudaStreamCreate(&dev[n].stream);

		cudaMalloc( (void**)&dev[n].xa, (xarr+1)*sizeof(double) );
		cudaMalloc( (void**)&dev[n].xv, xarr*sizeof(double) );

		cudaMalloc( (void**)&dev[n].ya, (yarr+1)*sizeof(double) );
		cudaMalloc( (void**)&dev[n].yv, yarr*sizeof(double) );

		cudaMalloc( (void**)&dev[n].za, (zarr+1)*sizeof(double) );
		cudaMalloc( (void**)&dev[n].zv, zarr*sizeof(double) );

		cudaMalloc( (void**)&dev[n].dt, sizeof(double) );
		cudaMalloc( (void**)&dev[n].Buff, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(double) );

		cudaMalloc( (void**)&dev[n].C, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(Cell) );
		cudaMalloc( (void**)&dev[n].T, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(Cell) );
		//cudaMalloc( (void**)&dev[n].F, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(Cell) );

		cudaMalloc( (void**)&dev[n].fx, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(double) );
		cudaMalloc( (void**)&dev[n].fy, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(double) );
		cudaMalloc( (void**)&dev[n].fz, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(double) );

		#ifdef dust_flag
		cudaMalloc( (void**)&dev[n].CD, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(Dust) );
		cudaMalloc( (void**)&dev[n].TD, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(Dust) );

		cudaMalloc( (void**)&dev[n].fx_d, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(double) );
		cudaMalloc( (void**)&dev[n].fy_d, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(double) );
		cudaMalloc( (void**)&dev[n].fz_d, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(double) );
		#endif

		#ifdef visc_flag
		#if ndim==1
		cudaMalloc( (void**)&dev[n].vis_tensor, ndim*dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(double) );
		#elif ndim==2
		cudaMalloc( (void**)&dev[n].vis_tensor, ndim*dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(double) );
		#elif ndim==3
		cudaMalloc( (void**)&dev[n].vis_tensor, ndim*dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(double) );
		#endif
		#endif

		#ifdef RadPres_flag
		cudaMalloc( (void**)&dev[n].ext, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(double) );
		#endif

		cudaMalloc( (void**)&dev[n].BuffL, xpad*dev[n].yarr*dev[n].zarr*sizeof(Cell) );
		cudaMalloc( (void**)&dev[n].BuffR, xpad*dev[n].yarr*dev[n].zarr*sizeof(Cell) );

		#ifdef dust_flag
		cudaMalloc( (void**)&dev[n].BuffLD, xpad*dev[n].yarr*dev[n].zarr*sizeof(Dust) );
		cudaMalloc( (void**)&dev[n].BuffRD, xpad*dev[n].yarr*dev[n].zarr*sizeof(Dust) );
		#endif

		cudaMalloc( (void**)&dev[n].orb_rot, dev[n].xarr*dev[n].zarr*sizeof(double) );
		cudaMalloc( (void**)&dev[n].orb_res, dev[n].xarr*dev[n].zarr*sizeof(double) );
		cudaMalloc( (void**)&dev[n].orb_shf, dev[n].xarr*dev[n].zarr*sizeof(int) );
		cudaMalloc( (void**)&dev[n].orb_inc, dev[n].xarr*dev[n].zarr*sizeof(int) );

		cudaMalloc( (void**)&dev[n].shfL, xpad*dev[n].zres*sizeof(int) );		
		cudaMalloc( (void**)&dev[n].shfR, xpad*dev[n].zres*sizeof(int) );		

		cudaMalloc( (void**)&dev[n].planets, n_planet*sizeof(body) );
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);

	///////////////////////////////////////////////////////////
	
	init(hst);
	
	if (narg==3) 
	{
		fname = args[1];
		sstep = atoi(args[2]);
		if ( atoi(args[2]) < 0 ) sstep = -sstep;
		sta_time = load_grid(hst,fname);
		cur_time = sta_time;
		printf("restarting from t = %f.\n",cur_time);
		sstep++;
	}
	else
	{
		fname = path+"binary_"+label+"_"+frame_num(sstep);
		check_point.open(fname.c_str(), ios::out | ios::binary);
		write_check_point(check_point, cur_time, hst);
		check_point.close();
		#ifdef dump_flag
		sstep++;
		#endif
	}

	///////////////////////////////////////////////////////////

	init_planet(hst,cur_time);

	cpy_grid_HosttoDevice(hst, dev);

	init_OrbAdv(dev);

	///////////////////////////////////////////////////////////

	auto clock = chrono::high_resolution_clock::now();
	auto begin = chrono::high_resolution_clock::now();
	double duration;
	double speed;
	double old_dt;

	bool print = false;
	bool savep = false;
	bool contu = true;

	while (contu)
	{
		dt = global_dt(hst, dev, dt);
		if (cur_time+dt > end_time)
		{
			dt = end_time-cur_time;
			print = true;
			#ifdef dump_flag
			savep = true;
			#endif
			contu = false;
		}
		#ifdef dump_flag
		else if (sav_time+dt > sav_interval)
		{
			old_dt = dt;
			dt = sav_interval-sav_time;
			savep = true;
		}
		#endif

		////////////////////////////////////////////////////////////////////////////////////

		solve(dev, cur_time, dt);

		cur_time += dt;
		sav_time += dt;
		nstep++;
		#ifndef silence_flag
		if (nstep%prt_interval==0) print = true;
		#endif
		if (nstep>=max_step) {contu = false; print = true; savep = true;}

		////////////////////////////////////////////////////////////////////////////////////

		if (print)
		{
			clock = chrono::high_resolution_clock::now();
			tmp = chrono::duration_cast<std::chrono::milliseconds>(clock-begin).count();
			duration = double(tmp)/1000.0;
			speed = duration/(double)nstep;
			duration /= 60.0;

			printf("Step %i:\n", nstep);
			printf("  Current time %f | End time %f | Time step %e\n", cur_time, end_time, dt);
			printf("  ETD %f / %f minutes\n", duration, duration*(end_time-sta_time)/(cur_time-sta_time));
			printf("  %f seconds per step\n", speed);
			printf("  %f Mcells per second \n\n", xarr*yarr*zarr/speed/1000000.0);

			print = false;
		}

		if (savep)
		{
			printf("Check point!\n");	
			//printf("  Current time %f | End time %f | Time step %e\n", cur_time, end_time, dt);

			cpy_grid_DevicetoHost(hst, dev);
			fname = path+"binary_"+label+"_"+frame_num(sstep);

			check_point.open(fname.c_str(), ios::out | ios::binary);
			write_check_point(check_point, cur_time, hst);
			check_point.close();

			printf("  %s is saved. \n\n", fname.c_str());

			sstep++;
			sav_time = 0.0;
			savep = false;
			dt = old_dt;
		}
	}

	///////////////////////////////////////////////////////////
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		cudaFreeHost(hst[n].xa);
		cudaFreeHost(hst[n].xv);

		cudaFreeHost(hst[n].ya);
		cudaFreeHost(hst[n].yv);

		cudaFreeHost(hst[n].za);
		cudaFreeHost(hst[n].zv);

		cudaFreeHost(hst[n].C);

		#ifdef dust_flag
		cudaFreeHost(hst[n].CD);
		#endif

		cudaFreeHost(hst[n].orb_shf);
		cudaFreeHost(hst[n].planets);
		cudaFreeHost(hst[n].dt);


		cudaFree(dev[n].xa);
		cudaFree(dev[n].xv);

		cudaFree(dev[n].ya);
		cudaFree(dev[n].yv);

		cudaFree(dev[n].za);
		cudaFree(dev[n].zv);

		cudaFree(dev[n].dt);
		cudaFree(dev[n].Buff);

		cudaFree(dev[n].C);
		cudaFree(dev[n].T);
		//cudaFree(dev[n].F);

		cudaFree(dev[n].fx);
		cudaFree(dev[n].fy);
		cudaFree(dev[n].fz);

		#ifdef dust_flag
		cudaFree(dev[n].CD);
		cudaFree(dev[n].TD);

		cudaFree(dev[n].fx_d);
		cudaFree(dev[n].fy_d);
		cudaFree(dev[n].fz_d);
		#endif

		#ifdef visc_flag
		cudaFree(dev[n].vis_tensor);
		#endif

		#ifdef RadPres_flag
		cudaFree(dev[n].ext);
		#endif

		cudaFree(dev[n].BuffL);
		cudaFree(dev[n].BuffR);

		#ifdef dust_flag
		cudaFree(dev[n].BuffLD);
		cudaFree(dev[n].BuffRD);
		#endif

		cudaFree(dev[n].orb_rot);
		cudaFree(dev[n].orb_res);
		cudaFree(dev[n].orb_shf);
		cudaFree(dev[n].orb_inc);

		cudaFree(dev[n].shfL);
		cudaFree(dev[n].shfR);

		cudaFree(dev[n].planets);
	}

	cudaFreeHost(hst);
	cudaFreeHost(dev);

	///////////////////////////////////////////////////////////

	return 0;
}
