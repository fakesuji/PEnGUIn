#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <chrono>

using namespace std;

#include "structs.h"
#include "util.h"
#include "parameters.h"
#include "output.h"
#include "timestep.h"
#include "geom.h"
#include "init.h"
#include "orbit.h"
#include "solver.h"

void cpy_grid_DevicetoHost(Grid* hst, Grid* dev)
{
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);
		cudaMemcpyAsync( hst[n].C, dev[n].C, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(Cell), cudaMemcpyDeviceToHost, dev[n].stream );
		cudaMemcpyAsync( hst[n].orb_shf, dev[n].orb_shf, dev[n].xarr*dev[n].zarr*sizeof(int), cudaMemcpyDeviceToHost, dev[n].stream );
		cudaMemcpyAsync( hst[n].planets, dev[n].planets, n_planet*sizeof(body), cudaMemcpyDeviceToHost, dev[n].stream );
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
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

		cudaMemcpyAsync( dev[n].planets, hst[n].planets, n_planet*sizeof(body), cudaMemcpyHostToDevice, dev[n].stream );
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
	return;
}

void mem_allocation(Grid* hst, Grid* dev)
{
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
		cudaMallocHost( (void**)&hst[n].orb_shf, dev[n].xarr*dev[n].zarr*sizeof(int) );

		cudaMallocHost( (void**)&hst[n].planets, n_planet*sizeof(body) );


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
		cudaMalloc( (void**)&dev[n].F, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(Cell) );

		cudaMalloc( (void**)&dev[n].fx, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(double) );
		cudaMalloc( (void**)&dev[n].fy, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(double) );
		cudaMalloc( (void**)&dev[n].fz, dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(double) );

		#ifdef visc_flag
		#if ndim==1
		cudaMalloc( (void**)&dev[n].vis_tensor, 1*dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(double) );
		#elif ndim==2
		cudaMalloc( (void**)&dev[n].vis_tensor, 3*dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(double) );
		#elif ndim==3
		cudaMalloc( (void**)&dev[n].vis_tensor, 6*dev[n].xarr*dev[n].yarr*dev[n].zarr*sizeof(double) );
		#endif
		#endif

		cudaMalloc( (void**)&dev[n].BuffL, xpad*dev[n].yres*dev[n].zres*sizeof(Cell) );
		cudaMalloc( (void**)&dev[n].BuffR, xpad*dev[n].yres*dev[n].zres*sizeof(Cell) );

		cudaMalloc( (void**)&dev[n].orb_rot, dev[n].xarr*dev[n].zarr*sizeof(double) );
		cudaMalloc( (void**)&dev[n].orb_res, dev[n].xarr*dev[n].zarr*sizeof(double) );
		cudaMalloc( (void**)&dev[n].orb_shf, dev[n].xarr*dev[n].zarr*sizeof(int) );
		cudaMalloc( (void**)&dev[n].orb_inc, dev[n].xarr*dev[n].zarr*sizeof(int) );

		cudaMalloc( (void**)&dev[n].shfL, xpad*dev[n].zres*sizeof(int) );		
		cudaMalloc( (void**)&dev[n].shfR, xpad*dev[n].zres*sizeof(int) );		

		cudaMalloc( (void**)&dev[n].planets, n_planet*sizeof(body) );
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
}

void save_check_point(ofstream &check_point, string fname, int &sstep, double cur_time, Grid* hst)
{
	open_binary_file(check_point,fname);
	write_check_point(check_point, cur_time, hst);
	close_output_file(check_point);
	return;
}

int main(int narg, char *args[])
{
	string label=create_label();
	//string path =path_to_cwd()+"/binary/";
	string path = "/mnt/penguin/fung/p2/";
	string fname;

	Grid* hst = new Grid[ndev];
	Grid* dev = new Grid[ndev];

	init_grid_dimensions(hst);
	init_grid_dimensions(dev);
	mem_allocation(hst, dev);

	init(hst);
	
	for (int i=0; i<ndev; i++)
	for (int n=0; n<n_planet; n++)
	{
		hst[i].planets[n].m = planet_mass;
		hst[i].planets[n].x = planet_radius;
		hst[i].planets[n].y = pi;
		hst[i].planets[n].z = hpi;

		hst[i].planets[n].vx = 0.0;
		hst[i].planets[n].vy = pow(hst[i].planets[n].x,-1.5);
		hst[i].planets[n].vz = 0.0;
		#if ndim==2
		hst[i].planets[n].rs = 0.5*sc_h;
		#else
		hst[i].planets[n].rs = 0.5*pow(planet_mass,1.0/3.0);
		#endif
	}

	cpy_grid_HosttoDevice(hst, dev);

	///////////////////////////////////////////////////////////

	init_OrbAdv(dev);

	///////////////////////////////////////////////////////////

	double cur_time = sta_time;
	double sav_time = 0.0;
	double dt = end_time-cur_time;
	int nstep = 0;
	int sstep = 0;
	int tmp;

	ofstream check_point;
	fname = path+"binary_"+label+"_"+frame_num(sstep);
	#ifdef dump_flag
	cpy_grid_DevicetoHost(hst, dev);
	save_check_point(check_point, fname, sstep, cur_time, hst);
	sstep++;
	#endif

	auto clock = chrono::high_resolution_clock::now();
	auto begin = chrono::high_resolution_clock::now();
	double duration;
	double speed;
	double old_dt;

	bool print = false;
	bool savep = false;
	bool contu = true;

	#ifdef visc_flag
	viscosity_tensor_evaluation(dev);
	#endif

	while (contu)
	{
		dt = global_dt(dev, dt);
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
		if (nstep%prt_interval==0) print = true;
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
			cpy_grid_DevicetoHost(hst, dev);
			fname = path+"binary_"+label+"_"+frame_num(sstep);
			save_check_point(check_point, fname, sstep, cur_time, hst);

			printf("Check point!\n");
			printf("  Current time %f | End time %f | Time step %e\n", cur_time, end_time, dt);
			printf("  %s is saved. \n\n", fname.c_str());

			sstep++;
			sav_time = 0.0;
			savep = false;
			dt = old_dt;
		}
	}

	return 0;
}
