
//////////////////////////////////////////////////////////
//                                                      //
//                                                      //
//     This README documents the variables in the       //
//     "parameters.h" file.
//                                                      //
//                                                      //
//////////////////////////////////////////////////////////

	dump_flag:
		-options: def or not def
		-def: save snapshots in "path_for_dumnp" at intervals of "sav_interval"
		-not def: snapshots not saved
	kill_flag:
		-options: 0 or 1
		-0: killing zones not included
		-1: killing zones included at both the inner and outer x-boundary with widths of "kill_width"
	OrbAdv_flag:
		-options: def or not def
		-def: orbital advection applied
		-not def: orbital advection not applied
	visc_flag:
		-options: 0 or 1
		-0: no explicit viscosity included
		-1: viscosity tensor is calculated and applied using paramter "ss_alpha", the Sunyaev-Shakura alpha prescription.
	cool_flag:
		-options: 0 or 1
		-0: simulation is strictly isothermal, adiabatic, or isentropic (see EOS_flag)
		-1: a relaxation term is included with a relaxation time determined by "beta_cool"
	twobd_flag:
		-options: 0 or 1
		-0: simulation centers on the star
		-1: simulation is barycentric
	silence_flag:
		-options: def or not def
		-def: turns off terminal output 
		-not def: terminal output at intervals defined by "prt_interval"
	advec_flag:
		-options: def or not def
		-def: solves advection equations only
		-not def: solves hydrodynamics equations

//=======================================================================
// Save path
//=======================================================================

	path_for_dump:
		-location of snapshot files
	
//=======================================================================
// Constants
//=======================================================================
	
	pi:
		-pi
	hpi:
		-0.5*pi
	twopi:
		-2.0*pi
	sqrt_hpi:
		-square root of pi
	third:
		-1/3
	fourthd:
		-4/3
	EarthMass:
		-one Earth mass in units of solar mass
	NeptuneMass:
		-one Naptune mass in units of solar mass
	SaturnMass:
		-one Saturn mass in units of solar mass
	JupiterMass:
		-one Jupiter mass in units of solar mass
	MMSN_1AU:
		-1700g per cm^2 in units of solar mass per au^2
	smallp:
		-pressure floor in code units
	smallr:
		-density floor in code units

//=======================================================================
// Geometric parameters
//=======================================================================

	recon_flag:
		-options: 0,1,2,3,4,5,6
		-0: PLM with van leer limiter
		-1: PLM with monotonized central limiter
		-2: PEM 2nd order boundary approximation
		-3: PEM 3rd order boundary approximation
		-4: PPM 3rd order boundary approximation
		-5: PEM 4th order boundary approximation
		-6: PPM 4th order boundary approximation
	std_thd:
		-standard number of threads per block
		-default 1024
	ndim:
		-options: 1,2,3
		-number of spatial dimensions
	npad:
		-number of padding zones
		-determined by "recon_flag"
	xpad, ypad, zpad:
		-number of padding zones in the x, y, and z direcitons
		-determined by "ndim" and "npad" 
	xres, yres, zres:
		-number of cells in the x, y, and z directions
		-yres=1 if ndim<2
		-zres=1 if ndim<3
	xarr, yarr, zarr:
		-array size in the x, y, and z directions
		-always equal x/y/zres + 2 * x/y/zpad

//=======================================================================
// Domain division parameters
//=======================================================================

	ndev:
		-number of GPUs being used
	x_xdiv, x_ydiv, x_zdiv:
		-the x, y, and z dimensions of the simulation domain in each block during a x-dimension sweep
	x_xthd:
		-the array size in each block during a x-dimension sweep
	y_xdiv, y_ydiv, y_zdiv:
		-the x, y, and z dimensions of the simulation domain in each block during a y-dimension sweep
	y_xthd:
		-the array size in each block during a y-dimension sweep
	z_xdiv, z_ydiv, z_zdiv:
		-the x, y, and z dimensions of the simulation domain in each block during a z-dimension sweep
	z_xthd:
		-the array size in each block during a z-dimension sweep

