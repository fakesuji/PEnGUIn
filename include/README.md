
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

//=======================================================================
// Temporal parameters
//=======================================================================

	frame_omega:
 		-frame rotation frequency in code units
	sav_interval:
 		-time interval between saving snapshots (see dump_flag) in code units
   	end_time:
    		-end time of the simulation in code units
	prt_interval:
 		-number of timesteps between producing terminal output (see silence_flag)
	max_step:
 		-maximum number of timestep

//=======================================================================
// Hydro parameters
//=======================================================================

	EOS_flag:
 		-options: 0,1,2
   		-0: isothermal
     		-1: isentropic
       		-2: adiabatic
	internal_e_flag
 		-options: 0, 1
   		-0: track full energy when EOS is adiabatic
     		-1: track internal energy when EOS is adiabatic
       	gam:
		-the adiabatic index "gamma"
  	gamm:
   		-gamma minus one
     	gamp:
      		-gamma plus one
	gammfac:
 		-(gamma minus one)/(2 * gamma)
   	gampfac:
    		-(gamma plus one)/(2 * gamma)
	CFL:
 		-the Courant number
   
//=======================================================================
// Dust parameters
//=======================================================================

	D_G_ratio:
 		-dust-to-gas mass ratio
	Stokes:
 		-Stokes number

//=======================================================================
// boundary parameters
//=======================================================================

	init_flag:
 		-initial conditions
 		-options: 0,1,2,3,4,5,6,7,8
   		-0: one-dimensional shock tube test
     		-1: one-dimensional strong shock test
       		-2: two-dimensional protoplanetary disk in cylindrical coordinates
	 	-3: three-dimensional protoplanetary disk in spherical coordinates
   		-4: two-dimensional shock tube test
     		-5: two-dimensional Kelvin-Helmholtz test
       		-6: one-dimensional linear wave test
	 	-7: two-dimensional modified Kelvin-Helmholtz test
   		-8: one-dimensional divergence test
	bound_lft:
 		-lower x boundary conditions
   		-options: 0,1,2,3
     		-0: initial conditions
       		-1: zero gradient
	 	-2: reflective
   		-3: periodic
	bound_rgh:
  		-upper x boundary conditions
   		-options: 0,1,2,3
     		-0: initial conditions
       		-1: zero gradient
	 	-2: reflective
   		-3: periodic
	bound_bak:
  		-lower y boundary conditions
   		-options: 0,1,2,3
     		-0: initial conditions
       		-1: zero gradient
	 	-2: reflective
   		-3: periodic
	bound_frn:
  		-upper y boundary conditions
   		-options: 0,1,2,3
     		-0: initial conditions
       		-1: zero gradient
	 	-2: reflective
   		-3: periodic
	bound_bom:
  		-lower z boundary conditions
   		-options: 0,1,2,3
     		-0: initial conditions
       		-1: zero gradient
	 	-2: reflective
   		-3: periodic
	bound_top:
  		-upper z boundary conditions
   		-options: 0,1,2,3
     		-0: initial conditions
       		-1: zero gradient
	 	-2: reflective
   		-3: periodic
     
//=======================================================================
// Disk parameters
//=======================================================================

	p_beta:
 		-disk temperature follows r^-p_beta
   	p_alpha:
    		-In 2D simulations, dist surface density follows r^-p_alpha
      		-In 3D simulations, dist midplane density follows r^-p_alpha
	ss_alpha:
 		-the Sunyaev-Shakura alpha viscosity paramter
   	sc_h:
    		-disk aspect ratio at r=planet_radius
	Sigma_0:
 		-In 2D simulations, disk surface density at r=1 in units of solar mass per au^2
   		-In 3D simulations, disk midplanet density at r=1 in units of solar mass per au^3
     	kill_width:
      		-width of the killing wave in units of the local disk scale height
	beta_cool:
 		-the thermal relaxation time in units of the local dynamical time

//=======================================================================
// planet parameters
//=======================================================================

	n_planet:
 		-number of planets (planet + star if twobd_flag is set to 1)
	planet_mass:
 		-planet-to-star mass ratio
   	planet_radius:
    		-semi-major axis of the planet's orbit in code units
	planet_ecc:
 		-eccentricity of the planet's orbit
	ramp_time:
 		-initial time, in code units, taken to gradually bring the planet's mass to "planet_mass"
   
//=======================================================================
// Grid parameters
//=======================================================================

	xmin:
 		-location of the lower x boundary
	xmax:
 		-location of the upper x boundary
	ymin:
 		-location of the lower y boundary
	ymax:
 		-location of the upper y boundary
	zmin:
 		-location of the lower z boundary
	zmax:
 		-location of the upper z boundary
	geomx, geomy, geomz:
 		-grid geometry: cartesian, cylindrical, spherical
   		-cartesian: geomx=0, geomy=0, geomz=0
     		-cylindrical: geomx=1, geomy=3, geomz=0
     		-spherical: geomx=2, geomy=4, geomz=5
	gridx, gridy, gridz:
		-grid spacing options:0,1,2,3,4
  		-0: uniform
    		-1: logarithmic
      		-2: non-uniform spacing that concentrates cells near the center of the grid
		-3: non-uniform spacing that concentrates cells near the top of the grid
  		-4: mixed non-uniform and uniform spacing that concentrates cells near the center of the grid
