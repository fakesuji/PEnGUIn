//////////////////////////////////////////////////
//                                              //
//              PEnGUIn Data File               //
//                                              //
//////////////////////////////////////////////////

File name:
	PEnGUIn simulation snapshot files are all named in the following format:
		binary_"resolution"_"label"_"frame number"
	where "resolution" is either a single number indicating the number of cells in 1D simulation, two number indicating the number of cells in the x and y directions, or three numbers indicating the number of cells in x, y, and z directions;
	"label" is a set of string assigned to the simulation that is primarily used as a identifier; and
	"frame number" is an index assigned to each snapshot, starting from zero.

Reading files:
	Output files are in binary format and can be read using "read_penguin.py". 
	Use functions "load_1D_data", "load_2D_data", and "load_3D_data" to load 1/2/3D simulation snapshots. 
	For example, loading the 2D simulation file "binary_600x1080_h50_2p1J_e50_a-30_b-20_OA_PPM4_00100" located in the directory "/home/xxx/" one should use:
		data = load_2D_data("/home/xxx/", 600, 1080, "h50_2p1J_e50_a-30_b-20_OA_PPM4", 100)
	where the first input is the path to the file;
	the second input is the x resolution;
	the third input is the y resolution;
	the fourth input is the label;
	and the fifth input is the frame number.
	"data" is the ouput python array that contains the snapshot.

Data structure:
	After reading the file, the output python array has the following format.
	For 1D simulations:
		data[0] is a double that contains the simulation time of the snapshot in code units
		data[1] is an array of doubles that contains the locations of all cell boundaries
		data[2] is an array of doubles that contains the density in each cell
		data[3] is an array of doubles that contains the pressure in each cell
		data[4] is an array of doubles that contains the velocity in each cell
	For 2D simulations:
		data[0] is a double that contains the simulation time of the snapshot in code units
		data[1] is an array of doubles that contains the locations of all cell boundaries in the x direction
		data[2] is an array of doubles that contains the locations of all cell boundaries in the y direction
		data[3] is a 2D array of doubles that contains the density in each cell with dimensions of (y-resolution, x-resolution)
		data[4] is a 2D array of doubles that contains the pressure in each cell with dimensions of (y-resolution, x-resolution)
		data[5] is a 2D array of doubles that contains the x-velocity in each cell with dimensions of (y-resolution, x-resolution)
		data[6] is a 2D array of doubles that contains the y-velocity in each cell with dimensions of (y-resolution, x-resolution)
	For 3D simulations:
		data[0] is a double that contains the simulation time of the snapshot in code units
		data[1] is an array of doubles that contains the locations of all cell boundaries in the x direction
		data[2] is an array of doubles that contains the locations of all cell boundaries in the y direction
		data[3] is an array of doubles that contains the locations of all cell boundaries in the z direction
		data[4] is a 3D array of doubles that contains the density in each cell with dimensions of (z-resolution, y-resolution, x-resolution)
		data[5] is a 3D array of doubles that contains the pressure in each cell with dimensions of (z-resolution, y-resolution, x-resolution)
		data[6] is a 3D array of doubles that contains the x-velocity in each cell with dimensions of (z-resolution, y-resolution, x-resolution)
		data[7] is a 3D array of doubles that contains the y-velocity in each cell with dimensions of (z-resolution, y-resolution, x-resolution)
		data[8] is a 3D array of doubles that contains the z-velocity in each cell with dimensions of (z-resolution, y-resolution, x-resolution)

	
