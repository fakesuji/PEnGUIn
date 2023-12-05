# PEnGUIn

** input parameters are located in include/parameters.h **

command lines for compiling instructions:

	1. If compiling on a machine with GPUs installed, go to the directory where the PEnGUIn code is, and enter:
		make
	Note: Re-make if the parameter file is changed.

	2. If compiling on a machine without GPUs installed, the Makefile must be modified. 
 	Go to Line 5 of Makefile, replace:
		-arch=native
	with the architecture of your choice. For example, choose
		-arch=sm_70
	for V100 GPUs.
