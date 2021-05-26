# SWeVDF_stability_inputs
Input files for LEOPARD code to produce data shown in Schroeder et. al. 2021 (in prep).

The code used, LEOPARD, is publically available from https://github.com/pastfalk/LEOPARD.

Files here include a python script to compute the analytical electron distribution function composed of a thermal core component and a superthermal, field-aligned beam component called the "strahl". The python file is written to be integrated with included shell script files for setting up and running scans over radius and angle of propagation.

To set up a LEOPARD scan, first compile LEOPARD according to instructions in its repository. One must also have a working installion of python3, numpy, and matplotlib. In a bash shell, cd to a directory (e.g. this git repository) containing the python file, shell file, and input file with parameters for the k range and initial guess (specifies wave branch) specified. Then simply type "sh setup_input_xradius_yangle.sha b c d e f g ..." where (x,y) is either (n,1) or (1,n) depending on if one wants to scan n radii at fixed angle, or n angles at fixed radii, and (a,b,c,...) is a series of floats where "a" is always the fixed quantity and (b,c,...) is a series of floats of any length corresponding to the n radii/angles requested. Here angles are given in degrees between 0 and 90, and radii are given in 0.01 AU units (e.g. 1 AU is specified as "100" -- this choice was made such that directories are created with integer names).

For example, if one wants LEOPARD to scan the FM-whistler branch at a fixed distance of 1.5 AU and for integer angles between 20 and 40 degress, ensure the input.dat file has the right initial omega guesses for the FM-whist branch, then run "sh setup_input_1radius_nangles.sh 150 {20..40}". This will set up a directory tree such that the outermost directories are the distance, and subdirectories are angles, and in each angle directory there will be an input file and distribution directory containing a distribution file for electrons. 

After setting up files, to actually run LEOPARD the "runLEOPARD_xradius_yangle.sh" behaves similarly to the setup file. For the example given above, simply type "sh runLEOPARD_1radius_nangle.sh 150 {20..40}" and the script will iterate through each angle at the given distance and LEOPARD sequentially for each case.   

Other examples: 

run at 50 degrees for radii 0.1 AU, 0.5AU, 2AU, 1AU --- 
	1. make sure input.dat file has good initial guess for the desired wave branch and angles
	2. run "sh setup_input_nradius_1angle.sh 50 10 50 200 100"
	3. run "sh runLEOPARD_nradius_1angle.sh 50 10 50 200 100"

Note: current input.dat file here is specified for FM-whistler branch. Next commit will include other input.dat files for KAW branch at various angles of propagation.
