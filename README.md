# gaiaerror_py

This is a simple Python wrapper for Merce Romero-Gomez's Gaia error code.

**DESCRIPTION**

This version of Merce Romero-Gomez's code includes the evaluation of the Drimmel extinction map.

Codes and files provided:

	* gaiaerr.py: Simulates Gaia errors 
	* gaiaerr_photerr.py: Simulates Gaia errors for all quantities except distance. For the distance, (gaussian) errors are simulated with a fixed relative error. 
	* example_data/myfile.ne.dat: Example input file


**INSTALL**

- First run the compile script inside the gaia_errors_color_tmission folder
- Edit the gaiaerr.py and gaiaerr_photerr.py files and substitute the appropriate path in the gerr_path variable at the very beginning of the file

**REQUIREMENTS**

The numpy and scipy libraries are required

**QUICK GUIDE**

The codes take input files (ascii) containing  X Y Z VX VY VZ Mv V-I, spatial coordinates must be in kpc and velocities in km/s. Mv and V-I are the absolute magnitude and V-I intrinsic colour.

You can use the example input file provided to run the Gaia error simulation as:

	gaiaerr.py myfile.ne.dat

The output file will be myfile.pe.dat.

To simulate 25% distance errors plus Gaia errors for the remaining quantities, run the gaiaerr_photerr.py code as:

	gaiaerr_photerr.py myfile.ne.dat 0.25



