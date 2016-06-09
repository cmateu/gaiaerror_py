# gaiaerror_py

This is a simple Python wrapper for Merce Romero-Gomez's Gaia error code.

**DESCRIPTION**

This version of Merce Romero-Gomez's code includes the evaluation of the Drimmel extinction map.

Codes provided:

	- gaiaerr.py: Simulates Gaia errors 
	- gaiaerr_photerr.py: Simulates Gaia errors for all quantities except distance. For the distance, (gaussian) errors are simulated with a fixed relative error. 

**INSTALL**

- First run the compile script inside the gaia_errors_color_tmission folder
- Edit the gaiaerr.py and gaiaerr_photerr.py files and substitute the appropriate path in the gerr_path variable at the very beginning of the file

**REQUIREMENTS**

The numpy and scipy libraries are required
