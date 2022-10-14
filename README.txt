Matlab code and synthetic data used in Albright and Gregg (2022) to test the effectiveness of different variants of the Ensemble Kalman Filter (EnKF) when applied to
geodetic displacement at simulated unresting volcanic systems.

***AssimilateClean_11Oct22.m*** - Matlab function for completing an EnKF analysis of the synthetic data. Used in most experiments conducted.

***AssimilateNoReparam_11Oct22.m***: - Matlab function for completing an EnKF analysis of the synthetic data. Used for the variant without parameter reformulation.


These functions require the following dependencies in the same folder:

***Synth_P_06Mar21.mat and Synth_R_06Mar21.mat***: These files each contain the synthetic GPS and InSAR data used in this study as well as the parameter time series
used to generate them. The "P" variant corresponds to a simulation in which magma inflation is driven by a pressurizing reservoir, while the "R" variant uses a 
laterally expanding reservoir under constant pressure. Their contents are listed below:

	- "GPS": Matlab cell structure containing synthetic GPS data. Each row of the cell array corresponds to a different GPS station. 
	   In order from left to right, the columns contain:
		1. Station xy coordinates relative to model space origin
		2. List of times with available data at that station, in seconds relative to the start of the simulation
		3-5. East-West, North-South, and Up-Down displacement of the station at each time listed in 2 [m].
		6-8. Uncertainties in (3)-(5), respectively [m].
		9. Station name

	- "INSAR": Matlab cell structure containing synthetic InSAR data. Each row corresponds to a different downsampled, cumulative,
	   unwrapped InSAR image. In order from left to right, the columns contain (d = # of observation points in a given image):
		1. 3x1 matrix containing: Decimal date of image collection, satellite look angle [deg], satellite track azimuth [deg]
		2. d-by-2 matrix containing the xy coordinates of each observation point [m,m]
		3. Line-of-sight displacement at each point listed in (2) [m]
		4. Uncertainty in (3) [m]
		5. d-by-10 matrix containing the x (columns 1-5) and y (columns 6-10) coordinates of the corners of the quadrilateral area
		   used to downsample the given data point.
			- Corners are listed in the following order: NW, NE, SE, SW, NW
		6. The number of original data points condensed into a given observation as part of the downsampling process

	- "Quads": Quadtrees quadrilateral corner coordinates, same as column (5) in "INSAR". The same quads derived from the last image were 
	   used to downsample all original synthetic images.

	- "TFs": Tensile failure properties of the original synthetic model. Each entry of the cell array is a 4-by-q matrix, each column 
	   corresponding to one of q finite element cells along the wall of the simulated magma reservoir. The rows of this matrix 
	   contain:
		1. The tensile stress experienced by the given model element [Pa, tension positive]
		2-4. [Unused, related to position of each element]

	- "steps_s": Vector containing the time steps at which the synthetic model space was evaluation [s]

	- "r1s, r2s, dPs, dXs, dYs, dZs" - Synthetic values of each model parameter at each time step [m,m,Pa,m,m,m]

	- "tol" - variance added as gaussian noise to InSAR observations [m]



**Elastic2D_17Jun19.m**: Given reservoir parameters, runs a COMSOL FEM and extracts surface displacements

**2DElasticEllip_17Jun19.mph**: COMSOL model file containing pre-built magma reservoir FEM
	- **ElasticEllip2D.mat** generates the same COMSOL FEM

**EnKF_14Jun21.m**: Updates an ensemble of parameter states based on available data using one of several variants on the EnKF algorithm

**EnKFParamShift.m**: Converts base reservoir parameters to reformulated unitless parameters and vice versa.
	- Dependencies: **VlimsTheta.m** **logist_18Jun21.m**

**imp_redist_09Apr21.m**: Corrects physically impossible values produced by the EnKF update. Only used in AssimilateNoReparam_11Oct22.m

**qfilt_rows.m**: Interquartile filter to remove outliers

**YangEllipsoid**: Folder of functions to conduct an analytical solution for surface deformation, after Yang et al (1988). Not used in published experiments and cannot 
account for reservoir stress, but allows for faster calculation while testing workflow.


***ParamExtract_23Feb22.m***: This function processes the large output files of the "Assimilate_*.m" functions above and calculates the filter's performance at the
final time step relative to the original synthetic model through 4 quantitative criteria: RMSE, reservoir wall tensile stress misfit, stable parameter misfit, and 
non-unique parameter misfit. When multiple files are loaded into this function, their results are additionally analyzed as a group for the purpose of calculating 
the mean and percentile values.
