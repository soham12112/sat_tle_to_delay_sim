# Two Line Element Reader / Orbit Plotter #

Given multiple NORAD Two-Line-Element (TLE) files, this matlab code plots the orbits of the satellites around Earth. 
[1]	T. G. R. Reid, "Orbital Diversity for Global Navigation Satellite Systems," Doctor of Philosophy, Aeronautics and Astronautics, Stanford University, Stanford, CA, 2017.

This thesis is available at the following link: https://purl.stanford.edu/dc409wn9227

## How to Use ##

'MAIN_TLE_READ_Multifile.m' reads in specified TLE files and plots their orbits around Earth. TLE orbit data be downloaded from the Celestrak website: http://www.celestrak.com/NORAD/elements/master.asp

'MAIN_TLE_READ_Multifile_ALL_TLE.m' does the same as 'MAIN_TLE_READ_Multifile.m' only is preset to read the TLE data for all operational orbits and plot them. 

## Position and Navigation use ## 

Data File: satellite_simulation_data.mat

This file contains the simulation output in a MATLAB struct called exportData.

Data Structure:

exportData.time_utc_datenum: A 1D array of time steps in MATLAB's datenum format (which is UTC).

exportData.satellite_names: A list of all satellite names (e.g., "IRIDIUM 140"). The index of a satellite in this list corresponds to its data in the arrays below.

State Vectors (3D Arrays): [Time, [X,Y,Z], Satellite_Index]

exportData.position_eci_m: Satellite position [X, Y, Z] in meters (ECI frame).

exportData.velocity_eci_m_s: Satellite velocity [X, Y, Z] in m/s (ECI frame).

Observation Data (2D Arrays): [Time, Satellite_Index]

exportData.elevation_ue1_deg: Elevation angle (degrees) from UE 1.

exportData.doppler_ue1_hz: Doppler shift (Hz) at UE 1.

exportData.elevation_ue2_deg: Elevation angle (degrees) from UE 2.

exportData.doppler_ue2_hz: Doppler shift (Hz) at UE 2.

Ground Station Info:

exportData.ue1_location_deg_m: [Lat, Lon, Alt] for UE 1.

exportData.ue2_location_deg_m: [Lat, Lon, Alt] for UE 2.