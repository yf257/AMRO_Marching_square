# AMRO_Marching_square
This program is used to simulate the AMRO data by numerically solving the boltzmann transport equation. 

## AMRO(Angle-dependent magnetoresistance):
Angle-dependent magnetoresistance is a measurement of the magnetoresistance of solidstate sample as a function of the angle between the crystal axis and the magnetic field, which provides information about the Fermi surface shape of the material.  
## FindFermi:
The first step of the AMRO calculation is to find a uniform grid of the Fermi surface given the bandstructure of the material. This version of the program use the self-developed marching square algorithm. The marching square algorithm is a computer graphics algorithm that generates contours for a two-dimensional scalar field. In our case, the two-dimensional scalar field is the energy band structure with a fixed momentum in z direction(kz). By iterating through a evenly spaced set of kz, we obtain the 3D fermi surface grid.  
