# DASBEM
This is a program for dynamic analysis of plane scalar time domain problems using half-plane boundary element method (BEM). This program is developed for analysis of 2D subsurface inclusions embedded in the elastic half-plane.

•	INPUT DATA: Read the input file created by WRITE DATA sub-problem. These data are including shear wave velocity, number of time steps, number of boundary nodes and their coordinates, boundary conditions of each node, number of internal points and loading data if it exists.
•	COMPUTE GHMAT: Compute the elements of G and H matrices using the numerical integration.
•	NONSING: Compute the non-singular elements of mentioned matrices.
•	SING: Determine the singular elements of G matrix using special numerical integration.
•	FUNSOLE: Including the half-plane elastodynamic fundamental solutions. This part is called by NONSING and SING in each numerical integration.
•	HDIAG: Compute the diagonal elements of H matrix. 
•	COMPUTE GHARG: Assembling the matrices which are calculated in previous steps.
•	COMPUTE AFMAT: Create the solvable form of equations based on the boundary conditions of each boundary node.
•	FUNCVAL: This sub-program is calls by COMPUTE AFMAT when the loading is in the form of external vibrations. 
•	COMPUTE UTMAT: In this part of program, the borders of u and q are formed based on the boundary conditions.
•	COMPUTE PREST: Considering the past dynamic time history and calling the border of free-field motion in current time step.
•	COMPUTE INTP: Compute the responses of internal points.
•	OUTPUT DATA: Print the responses of boundary nodes and internal points.
---------------------------------------------------------------------------------------------------------------------------
How to use:
1. Open the WRITE_DATA and set the following values, respectively.
a) The shear modulus of domain and inclusion. (The assumed values are '16e5' and '3.24e5', respectively).
b) The density of domain and inclusion. (The assumed values are '1' and '0.6667', respectively).
c) The number of time steps and time increments. (The assumed values are '400' and 0.0175', respectively).
d) The radius (b) and depth of inclusion (D). (The assumed values are '500' and 750', respectively).
e) The shape of inclusion (circular or elliptical). (The assumed shape is circular).
f) The size of half-element. (The assumed value is '30').
g) The predominant frequency. (The assumed value is '3').
h) The time shift parameter. (The assumed value is '2.4').
i) The angle of incident wave. (The assumed angle is '0').
j) The range of ground surface "internal points". (The assumed value is '-3b ~ 3b').
2. Open the MAIN_PROGRAM and run it and wait until finishing the calculations.
3. See CONVERGENCE_DIAGRAM to check "does the convergence is achieved in responses or not?" 
4. Use FFT to show 2D normalized displacement amplitude (NDA) of the surface.
5. Use the DIAGRAM_3D to show 3D Time and frequency domain responses.
