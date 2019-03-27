# DASBEM
This is a program for dynamic analysis of plane scalar time domain problems using half-plane boundary element method (BEM). This program is developed for analysis of 2D subsurface inclusions embedded in the elastic half-plane.

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
3. See DAGRAM to check "does the convergence is achieved in responses or not?" 
4. Use FFT to show 2D normalized displacement amplitude (NDA) of the surface.
5. Use the DIAGRAM_3D to show 3D Time and frequency domain responses.
