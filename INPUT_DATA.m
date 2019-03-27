                     function  INPUT_DATA
global X Y NE NFUNC NN CS1 CS2 NSTEP 
global AT GE1 RO1 GE2 RO2 NIE
global radif fid out name start L XI YI SEIS 
%
% READ NAME OF INPUT FILE
file='input.txt'; 
out = 'output.txt';
start = clock;
file = lower(file);
% OPEN INPUT FILE ID
fid = fopen(file);
name = fscanf(fid,'%s ',1); % Name of Example
%
GE1 = fscanf(fid,'%f',1); % SHEAR MODULUS OF DOMAIN
RO1 = fscanf(fid,'%f',1); % DENSITY OF DOMAIN
GE2 = fscanf(fid,'%f',1); % SHEAR MODULUS OF LINING
RO2 = fscanf(fid,'%f',1); % DENSITY OF LINING
%
CS1 = sqrt(GE1/RO1); % WAVE VELOCITY OF DOMAIN
CS2 = sqrt(GE2/RO2); % WAVE VELOCITY OF LINING
%
NSTEP = fscanf(fid,'%i',1); % NUMBER OF TIME INTERVALS
AT = fscanf(fid,'%f',1);    % SIZE OF TIME INTERVALS
NE = fscanf(fid,'%i',1);    % NUMBER OF ELEMENTS
NIE = fscanf(fid,'%i',1);   % NUMBER OF INTERFACE ELEMENTS (OUTER BOUNDARY)
NN = 2*NE; % NN IS NUMBER OF NODES
for i = 1 : NN
    radif = fscanf(fid,'%i',1);
    X(i) = fscanf(fid,'%f',1); % X DIR. OF NODE
    Y(i) = fscanf(fid,'%f',1); % Y DIR. OF NODE
end
% READ EARTHQUAKE CODE
SEIS = fscanf(fid,'%i',1);
%  READ DATA FOR LOADING FUNCTIONS
NFUNC = fscanf(fid,'%i',1); % NUMBER OF TIME FUNCTIONS
%
% READ DATA FOR INTERNAL POINTS
L = fscanf(fid,'%i',1); % NUMBER OF INTERNAL POINTS
%
for i = 1 : L
    radif = fscanf(fid,'%i',1);
    XI(i) = fscanf(fid,'%f',1);
    YI(i) = fscanf(fid,'%f',1);
end
fclose(fid);