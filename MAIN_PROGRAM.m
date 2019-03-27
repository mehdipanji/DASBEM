clear all;
clc
%                        PROGRAM   DYNA PSTDQBE

%                        Edited By :   Panji M. 

%  This is a Program for (DYN)amic (A)nalysis of (P)lane (S)calar (T)ime
%       (D)omain Problems using (Q)uadratic (B)oundary (E)lements.

%   This program is developed for analysis of sub-surface inclusions in 
%            elastic half-plane by half-plane time domain BEM

global NSTEP ISTEP NN U start
%                           " MAIN PROGRAM "
WRITE_DATA  % Write data into input file
INPUT_DATA  % Read data
for ISTEP = 1 : NSTEP
    tic;
    COMPUTE_GHMAT  % Form matrices G and H 
    COMPUTE_GHARG  % Arrange|Assemble G and H Matrices
    COMPUTE_AFMAT  % Form A1 and B1 matrices correspondig of B.C.
    COMPUTE_PREST  % Form right hand side vector for solution, adding the
                     % boundary condition vector and the previous step
                     % effect vector.
    file = ['Ti.St.',num2str(ISTEP),'.txt'];               
    nod = fopen(file,'wt');
    for ii = 1 : NN
        fprintf(nod,' %-10s\n',num2str(U(ii,ISTEP)));
    end
    fclose(nod);
    disp(['TIME STEP ( ',num2str(ISTEP),' ) : ',...,
        num2str(toc),'  sec']);
    disp(['TOTAL TIME UNTIL CURRENT STEP :',num2str(clock-start),'  sec']);
end
COMPUTE_INTP       % Compute Displacements at internal points
disp(['TOTAL TIME OF ANALYSIS : ',num2str(clock-start)]);