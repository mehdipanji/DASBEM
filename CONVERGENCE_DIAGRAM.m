global NN NN
global nstep delt U UI 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coding = U for checking the convergence of Boundary Elements (Surface of Inclusion)
% coding = UI for checking the convergence of Internal Points (Gorund Surface)
coding = U;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = linspace(0,nstep*(delt),nstep);
pp =1:length(coding);
A=size(coding);
A(1,1);
    plot(t,(coding(A(1,1),:)/0.001)/1,'k-')
    text(1,1,num2str(size(coding)))