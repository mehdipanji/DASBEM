function  COMPUTE_GHARG
global ISTEP NSTEP NN NIE Gd1 Hd1 Gd2 Hd2 H GE1 GE2 
NIN = 2*NIE;
if ISTEP == 1
    H = zeros(NN,NIN+3*NIE,NSTEP);
end
H(1:NIN,1:NIN,ISTEP) = H(1:NIN,1:NIN,ISTEP) + Hd1(:,:);
H(1:NIN,NIN+1:end,ISTEP) = H(1:NIN,NIN+1:end,ISTEP) - (1/GE1)*Gd1(:,:);
H(NIN+1:end,1:NIN,ISTEP) = H(NIN+1:end,1:NIN,ISTEP) + Hd2(:,:);
H(NIN+1:end,NIN+1:end,ISTEP) = H(NIN+1:end,NIN+1:end,ISTEP) + (1/GE2)*Gd2(:,:);
clear Gd1 Hd1 Gd2 Hd2
