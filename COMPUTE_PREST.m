                    function  COMPUTE_PREST
global U A1 ISTEP DFI H NN NIE EQ1 SEIS NSTEP
%
NIN = 2*NIE;
%SEISMIC LOADING
EQ = [EQ1,zeros(NSTEP,NIN)];
% PRELOCATION FOR RIGHT HAND SIDE VECTOR "A1*X=FI"
FI = zeros(NN,1);
% SOLVE EQUATION IN FIRST TIME STEP
if ISTEP == 1
    if SEIS == 1
        FI = FI + EQ(ISTEP,:)';
    end
    DFI = A1\FI;
end
% CONTRIBUTION OF PREVIOUS STEPS ON RIGHT HAND SIDE VECTOR IN TIME STEP>1
if ISTEP > 1
        for imat = 1 : ISTEP-1
            Un = U(:,imat);
            iut = ISTEP-imat+1;
            Hn = H(:,:,iut);
            FI = FI - Hn*Un ;
        end
        if SEIS == 1
            FI = FI + EQ(ISTEP,:)';
        end
    DFI = A1\FI;
end
U(:,ISTEP) = DFI;
clear DFI FI
