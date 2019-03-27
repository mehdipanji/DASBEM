                        function  COMPUTE_AFMAT
global ISTEP H A1
% FORM A1 AND B1 MATS IN FIRST TIME STEP. A1*X=B1*(B.C.)
if ISTEP == 1
    A1 = H(:,:,1);
end
