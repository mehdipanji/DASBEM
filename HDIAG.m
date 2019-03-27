function  [ci] = HDIAG(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5)
%
zeta1 = atan2(-0.5*x1+2*x2-1.5*x3,0.5*y1-2*y2+1.5*y3);
zeta2 = atan2(1.5*x3-2*x4+0.5*x5,-1.5*y3+2*y4-0.5*y5);
zet = (zeta2-zeta1)/pi;
if abs(abs(zet)-1) < 1e-4
    disp(' ERROR IN HDIAG FUNCTION, CONCIDENT ELEMENTS ');
end
const = real(zet);
if const < 0
    ci = 0.5*(-1-zet);
elseif const == 0
    ci = 0.5*(1-zet);
elseif const > 0
    ci = 0.5*(1-zet);
end
if ci < 0
    ci = -ci;
end