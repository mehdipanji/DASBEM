function [HT,GT] = NONSING(xp,yp,x1,y1,x2,y2,x3,y3,ISTEP,AT,CS)
% NUMERICAL INTEGRATION FOR CALCULATION OF 2 X 2 SUBMAT. G & H
kes(1) = 0.9739065285;       kes(2) = -0.9739065285;
kes(3) = 0.8650633666;       kes(4) = -0.8650633666;
kes(5) = 0.679409568;        kes(6) = -0.679409568;
kes(7) = 0.4333953941;       kes(8) = -0.4333953941; 
kes(9) = 0.1488743389;       kes(10) = -0.1488743389;
%
w(1) = 0.0666713443;         w(2) = 0.0666713443;
w(3) = 0.1494513491;         w(4) = 0.1494513491;
w(5) = 0.2190863625;         w(6) = 0.2190863625;
w(7) = 0.2692667193;         w(8) = 0.2692667193;
w(9) = 0.2955242247;         w(10) = 0.2955242247;
%
HT = zeros(1,3);
GT = zeros(1,3);
%
a = x3-2*x2+x1;
b = (x3-x1)/2;
c = y3-2*y2+y1;
d = (y3-y1)/2;

%
cst = CS*AT*ISTEP;
for i = 1 : 10
    f1 = kes(i)*(kes(i)-1)*0.5;
    f2 = 1-kes(i)^2;
    f3 = kes(i)*(kes(i)+1)*0.5;
    % calculation of real part of fundamental solution
    xco = (x1)*f1+(x2)*f2+(x3)*f3-xp;
    yco = (y1)*f1+(y2)*f2+(y3)*f3-yp;
    ra = sqrt(xco^2+yco^2);
    if ra > cst
        continue
    end
    %
    [hr,gr] = FUNDSOLE(ra,ISTEP,AT,CS);
    %
    rd1 = xco/ra;
    rd2 = yco/ra;
    jaw = sqrt((kes(i)*a+b)^2+(kes(i)*c+d)^2)*w(i);
    %
    eta1 = c*kes(i)+d;
    eta2 = -(a*kes(i)+b);
    rdn = rd1*eta1+rd2*eta2;
    %
    GT(1) = GT(1) + gr*jaw*f1;
    GT(2) = GT(2) + gr*jaw*f2;
    GT(3) = GT(3) + gr*jaw*f3;
    %
    HT(1) = HT(1) + hr*rdn*f1*w(i);
    HT(2) = HT(2) + hr*rdn*f2*w(i);
    HT(3) = HT(3) + hr*rdn*f3*w(i);
    % calculation of image part of fundamental solution
    ypco = (y1)*f1+(y2)*f2+(y3)*f3+yp;
    rpa = sqrt(xco^2+ypco^2);
    if rpa > cst
        continue
    end
    %
    [hi,gi] = FUNDSOLE(rpa,ISTEP,AT,CS);
    %
    rd1 = xco/rpa;
    rd2 = ypco/rpa;
    rpdn = rd1*eta1+rd2*eta2;
    %
    GT(1) = GT(1) + gi*jaw*f1;
    GT(2) = GT(2) + gi*jaw*f2;
    GT(3) = GT(3) + gi*jaw*f3;
    %
    HT(1) = HT(1) + hi*rpdn*f1*w(i);
    HT(2) = HT(2) + hi*rpdn*f2*w(i);
    HT(3) = HT(3) + hi*rpdn*f3*w(i);
end
%