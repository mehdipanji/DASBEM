function  [HT,GT] = SING(xp,yp,x1,y1,x2,y2,x3,y3,AT,CS,nodo)
%
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
kesl(1) = 0.0090426309;      kesl(2) = 0.0539712662;
kesl(3) = 0.1353118246;      kesl(4) = 0.2470524162;
kesl(5) = 0.3802125396;      kesl(6) = 0.5237923179;
kesl(7) = 0.6657752055;      kesl(8) = 0.7941904160;
kesl(9) = 0.8981610912;      kesl(10) = 0.9688479887;
%
wl(1) = 0.1209551319;        wl(2) = 0.1863635425;
wl(3) = 0.1956608732;        wl(4) = 0.1735771421;
wl(5) = 0.1356956729;        wl(6) = 0.0936467585;
wl(7) = 0.0557877273;        wl(8) = 0.0271598109;
wl(9) = 0.0095151826;        wl(10) = 0.0016381576;
%
HT = zeros(1,3);
GT = zeros(1,3);
%
a = x3-2*x2+x1;
b = (x3-x1)/2;
c = y3-2*y2+y1;
d = (y3-y1)/2;
%
switch nodo
    case 1
        a1 = 0.5;
        a0 = 0.5;
        jal = 2.0;
    case 2
        a1 = 1.0;
        a0 = 0.0;
        jal = 1.0;
    case 3
        a1 = -0.5;
        a0 = 0.5;
        jal = 2.0;
end
%
%
for i = 1:10
    f1 = kes(i)*(kes(i)-1)*0.5;
    f2 = 1-kes(i)^2;
    f3 = kes(i)*(kes(i)+1)*0.5;
    %   Real Part of GT
    xco = (x1)*f1+(x2)*f2+(x3)*f3-xp;
    yco = (y1)*f1+(y2)*f2+(y3)*f3-yp;
    ra = sqrt(xco^2+yco^2);
    [hr,gr] = FUNDSOLE(ra,1,AT,CS);
    %
    ja = sqrt((kes(i)*a+b)^2+(kes(i)*c+d)^2);
    eta = abs(a1*kes(i)+a0);
    gr = (2*pi*gr+log(eta))*ja*w(i);
    kesa = (kesl(i)-a0)/a1;
    ja = jal*sqrt((kesa*a+b)^2+(kesa*c+d)^2)*wl(i);
    s1 = 0.5*kesa*(kesa-1)*ja + f1*gr;
    s2 = (1-kesa^2)*ja + f2*gr;
    s3 = 0.5*kesa*(kesa+1)*ja+ f3*gr;
    if nodo == 2
        kesa = -kesl(i);
        ja = jal*sqrt((kesa*a+b)^2+(kesa*c+d)^2)*wl(i);
        s1 = s1 + 0.5*kesa*(kesa-1)*ja;
        s2 = s2 + (1-kesa^2)*ja;
        s3 = s3 + 0.5*kesa*(kesa+1)*ja;
    end
    GT(1) = GT(1) + s1 / (2*pi);
    GT(2) = GT(2) + s2 / (2*pi);
    GT(3) = GT(3) + s3 / (2*pi);
    %   Image Part of GT
    xco = (x1)*f1+(x2)*f2+(x3)*f3-xp;
    ypco = (y1)*f1+(y2)*f2+(y3)*f3+yp;
    rpa = sqrt(xco^2+ypco^2);
    [hi,gi] = FUNDSOLE(rpa,1,AT,CS);
    %
    ja = sqrt((kes(i)*a+b)^2+(kes(i)*c+d)^2);
    eta = abs(a1*kes(i)+a0);
    gi = (2*pi*gi+log(eta))*ja*w(i);
    kesa = (kesl(i)-a0)/a1;
    ja = jal*sqrt((kesa*a+b)^2+(kesa*c+d)^2)*wl(i);
    s1 = 0.5*kesa*(kesa-1)*ja + f1*gi;
    s2 = (1-kesa^2)*ja + f2*gi;
    s3 = 0.5*kesa*(kesa+1)*ja+ f3*gi;
    if nodo == 2
        kesa = -kesl(i);
        ja = jal*sqrt((kesa*a+b)^2+(kesa*c+d)^2)*wl(i);
        s1 = s1 + 0.5*kesa*(kesa-1)*ja;
        s2 = s2 + (1-kesa^2)*ja;
        s3 = s3 + 0.5*kesa*(kesa+1)*ja;
    end
    GT(1) = GT(1) + s1 / (2*pi);
    GT(2) = GT(2) + s2 / (2*pi);
    GT(3) = GT(3) + s3 / (2*pi);
    % Calculate Real Part of HT
    cst = CS*AT;
    if ra >= cst
        continue
    end
    rd1 = xco/ra;
    rd2 = yco/ra;
    eta1 = c*kes(i)+d;
    eta2 = -(a*kes(i)+b);
    rdn = rd1*eta1+rd2*eta2;
    hr = hr*rdn*w(i);
    HT(1) = HT(1) + f1*hr;
    HT(2) = HT(2) + f2*hr;
    HT(3) = HT(3) + f3*hr;
    % Calculate Image Part of HT
    %
    cst = CS*AT;
    if rpa >= cst
        continue
    end
    rd1 = xco/rpa;
    rd2 = ypco/rpa;
    rpdn = rd1*eta1+rd2*eta2;
    hi = hi*rpdn*w(i);
    HT(1) = HT(1) + f1*hi;
    HT(2) = HT(2) + f2*hi;
    HT(3) = HT(3) + f3*hi;
end
