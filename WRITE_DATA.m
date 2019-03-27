function WRITE_DATA
global EQ1 fp t0 amax teta 
%------------------------------  input  -----------------------------------
g1 = 1*(16e5);  % shear modulus of domain
ro1 = 1;     % density of domain
g2 = 2/3*(3.24e5) ;  % shear modulus of inclusion
ro2 = 2/3;     % density of inclusion
nstep = 400; % number of time step
delt = (7/400); % time increment
% mesh generate
blim = 500;
% height of inclusion
h = 750;
%--------------calculation of cavity surface elements----------------------
tot = 1; 
%          
switch  tot
    case 1
        % circular inclusion
        % Part 1 (Outer boundary of inclusion)
        r = blim; % radius
        soe = 30;  % size of half-element
        non = r*pi/soe; % number of node
        non = round(non);
        tin = pi/non;
        tet = pi/2 : tin : 5*pi/2;
        tol = 0.000;
        while rem(length(tet),2) == 0
            soe = soe + tol; 
            non = r*pi/soe; 
            non = round(non);
            tin = pi/non;
            tet = pi/2 : tin : 5*pi/2;
            tol = tol + 0.001;
        end
        xt1(1:length(tet)) = 0; 
        yt1(1:length(tet)) = 0;
        for i = 1 : length(tet)
            xt1(i) = r * cos(tet(i));
            yt1(i) = -r * sin(tet(i))-h;
        end
        % Part 2 (Iner boundary of inclusion)
        r = blim; % radius
        soe = 30;  % size of half-element
        non = r*pi/soe; % number of node
        non = round(non);
        tin = pi/non;
        tet = 5*pi/2 : -tin : pi/2;
        tol = 0.000;
        while rem(length(tet),2) == 0
            soe = soe + tol; 
            non = r*pi/soe; 
            non = round(non);
            tin = pi/non;
            tet = 5*pi/2 : -tin : pi/2;
            tol = tol + 0.001;
        end
        xt2(1:length(tet)) = 0; 
        yt2(1:length(tet)) = 0;
        for i = 1 : length(tet)
            xt2(i) = r * cos(tet(i));
            yt2(i) = -r * sin(tet(i))-h;
        end
        xt = [xt1(1:end-1),xt2(1:end-1)];
        yt = [yt1(1:end-1),yt2(1:end-1)];
    case 2
        % ellipse inclusion
        % Part 1 (Outer boundary of inclusion)
        blim1 = 200;
        blim2 = 100;
        soe = 20;  % size of half-element
        non = ((pi/2)*(3*(blim1+blim2)-sqrt((3*blim1+blim2)*(blim1+3*blim2))))/soe; % number of node
        non = round(non);
        tin = pi/non;
        tet = pi/2 : tin : 5*pi/2;
        tol = 0.000;
        while rem(length(tet),2) == 0
            soe = soe + tol; 
            non = ((pi/2)*(3*(blim1+blim2)-sqrt((3*blim1+blim2)*(blim1+3*blim2))))/soe; 
            non = round(non);
            tin = pi/non;
            tet = pi/2 : tin : 5*pi/2;
            tol = tol + 0.001;
        end
        xt1(1:length(tet)) = 0; 
        yt1(1:length(tet)) = 0;
        for i = 1 : length(tet)
            xt1(i) = blim1 * cos(tet(i));
            yt1(i) = -blim2 * sin(tet(i))-h;
        end
        % Part 2 (Iner boundary of inclusion)
        soe = 20;  % size of half-element
        non = ((pi/2)*(3*(blim1+blim2)-sqrt((3*blim1+blim2)*(blim1+3*blim2))))/soe; % number of node
        non = round(non);
        tin = pi/non;
        tet = 5*pi/2 : -tin : pi/2;
        tol = 0.000;
        while rem(length(tet),2) == 0
            soe = soe + tol; 
            non = ((pi/2)*(3*(blim1+blim2)-sqrt((3*blim1+blim2)*(blim1+3*blim2))))/soe; 
            non = round(non);
            tin = pi/non;
            tet = 5*pi/2 : -tin : pi/2;
            tol = tol + 0.001;
        end
        xt2(1:length(tet)) = 0; 
        yt2(1:length(tet)) = 0;
        for i = 1 : length(tet)
            xt2(i) = blim1 * cos(tet(i));
            yt2(i) = -blim2 * sin(tet(i))-h;
        end
        xt = [xt1(1:end-1),xt2(1:end-1)];
        yt = [yt1(1:end-1),yt2(1:end-1)];      
end
figure(1)
plot(xt,yt)
axis equal
%----------------------prepare of input file-------------------------------
% 
input = 'input.txt';
input = lower(input);
oput = fopen(input,'wt');
fprintf(oput,'Ground_Surface_Motion_under_Vertically_SH_Waves  \n');
fprintf(oput,' %-10s\n',num2str(g1));
fprintf(oput,' %-10s\n',num2str(ro1));
fprintf(oput,' %-10s\n',num2str(g2));
fprintf(oput,' %-10s\n',num2str(ro2));
fprintf(oput,' %-10s\n',num2str(nstep));
fprintf(oput,' %-10s\n',num2str(delt));
nn = fix(length(yt));
ne = fix(nn/2);
radif = 1:length(yt);
fprintf(oput,' %-10s\n',num2str(ne));
nie = (fix(length(xt1(1:end-1))))/2;
fprintf(oput,' %-10s\n',num2str(nie));
for ii = 1 : nn
    fprintf(oput,'  %2d %8.2f %8.2f\n',[radif(ii);xt(ii);yt(ii)]);
end
fprintf(oput,' \n');
%----------------------------  prepare ricker file ------------------------
%
% ricker wavelet inputs
cs = sqrt(g1/ro1); % wave velocity
fp = 3; % predominant frequency
t0 = 2.4; % time shift parameter
%   total timet = 3 sec
dt = delt; % time interval
nt = nstep; % number of time steps
amax = 0.001; % maximum amplitude
teta = (00/180)*pi; % angle of incident wave
% prepare node coordinates file
noden = 'noden.txt';
noden = lower(noden);
node = fopen(noden,'wt');
fprintf(node,' %-10s\n',num2str(nn));
for ii = 1 : nn
    fprintf(node,'  %2d %8.2f %8.2f\n',[radif(ii);xt(ii);yt(ii)]);
end
fclose(node);
% prepare and write : ricker file 
r =zeros(1,length(xt1)-1);
rr =zeros(1,length(xt1)-1);
for ii = 1 : length(xt1)-1
    r(ii) = abs(-sin(teta)*xt1(ii)+cos(teta)*yt1(ii));
    rr(ii) = abs(-sin(teta)*xt1(ii)-cos(teta)*yt1(ii));
end       
t=zeros(1,nt);
for jj = 1 : nt
    if jj == 1
        t(jj) = 0;
    else
        t(jj) = t(jj-1)+dt;
    end
    landa = cs/fp;
    for kk = 1 : length(xt1)-1
        if r(kk) > cs*t(jj)
            eq1 = 0;
        else
            alpha_inc = cs*(t(jj)-t0)-sin(teta)*xt1(kk)+cos(teta)*yt1(kk);
            cst = ((pi/landa)*alpha_inc)^2;
            eq1 = -amax*(2*cst-1)*exp(-cst);
        end
        if rr(kk) > cs*t(jj)
            eq2 = 0;
        else
            alpha_ref = cs*(t(jj)-t0)-sin(teta)*xt1(kk)-cos(teta)*yt1(kk);
            cst = ((pi/landa)*alpha_ref)^2;
            eq2 = -amax*(2*cst-1)*exp(-cst);
        end
        EQ1(jj,kk)=eq1+eq2;
    end
end
find(rr==max(rr))
figure(2)
plot(t,EQ1(:,1),'b-')
% if seismic = 1 earthquake wave is exist else seimic is zero
seismic = 1;
% traction loading. 
nfun = 0;
%
fprintf(oput,' %-10s\n',num2str(seismic));
fprintf(oput,' %-10s\n',num2str(nfun));
% internal points
xi=5*blim:-20:-5*blim;
yi=zeros(1,length(xi));
intp = length(xi);
radif = 1:length(xi);
if intp ~= 0
    fprintf(oput,' %-10s\n',num2str(length(xi)));
    for ii = 1 : length(xi)
        fprintf(oput,'  %2d %8.2f %8.2f\n',[radif(ii);xi(ii);yi(ii)]);
    end
else
    fprintf(oput,' %-10s\n',num2str(intp));
end
fclose(oput);