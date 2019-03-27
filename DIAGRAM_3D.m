% 3D diagram
global X CS1 L AT XI NIE
global U UI EQ1 EQ NSTEP
%
NIN = 2*NIE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coding = 1 for time-domain diagram & coding = 2 for freqency-domain diagram
coding = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if coding == 1
    tdr = U;
    tdrint = UI;
    
    tdr_tot= tdrint(:,:);
    eq_tot = EQ(:,:)';
    xnew = XI(:,:);
    time = AT:AT:NSTEP*AT;
    figure(1)
    hold on
    for i=1:1:length(xnew)
        plot3(time,xnew(i)*ones(1,length(time))/500,tdr_tot(i,:)/0.001,'k')
    end
    view(0,83)
    [tt,xxx] = meshgrid(time,xnew/X(1));
    mesh(tt,xxx,tdr_tot./0.001)
    return
end
nfft = 2^nextpow2(3000);
fftsignal = zeros(NIN,nfft);
fftinter = zeros(L,nfft);
ffteq1 = zeros(NIN,nfft);
ffteq = zeros(L,nfft);
for i = 1 : NIN
    fftsignal(i,:) = fft(U(i,:),nfft);
    ffteq1(i,:) = fft(EQ1(:,i),nfft);
end
for i = 1 : L
    fftinter(i,:) = fft(UI(i,:),nfft);
    ffteq(i,:) = fft(EQ(:,i),nfft);
end
% correction of sem-sine shaped valley at corner points
%Next, calculate the frequency axis, which is defined by the sampling rate
T = AT;
fs = 1/T;
f = fs/2*linspace(0,1,nfft/2+1);
% dimensionless frequency
df = 2*500*f./CS1;
% select dimensionless frequency " max "
sdf_max = 5;
sdf_min = 0.35;
vec_max = find(df>=sdf_max);
vec_min = find(df<=sdf_min);
frq_max = vec_max(1);
frq_min = vec_min(end);
dff = df(frq_min:frq_max);
%
xnew = XI;
%
[ff,xx] = meshgrid(dff,xnew/500);
%
ffteqn = zeros(NIN,nfft);
for i=1:NIN
    ffteqn(i,:) = ffteq1(end,:);
end
ynew = abs(fftinter(1:L,frq_min:frq_max))./(abs(ffteq(1:L,frq_min:frq_max))/1);
figure(2)
mesh(ff,xx,ynew)
view(-70,80)
axis([sdf_min sdf_max -5 5 min(min(ynew)) max(max(ynew))])
colorbar