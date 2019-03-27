% Calculate FFT
global U AT NN NIE EQ1 CS1 UI XI EQ L
NIN = 2*NIE;
nfft = 2^nextpow2(5000);
fftsignal = zeros(NN,nfft);
fftinter = zeros(L,nfft);
ffteq1 = zeros(NIN,nfft);
ffteq = zeros(L,nfft);
for i = 1 : NN
    fftsignal(i,:) = fft(U(i,:),nfft);
end
for i = 1:NIN
    ffteq1(i,:) = fft(EQ1(:,i),nfft);
end
for i = 1 : L
    fftinter(i,:) = fft(UI(i,:),nfft);
    ffteq(i,:) = fft(EQ(:,i),nfft);
end
%Next, calculate the frequency axis, which is defined by the sampling rate
T = AT;
fs = 1/T;
f = fs/2*linspace(0,1,nfft/2+1);
% dimensionless frequency
df = 2*500*f./CS1;
% select dimensionless frequency
sdf = 1;
vect = find(df>=sdf);
frq = vect(1);
plotxsig = zeros(1,NN);
plotxeq1 = zeros(1,NIN);
plotxint = zeros(1,L);
plotxeq = zeros(1,L);
for i = 1 : NN
    plotxsig(i) = fftsignal(i,frq);
end
plotxsig(NN+1)=plotxsig(NIN+1);
for i = 1 : NIN
    plotxeq1(i) = ffteq1(i,frq);
end
for i = 1 : L
    plotxint(i) = fftinter(i,frq);
    plotxeq(i) = ffteq(i,frq);
end
figure(1)
hold on
plot(XI./500,(abs(plotxint)./(abs(plotxeq(1))/2)),'k-') % ground surface