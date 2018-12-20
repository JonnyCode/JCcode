function Inew = addSpatialPhaseNoise(I, freqInterval, noiseAmp)
%freqInterval is an interval in log space in [0 1]
%for now, I must be square
%noise amp is a fraction of 2pi for s.d. of noise added

fI = fft2(I);
%fI = fftshift(fI);
R = abs(fI);
theta = angle(fI);

L = size(I,1); %L is also size of fI

logBins = logspace(log10(1),log10(L),100);
freqInterval = round(logBins(ceil(freqInterval*100)));

FreqLow = freqInterval(1);
FreqHigh = freqInterval(2);

%add noise to specific frequencies
Ind1 = FreqLow:FreqHigh;
Ind2 = L-[FreqLow:FreqHigh]+1;

freqBins = length(Ind1);

noiseAmp = noiseAmp*2*pi;

theta(Ind1,Ind1) = theta(Ind1,Ind1) + noiseAmp*randn(freqBins,freqBins);
theta(Ind2,Ind2) = theta(Ind2,Ind2) + noiseAmp*randn(freqBins,freqBins);

fInew = R.*exp(1i*theta);
%fInew = ifftshift(fInew);
Inew = abs(ifft2(fInew));
