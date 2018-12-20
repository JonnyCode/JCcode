function Inew = addSpatialAmplitudeNoise(I, freqInterval, noiseAmp)
%freqInterval is an interval in log space in [0 1]
%for now, I must be square
%noise amp is a fraction of variance across all freqs

fI = fft2(I);
%fI = fftshift(fI);
R = abs(fI);
theta = angle(fI);

L = size(I,1); %L is also size of fI

logBins = logspace(log10(1),log10(L),100);
freqInterval = round(logBins(ceil(freqInterval*100)));

%FreqStepSize = 1/L;
%FreqLow = ceil(freqInterval(1) / FreqStepSize);
%FreqHigh = ceil(freqInterval(2) / FreqStepSize);
FreqLow = freqInterval(1);
FreqHigh = freqInterval(2);

%add noise to specific frequencies
Ind1 = FreqLow:FreqHigh;
Ind2 = L-[FreqLow:FreqHigh]+1;

freqBins = length(Ind1);

%noisePerFreq = (sum(sum(R(Ind1,Ind1)))+sum(sum(R(Ind2,Ind2)))).*noiseAmp./freqBins;
noisePerFreq = sum(sum(R(Ind1,Ind1))).*noiseAmp./freqBins.^2;

R(Ind1,Ind1) = R(Ind1,Ind1)+noisePerFreq*randn(freqBins,freqBins);
R(Ind2,Ind2) = R(Ind2,Ind2)+noisePerFreq*randn(freqBins,freqBins);
fInew = R.*exp(1i*theta);
%fInew = ifftshift(fInew);
Inew = abs(ifft2(fInew));
