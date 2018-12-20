function ps = LogPowerSpectrum(x)
% ps = LogPowerSpectrum(x)
%
% This function returns the log of teh powerspectrum
% FFT is a complex function of frequency
% ps is given by:
% log(fftshift(x).*conj(fftshift(x)))
ps = log(fftshift(x).*conj(fftshift(x)));