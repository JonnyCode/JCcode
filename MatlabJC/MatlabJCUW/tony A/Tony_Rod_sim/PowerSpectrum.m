function ps = PowerSpectrum(x)
% ps = LogPowerSpectrum(x)
%
% This function returns the log of teh powerspectrum
% FFT is a complex function of frequency
% ps is given by:
% log(fftshift(x).*conj(fftshift(x)))
ps = (fftshift(x).*conj(fftshift(x)));