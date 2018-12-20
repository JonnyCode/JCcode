function [powerspec_xvalues, mean_powerspec] = PowerSpectrumFinder(signal,samplerate) ;

% this function will find the power spectrum of a signal
% input = signal in matrix with individual signals in a row, samplerate of
% signal in Hz (eg - 10000)
% output = xvalues for power spectrum, powerspectrum

% J Cafaro 9/24/07

points = length(signal) ; % number of sample points
maxfreq = samplerate/2 ; % Nyquist
powerspec_xvalues = [0:2*maxfreq/points:maxfreq]; % time vector for power spectrum

fft_signal = fft(signal,[],2) ;     % fourier transform of signal
tempps = fft_signal.*conj(fft_signal) ;  % power spectrum (almost)
tempps2 = real(tempps) ;   %take only real parts
powerspec = (2*tempps2*(1/samplerate))/points ; % power spectrum for each individual Gexc
mean_powerspec = mean(powerspec(:,1:length(powerspec_xvalues)),1) ;  % mean power spectrum 



%NewPowerSpec = SmoothPowerSpectrum(mean_powerspec, OriginalFreq, SmoothFact, SkipPts) % function to smooth powerspec for log scaling
