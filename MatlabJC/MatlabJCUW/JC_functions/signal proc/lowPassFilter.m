function lpFilteredSignal = lowPassFilter(signal, samplerate, freqcutoff) ;

% This function will low pass filter signals which are in rows in a
% matrix "signal", "samplerate" and "freqcutoff" should be points in Hz
% JC 3/31/08


% low pass filter currents to remove noise from instrumentation
freqcutoff_adjusted = floor(freqcutoff/(samplerate/length(signal))) ; % this adjusts the freq cutoff for the length

signal_fft = fft(signal,[],2) ; % take fft of exc signal
signalFiltered_fft = signal_fft ;
signalFiltered_fft(:,1+freqcutoff_adjusted:length(signal)-freqcutoff_adjusted) = 0 ; % cut out high frequencies in first and second half of fft
lpFilteredSignal = real(ifft(signalFiltered_fft,[],2)) ; % inverse fft

