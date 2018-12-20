function hpFilteredSignal = highPassFilter(signal, samplerate, freqcutoff) 

% This function will high pass filter signals which are in rows in a
% matrix "signal", "samplerate" and "freqcutoff" should be points in Hz
% JC 6/11/08


% high pass filter currents to remove drift 
freqcutoff_adjusted = floor(freqcutoff/(samplerate/length(signal))) ; % this adjusts the freq cutoff for the length

signal_fft = fft(signal,[],2) ; % take fft of exc signal
signalFiltered_fft = signal_fft ;
signalFiltered_fft(:,[1:freqcutoff_adjusted,(length(signalFiltered_fft)-freqcutoff_adjusted):length(signalFiltered_fft)]) = 0 ; % cut out low frequencies in first and second half of fft
hpFilteredSignal = real(ifft(signalFiltered_fft,[],2)) ; % inverse fft

