function [LinearFilter] = LinFilterFinder(signal,response, samplerate, freqcutoff)

% this function will find the linear filter that changes row vector "signal" 
% into a set of "responses" in rows.  "samplerate" and "freqcuttoff" 
% (which should be the highest frequency in the signal) should be in HZ.

% The linear filter is a cc normalized by the power spectrum of the signal
% JC 3/31/08 

% JC 8/18/10 added floor to freqcuttoff_adjusted

FilterFft = mean((fft(response,[],2).*conj(fft(signal,[],2))),1)./mean(fft(signal,[],2).*conj(fft(signal,[],2)),1) ;

freqcutoff_adjusted = floor(freqcutoff/(samplerate/length(signal))) ; % this adjusts the freq cutoff for the length

FilterFft(:,1+freqcutoff_adjusted:length(signal)-freqcutoff_adjusted) = 0 ; 

LinearFilter = real(ifft(FilterFft)) ;

end