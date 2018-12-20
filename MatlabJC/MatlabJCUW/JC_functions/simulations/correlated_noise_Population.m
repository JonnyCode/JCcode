% populations of cells +/- correlated noise with similar and different linear filters

% jc 6/11/10

% parameters
numNeurons = 2 ; % number of simulation neurons 
time = [0:.0001:10] ;  % sec

tpeak = [.025, .070, .020, .100] ; % time to peak of filters
peakRise = [.009, .02, .030, .040] ; % rise time of filters
peakAmp = [1,1,1,1] ;
ttrough = [.060, .09, .070, .300] ; % trough time of filters
troughDecay = [.060, .09, .070, .300] ;
troughAmp = [1,.5,.5,.5] ;

signal = normrnd(0,20,1,length(time)) ; % create signal
signal = lowPassFilter(signal,100000,60) ;

% create filters
for neuron=1:numNeurons ; % for each neuron 
    filter(neuron,:) = simFilter(time,tpeak(neuron),peakRise(neuron),peakAmp(neuron),ttrough(neuron),troughDecay(neuron),troughAmp(neuron)) ;
end

% create an indepednent noise source 
for neuron=1:numNeurons ; % for each neuron
    noise(neuron,:) = normrnd(0,10,1,length(time)) ; 
    noise(neuron,:) = lowPassFilter(noise(neuron,:),100000,100) ;
end

for filterTrial = 1:2 ;

    for neuron=1:numNeurons ; % for each neuron add common or ind noise, filter stimulus, and estimate signal 
        if filterTrial == 1 ; % different filters
            response_ind(neuron,:) = ifft(fft(signal+noise(neuron,:)).*fft(filter(neuron,:))) ; % convolve filter and signal for response
            response_corr(neuron,:) = ifft(fft(signal+noise(1,:)).*fft(filter(neuron,:))) ;
        else % the same filters
            response_ind(neuron,:) = ifft(fft(signal+noise(neuron,:)).*fft(filter(1,:))) ; 
            response_corr(neuron,:) = ifft(fft(signal+noise(1,:)).*fft(filter(1,:))) ;
        end

        signalEst_ind(neuron,:) = ifft(fft(response_ind(neuron,:))./fft(filter(neuron,:))) ; % estimate stimulus
        signalEst_corr(neuron,:) = ifft(fft(response_corr(neuron,:))./fft(filter(neuron,:))) ;
    end

    signalEst_ind_mean = mean(signalEst_ind) ;
    signalEst_corr_mean = mean(signalEst_corr) ;

    signalEst_ind_MSE(filterTrial) = sum((signal-signalEst_ind_mean).^2) ;
    signalEst_corr_MSE(filterTrial) = sum((signal-signalEst_corr_mean).^2) ;

end

FracInd = signalEst_ind_MSE./signalEst_corr_MSE



figure
plot(time,filter)

figure
plot(time,signal)
hold on
plot(time,signalEst_ind_mean,'r')
plot(time,signalEst_corr_mean,'g')
