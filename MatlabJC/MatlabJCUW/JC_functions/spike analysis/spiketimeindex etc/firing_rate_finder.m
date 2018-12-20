function [freq_train, binned_freq_train]=firing_rate_finder(SpikeTimeIndex)

% this function will get take the spike Time index and get binned freqency 
                               
for a=1:length(SpikeTimeIndex);            % for each spike train
    inst_freq{a}=(1./diff(SpikeTimeIndex{a}))*1000;  % find the instantaneous frequency for each spike in Hz  

    SpikeTimeIndex2{a}=SpikeTimeIndex{a}(2:end)*10 ;  % remove the first spike in spike time index and change into .1ms from 1ms
   
    freq_train(a,:) = zeros(1,95000);             % prepare a place for instantaneous frequencies to go 
    freq_train(a,SpikeTimeIndex2{a})=inst_freq{a};   %place instantaneous frequencies at appropriate time

    freq_train2(a,:) = DecimateWave(freq_train(a,:), 100);         % averages each spike train over 10ms binns
end


binned_freq_train = mean(freq_train2) ;              % takes the mean the binned frq_trains

end