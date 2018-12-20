function [BurstStats,FirstBurstSpike_time,LastBurstSpike_time] = BurstFinder3(SpikeTimeIndex,ibi) ;

% This function will identify the bursts and then get
% rid of those bursts with only one spike and pull out relavent stats.  This is
% better than burstfinder_JC because it keeps the first and last bursts and
% doesn't choke if there is only one burst.

% JC 10/2/07


for a=1:length(SpikeTimeIndex) ;                        % for each spiketrain
    NumSpikes(a)=length(SpikeTimeIndex{a}) ;            % find the total number of spikes in that train
    isi=diff(SpikeTimeIndex{a}) ;                       % find the interspike interval
    between_bursts = find(isi>ibi) ;                    % if the isi is greater than the quiet time than it could be preceding a burst
    between_bursts2 = [0,between_bursts,NumSpikes(a)] ;   % this includes the first and last burst in analysis
    NumSpikesPerPseudoBurst = diff(between_bursts2) ;    % the number of spikes in each burst
    
    first = between_bursts2 + ones(1,length(between_bursts2)) ; % add 1 to find the spike indexs from the isi indexs
    first = first(1:end-1) ;                                    % get rid of first spike from last burst
    last = first + NumSpikesPerPseudoBurst ;                           % get the last spike
    last = last - ones(1,length(last)) ;                        % subtract 1 to find spike indexs from the isi indexs 
    
    burst_check = find(NumSpikesPerPseudoBurst>1) ;    % find index of those bursts with more than one spike
    FirstBurstSpike = first(burst_check) ;      % find first spikes index of bursts with one or more spikes
    LastBurstSpike = last(burst_check) ;        % find last spikes index of bursts with one or more spikes
    
    FirstBurstSpike_time{a} = SpikeTimeIndex{a}(FirstBurstSpike) ; %the spike time of the first spike in each burst
    LastBurstSpike_time{a} = SpikeTimeIndex{a}(LastBurstSpike) ;   % the spike time of the last spike in each burst
    
    NumBursts(a) = length(FirstBurstSpike) ;                        % find the number of bursts 
    NumSpikesPerBurst = NumSpikesPerPseudoBurst(burst_check) ;      % find the number of spikes in each burst 
    DurationBursts = LastBurstSpike_time{a} - FirstBurstSpike_time{a} ;   % find the duration of each burst
    BurstRate = NumSpikesPerBurst./DurationBursts*1000 ;            % firing rate of each bursts in Hz
    PercentBurstSpikes(a) = (sum(NumSpikesPerBurst))/NumSpikes(a) ;   % find precentage of spikes counted in a burst
    
    for b=1:length(FirstBurstSpike);                                                                    % for each burst 
        c=find(SpikeTimeIndex{a}>FirstBurstSpike_time{a}(b)&SpikeTimeIndex{a}<FirstBurstSpike_time{a}(b)+20); % find all the spikes that are within 20ms of the first spike in the burst
        InitBurstRate(b) = length(c)*50;                                                                % count the number of spikes and change to Hz
    end
    
    AdaptationMetric = InitBurstRate./BurstRate ;   % the higher this number the more adaptation within each burst
    
    % means within spike trains
    mean_DurationBursts(a) = mean(DurationBursts) ;         % mean duration of bursts per spike train
    mean_NumSpikesPerBurst(a) = mean(NumSpikesPerBurst) ;   % mean number of spikes per burst
    mean_InitBurstRate(a) = mean(InitBurstRate) ;           % mean firing rate in the first 20ms of a burst
    mean_BurstRate(a) = mean(BurstRate) ;                   % mean firing rate in the entire burst
    mean_AdaptationMetric(a) = mean(AdaptationMetric) ;     % mean adaptation metric per spike train
    
    clear DurationBursts NumSpikesPerBurst InitBurstRate BurstRate AdaptationMetric  
    
end % ends the SpikeTimeIndex spiketrain loop

% means between different spike trains in the same SpikeTimeIndex 
mean_NumSpikes = mean(NumSpikes) ;                      % mean number of spikes 
mean_NumBursts = mean(NumBursts) ;                      % mean number of bursts
mean_PercentBurstSpikes = mean(PercentBurstSpikes) ;    % mean percentage of spikes that are couted as part of a burst
mean2_DurationBursts = mean(mean_DurationBursts) ;      % mean of the mean duration of bursts 
mean2_NumSpikesPerBurst = mean(mean_NumSpikesPerBurst) ;% mean of the mean number of spikes per burst 
mean2_InitBurstRate = mean(mean_InitBurstRate) ;        % mean of the mean firing rate in the first 20ms of a burst
mean2_BurstRate = mean(mean_BurstRate) ;                % mean of the mean firing rate within a burst
mean2_AdaptationMetric = mean(mean_AdaptationMetric) ;  % mean of the mean adaptation metic for each burst

% standard deviation between different spike trains in the same SpikeTimeIndex
std_NumSpikes = std(NumSpikes) ;                            % standard deviation number of spikes 
std_NumBursts = std(NumBursts) ;                            % std number of bursts
std_PercentBurstSpikes = std(PercentBurstSpikes) ;          % std percentage of spikes that are couted as part of a burst
std_mean_DurationBursts = std(mean_DurationBursts) ;        % std of the mean duration of bursts within a spiketrain
std_mean_NumSpikesPerBurst = std(mean_NumSpikesPerBurst)    % std of the mean number of spikes per burst within a spiketrain 
std_mean_InitBurstRate = std(mean_InitBurstRate) ;          % std of the mean firing rate in the first 20ms of a burst
std_mean_BurstRate = std(mean_BurstRate) ;                  % std of the mean firing rate within a spiketrain
std_mean_AdaptationMetric = std(mean_AdaptationMetric) ;    % std of the mean adaptation metric within each spike train

BurstStats = [mean_NumSpikes, mean_NumBursts,mean_PercentBurstSpikes, mean2_NumSpikesPerBurst,mean2_DurationBursts, mean2_InitBurstRate, mean2_BurstRate,mean2_AdaptationMetric,...
    std_NumSpikes,std_NumBursts,std_PercentBurstSpikes,std_mean_NumSpikesPerBurst,std_mean_DurationBursts, std_mean_InitBurstRate,std_mean_BurstRate,std_mean_AdaptationMetric] ;


end % end the function

    
    