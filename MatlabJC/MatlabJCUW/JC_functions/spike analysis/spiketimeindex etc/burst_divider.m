function burst_spikes = burst_divider(SpikeTimeIndex, ibi);

% ibi = interburst interval
% burst_divider will divide a set of spike times into a matrix of
% spiketimes with each burst (which can be only one spike) in its own row 
% J. Cafaro 6/25/07

for a=1:length(SpikeTimeIndex) ;              % for each spiketrain
    num_spikes=length(SpikeTimeIndex{a}) ;    % find the total number of spikes in that train
    isi=diff(SpikeTimeIndex{a}) ;            % find the interspike interval
    between_bursts = find(isi>ibi) ;   % if the isi is greater than the quiet time than it could be preceding a burst
    between_bursts2 = [0,between_bursts,num_spikes] ; % this includes the first and last burst in analysis
    num_burst_spikes = diff(between_bursts2) ;  % the number of spikes in each burst
    
    first = between_bursts2 + ones(1,length(between_bursts2)) ;       % add 1 to find the spike indexs from the isi indexs
    first = first(1:end-1) ;                       % get rid of first spike from last burst
    last = first + num_burst_spikes ;             % get the last spike
    last = last - ones(1,length(last)) ;            % subtract 1 to find spike indexs from the isi indexs                             
    
    burst_spikes{a} = NaNs(length(first),max(num_burst_spikes)) ; % make a NaN matrix to put in burst spike times
    for b = 1:length(first) ;% for each burst... 
        burst_spikes{a}(b,1:num_burst_spikes(b)) = SpikeTimeIndex{a}(first(b):last(b)) ; % put each spike time in a burst into a row 
    end
end
   
    