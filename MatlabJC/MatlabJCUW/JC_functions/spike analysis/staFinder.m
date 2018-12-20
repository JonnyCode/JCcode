function sta = staFinder(input,spikeTrain,staLength) 

% this function will get the average input waveform preceding a spike;
% spikeTrain must be a vector of 0s with 1s indicating spike.  staLength is
% the length of the desired spike triggered average.

% JC 9/22/11 
if size(input,1)~=size(spikeTrain,1) | size(input,2)~=size(spikeTrain,2) ;
    errormsg('size of input does not match size of spikeTrain')
end


spikeTrain(:,1:staLength) = 0 ; % spikes before length of sta cannot be used

totalNumSpikes = sum(spikeTrain(:)) ;
inputPreSpike = nan(totalNumSpikes,staLength) ;

[spikeRows,spikeColumns] = find(spikeTrain==1) ;
    
for a=1:length(spikeRows) ;
    inputPreSpike(a,:) = input(spikeRows(a),spikeColumns(a)-staLength+1:spikeColumns(a)) ;
end

sta = mean(inputPreSpike,1) ;
    