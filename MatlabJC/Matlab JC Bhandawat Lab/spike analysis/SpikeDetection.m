function SpikePnts = SpikeDetection(Data,threshold,samplerate) ;

% spike detection for cell attached data.  This function will look for
% points that cross the upper threshold within 4/10ms of crossing the lower
% threshold. 

% input
% threshold = lower and upper bound
% samplerate = hz  
% data = data samples in rows

% output
% SpikePnts = each cell is contains indicies of detected spikes

% JC 6/11/08

PntTime = 1000/samplerate ; % ms per sample point 
hpData = highPassFilter(Data, samplerate, 20) ; % remove drift below 20Hz
SpikePnts = cell(1,size(Data,1)) ;              % each epoch has a cell of spikepnts which will be filled or empty

for epoch = 1:size(Data,1) ; % for each data sample ...

    count = 1 ; % for spike counting 
    
    % find the upper and lower cross points
    lowerPnts = find(hpData(epoch,:)<-threshold) ;                                    % all points below threshold
    upperPnts = find(hpData(epoch,:)>threshold) ;                                     % all points above threshold
    
    if ~isempty(lowerPnts) & ~isempty(upperPnts) ;                                  % as long as there are both points above and below threshold
        lowerCrossPnts = lowerPnts([find(diff(lowerPnts)>1),length(lowerPnts)]) ;       % the last points in a set of points below threshold      
        upperCrossPnts = upperPnts([1,(find(diff(upperPnts)>1)+1)]) ;                   % first points in a set of points above threshold
   
    % Spike if upper cross point is preceded by lower cross point within 4/10 ms   
    for a = 1:length(upperCrossPnts) ; % for each upper cross point ...
        
        if ~isempty(find(lowerCrossPnts < upperCrossPnts(a) & lowerCrossPnts > upperCrossPnts(a) - ((4/10)/PntTime))) ;
            SpikePnts{epoch}(count) = upperCrossPnts(a) ;
            count = count + 1 ;
        end        
    end 
    end                                   
        
end

end % ends function

