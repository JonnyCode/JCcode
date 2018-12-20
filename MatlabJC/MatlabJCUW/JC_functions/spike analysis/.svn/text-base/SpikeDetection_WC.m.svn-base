function SpikePnts = SpikeDetection_WC(Data,threshold,samplerate) ;

% spike detection for current clamp data.  This function will look for
% points that cross an upper threshold more than 1 ms between spikes. 

% input
% threshold = upper bound
% samplerate = hz  
% data = data samples in rows

% output
% SpikePnts = each cell is contains indicies of detected spikes as labeled
% by at the first point above threshold

% JC 6/23/08

PntTime = 1000/samplerate ; % ms per sample point 
SpikePnts = cell(1,size(Data,1)) ;              % each epoch has a cell to be filled with spike indicies or left empty if there are no spikes

for epoch = 1:size(Data,1) ; % for each data sample ...
                       
    upperPnts = find(Data(epoch,:)>threshold) ;                                     % all points above threshold
    if ~isempty(upperPnts) ;                                                        % as long as there are points above threshold
        upperCrossPnts = upperPnts([1,(find(diff(upperPnts)>1)+1)]) ;                           % first points in a set of points above threshold
        SpikePnts{epoch} = upperCrossPnts([1,find(diff(upperCrossPnts)>(1/PntTime))+1]) ;       % pick only those first points that are preceded by at least 1ms of no other idetified first points
    end % if spikes exist loop
end % epoch loop
    

end % ends function
