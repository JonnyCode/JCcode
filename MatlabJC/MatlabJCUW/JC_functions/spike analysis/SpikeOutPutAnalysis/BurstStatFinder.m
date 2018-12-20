function [BurstPnts,LoneSpkPnts,NumSpikes,BurstISI,FracBurstSpikes] = BurstStatFinder(SpikePnts,QuietPnts,MinSpikeNumber) 

% This function finds the first spike which is preceded by a defined quiet
% period and followed by at least a minimum number of spikes

% BurstPnts = Cell structure of indicies which indicate the first spike times
% SpikePnts = Cell structure of indicies which indicate all spike times as
% found using "SpikeDetection"
% QuietPnts = Number of sample points without a spike
% MinSpikeNumber = Minimum number of spikes that define a burst

% JC 8/5/08

% edit JC 9/25/08 added stats for bursts below
% NumSpikes = number of spikes in each selected burst
% BurstISI = the mean interspike interval for spikes within each burst
% FracBurstSpikes = fraction of spikes in the epoch that are counted as part of a burst

% edit JC 10/26/08 added id of non burst spikes

BurstPnts = cell(size(SpikePnts)) ;                                 % prep cell structures for speed
FollowerSpikeNum = cell(size(SpikePnts)) ; 
NumSpikes = cell(size(SpikePnts)) ; 
BurstISI = cell(size(SpikePnts)) ; 
FracBurstSpikes = nans(1,length(SpikePnts)) ;

for a= 1:length(SpikePnts) ;                                        % for each epoch of spikes...
    isi = diff(SpikePnts{a}) ;                                              % interval between each spike
    StarterSpikes = [1,find(isi>QuietPnts)+1] ;                             % indicies of SpikePnts that are preceded by a quiet time (the first spike is automatic)
    FollowerSpikeNum{a} = diff([StarterSpikes,(length(SpikePnts{a})+1)]) ;      % number of spikes between each starter spike
    TrueStarters = StarterSpikes(find(FollowerSpikeNum{a}>=MinSpikeNumber)) ;   % indicies of SpikePnts that are starters with more than one spike between them and the next starter
    FalseStarters = StarterSpikes(find(FollowerSpikeNum{a}<MinSpikeNumber)) ;   % indicies of SpikePnts that are spikes outside of a defined burst
    BurstPnts{a} = SpikePnts{a}(TrueStarters) ;                                 % indicies of raw data that are the first spikes in the defined burst
    LoneSpkPnts{a} = SpikePnts{a}(FalseStarters) ;                              % indicies of raw data that are spikes outside of a defined burst
    
    NumSpikes{a} = FollowerSpikeNum{a}(find(FollowerSpikeNum{a}>=MinSpikeNumber)) ; % number of spikes in selected bursts
    NumPnts = SpikePnts{a}(TrueStarters+NumSpikes{a}-1) - BurstPnts{a} ;            % number of points between the first and last spike in each burst 
    BurstISI{a} = NumPnts./(NumSpikes{a}-1) ;                                       % the mean interspike interval for spikes within each burst
    FracBurstSpikes(a) = sum(NumSpikes{a})/length(SpikePnts{a}) ;                   % the fraction of spikes in the epoch that are counted as part of a burst
end

