function [BurstPnts,FollowerSpikeNum,NumSpikes,BurstISI] = BurstPntFinder(SpikePnts,QuietPnts,MinSpikeNumber) 

% This function finds the first spike which is preceded by a defined quiet
% period and followed by at least a minimum number of spikes

% BurstPnts = Cell structure of indicies which indicate the first spike times
% SpikePnts = Cell structure of indicies which indicate all spike times as
% found using "SpikeDetection"
% QuietPnts = Number of sample points without a spike
% MinSpikeNumber = Minimum number of spikes that define a burst

% JC 8/5/08

% edit JC 9/25/08 added NumSpikes and BurstISI
% NumSpikes = number of spikes in each selected burst
% BurstISI = the mean interspike interval for spikes within each burst


BurstPnts = Cell(size(SpikePnts)) ;                                 % prep cell structures for speed
FollowerSpikeNum = Cell(size(SpikePnts)) ; 
NumSpikes = Cell(size(SpikePnts)) ; 
BurstISI = Cell(size(SpikePnts)) ; 

for a= 1:length(SpikePnts) ;                                        % for each epoch of spikes...
    isi = diff(SpikePnts{a}) ;                                              % interval between each spike
    StarterSpikes = [1,find(isi>QuietPnts)+1] ;                             % indicies of SpikePnts that are preceded by a quiet time (the first spike is automatic)
    FollowerSpikeNum{a} = diff([StarterSpikes,(length(SpikePnts{a})+1)]) ;     % number of spikes between each starter spike
    TrueStarters = StarterSpikes(find(FollowerSpikeNum{a}>=MinSpikeNumber)) ;  % indicies of SpikePnts that are starters with more than one spike between them and the next stater
    BurstPnts{a} = SpikePnts{a}(TrueStarters) ;                                % indicies of raw data that are the first spikes in the defined burst
    
    NumSpikes{a} = FollowerSpikeNum{a}(find(FollowerSpikeNum{a}>=MinSpikeNumber)) ; % number of spikes in selected bursts
    NumPnts = SpikePnts{a}(TrueStarters+NumSpikes{a}-1) - BurstPnts{a} ;         % number of points between the first and last spike in each burst 
    BurstISI{a} = NumPnts./(NumSpikes{a}-1) ;                                    % the mean interspike interval for spikes within each burst
end

