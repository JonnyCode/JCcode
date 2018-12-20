function [PSTH,SpikeTrain] = PSTHmaker(SpikePnts,SmoothPnts,RawDataLength,SampleRate) 

% This function uses SpikePnt indicies created by "SpikeDetection" and
% creates a spiketrain (1 and 0) and a smoothed PSTH

% Output:
% SpikeTrain - Matrix of zeros and ones which is the same length as the
% "RawDataLength" from which the spikes were originaly detected (each row
% is a trial).
% PSTH - "SpikeTrain" averaged over a moving window length of "SmoothPnts" 
% it is assumed that there are no spikes before or after the data.

% Input:
% SpikePnts - indicies of spikes in a cell structure, output of "SpikeDetection". 
% SmoothPnts - number of indicies over which you want to average.
% RawDataLength - the length of the raw data from which the spikes were
% detected. 
% samplerate - Hz (so that smoothPnts is time relative)

% JC 8/5/08

% edited JC 8/26/08 so that there is no time shift, the resulting
% PSTH is in units of spikes per second, and there is no issue if spike is
% at very begining or ending of spiketrain

SpikeTrain = zeros(length(SpikePnts),RawDataLength) ;  % prep SpikeTain matrix
PSTH = nans(length(SpikePnts),RawDataLength) ;         % prep PSTH matrix

for a = 1:length(SpikePnts) ;                   % for each set of spike points 
    SpikeTrain(a,SpikePnts{a}) = 1 ;                                            % indicate a one when there was a spike
    prePSTH = smooth([zeros(1,floor(SmoothPnts/2)),SpikeTrain(a,:),zeros(1,floor(SmoothPnts/2))],SmoothPnts) ;    % smooth SpikeTrain by SmoothPnts - the zero padding prevents moving average from accetuating early or late spikes (smooth filter begins as size 1, grows to 3, and keeps growing until it is size specified- it than shrinks back the same way)
    PSTH(a,:) = prePSTH((floor(SmoothPnts/2)+1):end-floor(SmoothPnts/2))*SampleRate ;                             % make the PSTH the correct length and in correct units
end


