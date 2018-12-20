function [F,Fmean] = SpikeMatcherFast(SpikeTimeIndex, maxDeltaTVector) ;

% this function will calculate the fraction of spikes that can be matched
% for a given max time shift, for each possible pairs of spike trains.
% spiketimeindex is a cell array of spike times.  maxDeltaTVector, is a
% vector of the max time shifts.

% JC 4/28/09
F = nans(sum([1:length(SpikeTimeIndex)-1]),length(maxDeltaTVector)) ; % prepare matrix

pair = 0 ; 
for a = 1:length(SpikeTimeIndex) ; % for spike train
    for b = a+1:length(SpikeTimeIndex) ; % for each possible pairing
        pair = pair+1 ;
        D = abs(repmat(SpikeTimeIndex{a}',1,length(SpikeTimeIndex{b})) - repmat(SpikeTimeIndex{b},length(SpikeTimeIndex{a}),1)) ; % time difference between all spikes in two spike trains
        i = find(D<=maxDeltaTVector(end)) ;
        L = zeros(size(D)) ;
            for c = 1:length(maxDeltaTVector) ; % for each delta t
                L(i) = D(i)<= maxDeltaTVector(c) ;
                N = min(sum(sum(L,1)~=0),sum(sum(L,2)~=0)) ; % min(number of rows D<t, number of columns D<t) = number of possible matches
                F(pair,c) = N/max(length(SpikeTimeIndex{a}),length(SpikeTimeIndex{b})) ; % fraction of spikes matches
            end
    end
end

Fmean = mean(F,1) ; % mean of all pairs


% figure
% plot(maxDeltaTVector,F,'k')
% hold on
% plot(maxDeltaTVector,Fmean,'r')
