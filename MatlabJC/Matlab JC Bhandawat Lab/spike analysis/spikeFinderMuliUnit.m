function spikePnt = spikeFinderMuliUnit(data, sampleRate, NegDiffThreshStd, plotFigs)

% function will find spike times from scincillium multi unit recordings
% data is row vector, sampleRate in Hertz, NegDiffThreshStd is constant in
% delta mV, plotFigs is logical (true if you want a figure generated)
% JC 1/26/12
% JC edit 3/14/12 - changed NegDiffThresh into NegDiffThreshStd

smthTime = .0015 ; % sec over which data is smoothed 
HighPassFreq = 1 ; % hz high pass frequency 

smthPnts = sampleRate*smthTime ; 

dataHpf = highPassFilter(data,sampleRate,HighPassFreq) ; % high pass filter at 1 hz
dataHpfSmth = smooth(dataHpf,smthPnts) ; % average over approx maxRisePnts moving bin
dataHpfSmthDiff = diff(dataHpfSmth) ; % derivative

dataHpfSmthDiffStd = std(dataHpfSmthDiff(:)) ;
dataHpfSmthDiffMean = mean(dataHpfSmthDiff(:)) ;
NegDiffThresh = dataHpfSmthDiffMean - NegDiffThreshStd*dataHpfSmthDiffStd ;

[localMaxPnts,localMaxValues] = localMaxFinder(dataHpfSmth) ; % local maximum of data
localMaxPnts = [1,localMaxPnts] ; % include the very first point as a local max

spikePntTemp = nan(size(data)) ; % preallocate vector for speed
r=0 ;
for a=1:length(dataHpfSmthDiff) ; % for each diff data point
    if dataHpfSmthDiff(a)<NegDiffThresh ;
        r = r+1 ;
        localMaxi = find(localMaxPnts<=a+1,1,'last') ;
        spikePntTemp(r) = localMaxPnts(localMaxi) ;
    end
end    

spikePnt = unique(spikePntTemp(1:r)) ; % get rid of repeats and nans

% figures
if plotFigs ;
    figure
    plot(dataHpfSmthDiff)
    hold on 
    plot(ones(1,length(dataHpfSmthDiff))*NegDiffThresh,'k-')
    plot(ones(1,length(dataHpfSmthDiff))*dataHpfSmthDiffMean,'k:')    
    plot(ones(1,length(dataHpfSmthDiff))*(dataHpfSmthDiffMean - dataHpfSmthDiffStd),'k:')
end
