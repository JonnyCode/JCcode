function [spikePnt,SpikeData,NonSpikeData] = spikeFinder(data,sampleRate,spikeDetectionParameters)

% this function will look for spikes that are the largest local max,
% seperated by absRefTime, with a possible plateau, preceded a minRise and followed by a minFall 

%data in rows, sampeRate in Hz, minRise and minFall in mV,
%absRefTime in seconds, spikePnt in sample points

% JC 5/23/12 (earlier version of function with same name is entirely different)
% JC 12/14/12 revised to add spike detection params to input params (changed into a structure)

% hard coded parameters
absRefTime = spikeDetectionParameters.absRefTime ; % sec
minRise = spikeDetectionParameters.minRise ; % mV
minFall = spikeDetectionParameters.minFall ; % mV
filterOrder = spikeDetectionParameters.filterOrder ;
lpfCutOff = spikeDetectionParameters.lpfCutOff ; % hz
minRiseTime = spikeDetectionParameters.minRiseTime ; % sec
minFallTime = spikeDetectionParameters.minFallTime ; % sec
maxPlatTime = spikeDetectionParameters.maxPlatTime ; % sec 

% pnts from time
absRefPnts = round(absRefTime*sampleRate) ;
minRisePnts = round(minRiseTime*sampleRate) ;
minFallPnts = round(minFallTime*sampleRate) ;
maxPlatPnts = round(maxPlatTime*sampleRate) ;

% implement a butterworth to get rid of high frequ noise
Wn = lpfCutOff*2/sampleRate; %normalized frequency cutoff (fraction of nyquist)
[z, p, k] = butter(filterOrder,Wn,'low'); %high can be changed to low
[sos,g]=zp2sos(z,p,k); % convesion?
myfilt=dfilt.df2sos(sos,g);
dataLpf = filter(myfilt,data')'; %implementation

% all local maximum 
spikePntTemp = cell(1,size(data,1)) ;
for a=1:size(data,1) ;
    [pks,locs] = findpeaks(dataLpf(a,:),'minpeakdistance',absRefPnts) ;
    spikePntTemp{a} = locs ;
end

% are local max truly spikes?
spikePnt = cell(1,size(data,1)) ;
nonspikePnt = cell(1,size(data,1)) ;
numSpikes = nan(1,size(data,1)) ;
numNonSpikes = nan(1,size(data,1)) ;

for a=1:length(spikePntTemp) ;
    for b=1:length(spikePntTemp{a}) ;
        if spikePntTemp{a}(b)>minRisePnts+maxPlatPnts && spikePntTemp{a}(b)<length(data)-minFallPnts-maxPlatPnts % if theres enough space to assess spike
            if sum(dataLpf(a,spikePntTemp{a}(b)) - dataLpf(a,spikePntTemp{a}(b)-minRisePnts-maxPlatPnts:spikePntTemp{a}(b)-minRisePnts)>minRise)~=0 && ... % if the riseing slope is steep enough
            sum(dataLpf(a,spikePntTemp{a}(b)) - dataLpf(a,spikePntTemp{a}(b)+minFallPnts:spikePntTemp{a}(b)+minFallPnts+maxPlatPnts)>minFall)~=0 ; % if the falling slope is steep enough        
                
                spikePnt{a} = [spikePnt{a},spikePntTemp{a}(b)] ; % then it is a spike
            else
                nonspikePnt{a} = [nonspikePnt{a},spikePntTemp{a}(b)] ; % then its not a spike
            end
        end
    end
    numSpikes(a) = length(spikePnt{a}) ;
    numNonSpikes(a) = length(nonspikePnt{a}) ;
end

% ensemble of selected and unselected spikes (local max with abs ref between)
NonSpikeData = nan(sum(numNonSpikes),minFallPnts+minRisePnts+2*maxPlatPnts+1) ;
SpikeData = nan(sum(numSpikes),minFallPnts+minRisePnts+2*maxPlatPnts+1) ;

nonspikeR=0 ;
spikeR=0 ;

for a=1:length(spikePntTemp) ;
    for b=1:length(spikePnt{a}) ;
        spikeR = spikeR+1 ;
        SpikeData(spikeR,:) = dataLpf(a,spikePnt{a}(b)-minFallPnts-maxPlatPnts:spikePnt{a}(b)+minRisePnts+maxPlatPnts) ;
        SpikeData(spikeR,:) = SpikeData(spikeR,:)-SpikeData(spikeR,1) ;
    end
    
    for b=1:length(nonspikePnt{a}) ;
        nonspikeR = nonspikeR+1 ;
        NonSpikeData(nonspikeR,:) = dataLpf(a,nonspikePnt{a}(b)-minFallPnts-maxPlatPnts:nonspikePnt{a}(b)+minRisePnts+maxPlatPnts) ;
        NonSpikeData(nonspikeR,:) = NonSpikeData(nonspikeR,:)-NonSpikeData(nonspikeR,1) ;
    end    
end





