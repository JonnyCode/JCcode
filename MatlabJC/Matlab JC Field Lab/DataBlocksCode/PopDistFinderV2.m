function PopDist = PopDistFinderV2(SpikeTimesArray, TrialStartTimes, varargin)

% adapted from 'PopDistFinder' to utilize z score distance (calc across
% trials).

% JCafaro 2/27/17

% parse inputs
p = inputParser;
p.addParamValue('PsthBinTime', 0.01, @isnumeric) ; % (sec) time of psth bin
p.addParamValue('BinSearchNumber', 0, @isnumeric) ; % number of psth bins beyond the time point that can be searched for distance minumum
p.addParamValue('TrialDuration', 0, @isnumeric) ; % (sec) duration of each trial
p.addParamValue('TrialNumberMax',0, @isnumeric) ; % number of rep trials used to calculate AcrossStim Distance

p.parse(varargin{:});
params = p.Results;

minStd = .001 ; % minimum std can be 

if params.TrialDuration == 0; % if no length of trial has been specified
    TrialDuration = min(diff(TrialStartTimes)) ; % assume its the time between triggers
else
    TrialDuration = params.TrialDuration ;
end

if params.TrialNumberMax == 0; % if no number of trial has been specified
    trialNumberMax = length(TrialStartTimes) ; % use all trials
else
    trialNumberMax = params.TrialNumberMax ;
end

PsthBinTime = params.PsthBinTime ;
BinSearchNumber = params.BinSearchNumber ;

trialNumber = length(TrialStartTimes) ; % number of trials
numCells = length(SpikeTimesArray) ; % number of cells
numBins = floor(TrialDuration/PsthBinTime) ; % number of time bins

% responses 
for c = 1:numCells ; % for each cell
    for tr = 1:trialNumber ; % for each trial
        Bins = [TrialStartTimes(tr):PsthBinTime:TrialStartTimes(tr)+TrialDuration] ; % time bins for trial 
        r{c}(tr,:)=histcounts(SpikeTimesArray{c},Bins) ; % number of spikes each bin
    end
    r_mean(c,:) = mean(r{c},1) ; % average spike number across trials in a bin
    r_std(c,:) = std(r{c},[],1)+minStd ; % std of spike number across traisl in a bin
    r_mean_mean(c) = mean(r_mean(c,:)) ; % mean spike number across bins
    
    z(c,:) = (r_mean(c,:)-r_mean_mean(c))./(r_std(c,:)) ; % number of std bin mean deviates from mean of all bins
end

PopDist.r_mean = r_mean ; % psth(cell,time bin)

% response discriminability

for b = 1:numBins ; % for each time bin 
    % distance between time bins 
    for b2 = 1:numBins ; % for each time bin to be compared
        %d_acrossStim(b2) = sqrt(sum((z(:,b)-z(:,b2)).^2)) ; % euclidean distance between vectors     
        d_acrossStim(b2) = abs(r_mean(:,b)-r_mean(:,b2))/(sqrt((r_std(:,b).^2+r_std(:,b2).^2)/2)) ; % d'
    end
    d_acrossStim(b) = nan ; % no distance from itself
    PopDist.AcrossStim(b) = nanmean(d_acrossStim) ; % average across all bins
end

        
                
             
            
            
        
            
            
        
        

