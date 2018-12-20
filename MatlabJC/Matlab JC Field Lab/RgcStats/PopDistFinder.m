function PopDist = PopDistFinder(SpikeTimesArray, TrialStartTimes, varargin)

% this function will caculate a population vector distance between 1 time
% bin and all other time bins.  It will also find a population vector
% distance between the response and other trials for the 'same' stimulus.

% JCafaro 2/27/17

% parse inputs
p = inputParser;
p.addParamValue('PsthBinTime', 0.01, @isnumeric) ; % (sec) time of psth bin
p.addParamValue('BinSearchNumber', 0, @isnumeric) ; % number of psth bins beyond the time point that can be searched for distance minumum
p.addParamValue('TrialDuration', 0, @isnumeric) ; % (sec) duration of each trial
p.addParamValue('TrialNumberMax',0, @isnumeric) ; % number of rep trials used to calculate AcrossStim Distance

p.parse(varargin{:});
params = p.Results;
  
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
r1 = cell(1,trialNumber) ; % prep response vector
r2 = cell(1,trialNumber) ; % prep response vector

for tr = 1:trialNumber ; % for each trial
    Bins = [TrialStartTimes(tr):PsthBinTime:TrialStartTimes(tr)+TrialDuration] ; % time bins for trial 
    for c = 1:numCells ; % for each cell
        r{tr}(c,:)=histcounts(SpikeTimesArray{c},Bins) ; % number of spikes each bin
    end
    r_mean(tr,:) = mean(r{tr},1) ; % average number of spikes across all cells in the trial
end

PopDist.r_mean = mean(r_mean,1) ; % average spike number across cells and trials

% response discriminability
for tr = 1:trialNumberMax ; % for each trial to used to calculate Across Stim Distance
    for b = 1:numBins-1 ; % for each time bin 
        
        % distance between time bins within the same trial
        for b2 = 1:numBins-1 ; % for each time bin to be compared
            d_acrossStim(b2) = sqrt(sum((r{tr}(:,b)-r{tr}(:,b2)).^2)) ; % euclidean distance between vectors     
        end
        d_acrossStim(b) = nan ; % no distance from itself
        PopDist.AcrossStim(tr,b) = nanmean(d_acrossStim) ; % average across all bins
        
        % distance in same stimuli in different trials
        for tr2 = 1:trialNumber ; % for each trial to be compared
            d_temp = [] ; % prep mat
            for b2 = b-BinSearchNumber:b+BinSearchNumber ; % for each bin to be searched
                if b2>0 && b2<numBins ; 
                    d = sqrt(sum((r{tr}(:,b)-r{tr2}(:,b2)).^2)) ; % euclidean distance between vectors
                    d_temp = [d_temp,d] ;
                end
            end
            d_acrossTrials(tr2) = min(d_temp) ; % pic min distance within search range
        end
        d_acrossTrials(tr) = nan ; % no distance from itself
        PopDist.AcrossTrials(tr,b) = nanmean(d_acrossTrials) ; % average across all trials
        
        PopDist.AcrossStimMinusAcrossTrials(tr,b) = PopDist.AcrossStim(tr,b)-PopDist.AcrossTrials(tr,b) ;
        PopDist.AcrossStimOverAcrossTrials(tr,b) = PopDist.AcrossStim(tr,b)/PopDist.AcrossTrials(tr,b) ;
    end
end
        
                
             
            
            
        
            
            
        
        

