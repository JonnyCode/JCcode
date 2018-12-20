function[NumSpikesCell, StimComb] = get_spikescellstimJC(datarun,cellids,totaldur)

% adapted from 'get_spikescellstimJC' to look at statistics other than mean spike number
% JC 11/1/2016

%Calculate average total number of spikes for each stimulus type      

%Input: Drifting grating datarun structure

        %totaldur: Time over which to average spikes: 1 - average over entire trial, 0 - average over DG stimulus
        
        %cellids: cells being analyzed

%Outputs: NumSpikesCell: Has average total number of spikes for each stimulus type for each cell

          %StimComb: All the DG stimulus combinations (col 1: spatial,col 2: temporal, col 3: direction)

%Sneha Ravi 
%Last revision: 12-17-2012
          
          
%Initialize output matrices          
NumSpikesCell = zeros(length(cellids), length(datarun.stimulus.combinations));
StimComb = zeros(length(datarun.stimulus.combinations), length(fieldnames(datarun.stimulus.params)));
StimComb(:,1) = [datarun.stimulus.combinations.SPATIAL_PERIOD];
StimComb(:,2) = [datarun.stimulus.combinations.TEMPORAL_PERIOD];
StimComb(:,3) = [datarun.stimulus.combinations.DIRECTION];

%When last trial should stop
if(totaldur == 1)
    timedur = datarun.stimulus.triggers(1,length(datarun.stimulus.triggers)) + mean(diff(datarun.stimulus.triggers));
end
if(totaldur == 0)
    timestop = 8; %8s stimulus run
    timedur = datarun.stimulus.triggers(1,length(datarun.stimulus.triggers)) + timestop; %8s stimulus run
end

%Calculation of total spike number for each trial for each cell, then average calculated and placed in NumSpikesCell
for k = 1:length(cellids)
    i = get_cell_indices(datarun, cellids(1,k));
    NumSpikesAll = zeros(1, length(datarun.stimulus.triggers));
    [psth,psthTime] = get_smooth_psth(spike_times, datarun.stimulus.triggers, varargin)
    NumSpikesCell(k,:) = grpstats(NumSpikesAll,datarun.stimulus.trial_list,'mean')';
end

end

%-------------------
%Slower way to do the same thing using cellfun
%for j = 1: (length(datarun.stimulus.triggers)-1)
%           NumSpikesAll{:,j} = cellfun(@(x) length(x(x >= datarun.stimulus.triggers(1,j) & x < datarun.stimulus.triggers(1,j+1))), datarun.spikes, 'UniformOutput', false);
%end

%NumSpikesAll{:,length(datarun.stimulus.triggers)} = cellfun(@(x) length(x(x >= datarun.stimulus.triggers(1,j) & x < timedur)), datarun.spikes, 'UniformOutput', false);
%-------------------