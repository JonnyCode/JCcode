function[NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datarun,cellids,timedur, bin_size)



%Calculate average total number of spikes for each stimulus type      

%Input: Moving bar datarun structure

        %timedur: When last trial should stop
        
        %cellids: cells being analyzed

%Outputs: NumSpikesCell: Has average total number of spikes for each stimulus type for each cell

          %StimComb: All the DG stimulus combinations (col 1: spatial,col 2: temporal, col 3: direction, col 4: bar color, col 5: background color)

%Sneha Ravi 
%Last revision: 12-17-2012
%Revised by Xiaoyang Yao: 7-1-2015
          
          
%Initialize output matrices          
NumSpikesCell = zeros(length(cellids), length(datarun.stimulus.combinations));
MaxRate = zeros(length(cellids), length(datarun.stimulus.combinations));
StimComb = cell(length(datarun.stimulus.combinations), length(fieldnames(datarun.stimulus.params)));
StimComb(:,1) = {datarun.stimulus.combinations.BAR_WIDTH}';
StimComb(:,2) = {datarun.stimulus.combinations.DELTA}';
StimComb(:,3) = {datarun.stimulus.combinations.DIRECTION}';
StimComb(:,4) = {datarun.stimulus.combinations.RGB}';
StimComb(:,5) = {datarun.stimulus.combinations.BACK_RGB}';
edges = 0:bin_size:500;
%Calculation of total spike number for each trial for each cell, then average calculated and placed in NumSpikesCell
for k = 1:length(cellids)
    if ismember(cellids(k), datarun.cell_ids)
    i = get_cell_indices(datarun, cellids(1,k));
    NumSpikesAll = zeros(1, length(datarun.stimulus.triggers));
    MaxRateAll = zeros(1, length(datarun.stimulus.triggers));
    for j = 1: length(datarun.stimulus.triggers)
        if(j == length(datarun.stimulus.triggers))
%            Spikes_idx = datarun.spikes{i,1} >= datarun.stimulus.triggers(1,j) &datarun.spikes{i,1} < timedur;
             Spikes_idx = datarun.spikes{i,1} >= datarun.stimulus.triggers(j) &datarun.spikes{i,1} < timedur;
            Spikes_temp = datarun.spikes{i,1}(Spikes_idx);
%            Spikes_temp = Spikes_temp - datarun.stimulus.triggers(1,j);
             Spikes_temp = Spikes_temp - datarun.stimulus.triggers(j);
            NumSpikesAll(1, j) = sum(Spikes_idx);
            h = histc(Spikes_temp, edges);
            if isempty(h)
                h = zeros(size(h, 1),1);
            end
            MaxRateAll(1, j) = max(h)/bin_size;
        else
%            Spikes_idx = datarun.spikes{i,1} >= datarun.stimulus.triggers(1,j) &datarun.spikes{i,1} < datarun.stimulus.triggers(1,j+1);
             Spikes_idx = datarun.spikes{i,1} >= datarun.stimulus.triggers(j) &datarun.spikes{i,1} < datarun.stimulus.triggers(j+1);
            Spikes_temp = datarun.spikes{i,1}(Spikes_idx);
%            Spikes_temp = Spikes_temp - datarun.stimulus.triggers(1,j);
             Spikes_temp = Spikes_temp - datarun.stimulus.triggers(j);
            NumSpikesAll(1, j) = sum(Spikes_idx);
            h = histc(Spikes_temp, edges);
            if isempty(h)
                h = zeros(size(h, 1),1);
            end
            MaxRateAll(1, j) = max(h)/bin_size;
        end
    end
    NumSpikesCell(k,:) = grpstats(NumSpikesAll,datarun.stimulus.trial_list,'mean')';
    MaxRate(k,:) = grpstats(MaxRateAll,datarun.stimulus.trial_list,'mean')';
    end
end

end

%-------------------
%Slower way to do the same thing using cellfun
%for j = 1: (length(datarun.stimulus.triggers)-1)
%           NumSpikesAll{:,j} = cellfun(@(x) length(x(x >= datarun.stimulus.triggers(1,j) & x < datarun.stimulus.triggers(1,j+1))), datarun.spikes, 'UniformOutput', false);
%end

%NumSpikesAll{:,length(datarun.stimulus.triggers)} = cellfun(@(x) length(x(x >= datarun.stimulus.triggers(1,j) & x < timedur)), datarun.spikes, 'UniformOutput', false);
%-------------------