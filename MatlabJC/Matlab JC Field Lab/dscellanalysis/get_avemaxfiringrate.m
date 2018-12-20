function[NumSpikesCell, StimComb] = get_avemaxfiringrate(datarun)

NumSpikesCell = zeros(length(datarun.cell_ids), length(datarun.stimulus.combinations));
StimComb = zeros(length(datarun.stimulus.combinations), length(fieldnames(datarun.stimulus.params))); 
%timedur = datarun.stimulus.triggers(1,length(datarun.stimulus.triggers)) + mean(diff(datarun.stimulus.triggers));
timedur = datarun.stimulus.triggers(1,length(datarun.stimulus.triggers)) + 8;%8s stimulus run
timestop = 8; %8s stimulus run
%nbins = mean(diff(datarun.stimulus.triggers)); %??? bin size now = 1s
nbins = 8; %8s stimulus run
binsize = 1; % 1 second

StimComb(:,1) = [datarun.stimulus.combinations.SPATIAL_PERIOD];
StimComb(:,2) = [datarun.stimulus.combinations.TEMPORAL_PERIOD];
StimComb(:,3) = [datarun.stimulus.combinations.DIRECTION];



for i = 1:length(datarun.cell_ids)
    NumSpikesAll = zeros(1, length(datarun.stimulus.triggers));
    for j = 1: length(datarun.stimulus.triggers)
        if(j == length(datarun.stimulus.triggers))
            [nhist xouthist] = hist(datarun.spikes{i,1}(datarun.spikes{i,1} >= datarun.stimulus.triggers(1,j) &datarun.spikes{i,1} < timedur), nbins);
            NumSpikesAll(1, j) = max(nhist);
        else
            %[nhist xouthist] = hist(datarun.spikes{i,1}(datarun.spikes{i,1} >= datarun.stimulus.triggers(1,j) &datarun.spikes{i,1} < datarun.stimulus.triggers(1,j+1)), nbins);
            [nhist xouthist] = hist(datarun.spikes{i,1}(datarun.spikes{i,1} >= datarun.stimulus.triggers(1,j) &(datarun.spikes{i,1} < datarun.stimulus.triggers(1,j) + timestop)), nbins); %8s stimulus run
            NumSpikesAll(1, j) = max(nhist);
        end
        
    end
    NumSpikesCell(i,:) = grpstats(NumSpikesAll,datarun.stimulus.trial_list,'mean')';
end

