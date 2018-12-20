function[NumSpikesCell, StimComb] = get_maxavefiringrate(datarun)


NumSpikesCell = zeros(length(datarun.cell_ids), length(datarun.stimulus.combinations));
StimComb = zeros(length(datarun.stimulus.combinations), length(fieldnames(datarun.stimulus.params))); 
%timedur = datarun.stimulus.triggers(1,length(datarun.stimulus.triggers)) + mean(diff(datarun.stimulus.triggers));
timedur = datarun.stimulus.triggers(1,length(datarun.stimulus.triggers)) + 8;%8s stimulus run
timestop = 8; %8s stimulus run
%nbins = round(mean(diff(datarun.stimulus.triggers))); %??? bin size now = 1s
nbins = 8; %8s stimulus run
%binsize = mean(diff(datarun.stimulus.triggers)) / nbins; % 1 second
binsize = 0.5; % 1 second, 8s stimulus run




StimComb(:,1) = [datarun.stimulus.combinations.SPATIAL_PERIOD];
StimComb(:,2) = [datarun.stimulus.combinations.TEMPORAL_PERIOD];
StimComb(:,3) = [datarun.stimulus.combinations.DIRECTION];



for i = 1:length(datarun.cell_ids)
    nhist = zeros(length(datarun.stimulus.triggers), nbins);
    for j = 1: length(datarun.stimulus.triggers)
        if(j == length(datarun.stimulus.triggers))
            [nhist(j,:) xouthist] = hist(datarun.spikes{i,1}(datarun.spikes{i,1} >= datarun.stimulus.triggers(1,j) &datarun.spikes{i,1} < timedur), nbins);
        else
            %[nhist(j,:) xouthist] = hist(datarun.spikes{i,1}(datarun.spikes{i,1} >= datarun.stimulus.triggers(1,j) &datarun.spikes{i,1} < datarun.stimulus.triggers(1,j+1)), nbins);
            [nhist(j,:) xouthist] = hist(datarun.spikes{i,1}(datarun.spikes{i,1} >= datarun.stimulus.triggers(1,j) &datarun.spikes{i,1} < (datarun.stimulus.triggers(1,j)+timestop)), nbins); %8s stimulus run
        end
        
    end
    NumSpikesCell(i,:) = max(grpstats(nhist,datarun.stimulus.trial_list,'mean')');
end
