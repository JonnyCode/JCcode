function[NumSpikesCell, StimComb] = get_meanISI(datarun)


NumSpikesCell = zeros(length(datarun.cell_ids), length(datarun.stimulus.combinations));
StimComb = zeros(length(datarun.stimulus.combinations), length(fieldnames(datarun.stimulus.params))); 
%timedur = datarun.stimulus.triggers(1,length(datarun.stimulus.triggers)) + mean(diff(datarun.stimulus.triggers));
timedur = datarun.stimulus.triggers(1,length(datarun.stimulus.triggers)) + 8;%8s stimulus run
timestop = 8; %8s stimulus run

StimComb(:,1) = [datarun.stimulus.combinations.SPATIAL_PERIOD];
StimComb(:,2) = [datarun.stimulus.combinations.TEMPORAL_PERIOD];
StimComb(:,3) = [datarun.stimulus.combinations.DIRECTION];


for i = 1:length(datarun.cell_ids)
    ISI =[];
    ISInew = zeros(length(datarun.stimulus.combinations), 1);
    for j = 1: length(datarun.stimulus.triggers)
        isis=[];
        if(j == length(datarun.stimulus.triggers))
             isis = diff(datarun.spikes{i,1}(datarun.spikes{i,1} >= datarun.stimulus.triggers(1,j) &datarun.spikes{i,1} < timedur))';
             ISI(j, 1:numel(isis)) = isis;
        else
            %isis = diff(datarun.spikes{i,1}(datarun.spikes{i,1} >= datarun.stimulus.triggers(1,j) &datarun.spikes{i,1} < datarun.stimulus.triggers(1,j+1)))';
            isis = diff(datarun.spikes{i,1}(datarun.spikes{i,1} >= datarun.stimulus.triggers(1,j) &datarun.spikes{i,1} < datarun.stimulus.triggers(1,j)+timestop))';
            ISI(j, 1:numel(isis)) = isis;
        end
        for k = 1:length(datarun.stimulus.combinations)
            if(k == datarun.stimulus.trial_list(1,j))
                   ISInew(k,1:(numel(ISInew(k,:))+numel(isis))) = horzcat(isis, ISInew(k,:));
            end
        end
        
    end
    ISInew = ISInew';
    for k = 1:length(datarun.stimulus.combinations)
        NumSpikesCell(i,k) = mean(nonzeros(ISInew(:,k)));
    end

end