% Function to calculate the mean of an entity.

function calcMean(obj)

    % Call the appropriate function based on the input's type.
    if isa(obj, 'Epoch')
        calcEpochMean(obj);
    elseif isa(obj, 'EpochGroup')
        calcGroupMean(obj);
    elseif isa(obj, 'Population')
        calcPopMean(obj);
    end
end


function calcEpochMean(epoch)
    epoch.results.mean = mean(epoch.response);
end


function calcGroupMean(group)
    for i = 1:length(group.epochs)
        
        epoch = group.epochs(i);
        
        % If this epoch mean wasn't already calculated we need to calculate it 
        % now.
        if ~isfield(epoch.results, 'mean')
            calcEpochMean(epoch);
        end
        
        epochMeans(i) = epoch.results.mean;
    end
    
    % Store calculated mean in results struct.
    group.results.mean = sum(epochMeans) / length(epochMeans);
end


function calcPopMean(pop)
    for i = 1:length(pop.epochGroups)
        
        group = pop.epochGroups(i);
        
        % If this epoch groups mean wasn't already calculated we need to
        % calculate it now.
        if ~isfield(group.results, 'mean')
            calcGroupMean(group);
        end
        
        groupMeans(i) = group.results.mean;
    end

    % Store calculated mean in results struct.
    pop.results.mean = sum(groupMeans) / length(groupMeans);
end
