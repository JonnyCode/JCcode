function ForIgorNew = PopDataConcatAndStat(ForIgor,StatString) ;

% this function will take a structure, "ForIgor" and look for fields with
% similar names from different cells, concatinate them, and find the
% requested population statistics

% JC 7/2/11
ForIgorNew = ForIgor ;

% concatinate all population data into "All" fields
if sum(strcmpi({'all','average','sem'},StatString))>0 ; % if any concatinating is needed
    ForIgorFields = fieldnames(ForIgorNew) ; % field names in structure
    for a=1:length(ForIgorFields) ; % for every field
        Breaki = strfind(ForIgorFields{a},'cell') ; % look for the string 'cell'
        if ~isempty(Breaki) ; % if such a string exists
            AllField = [ForIgorFields{a}(1:Breaki-1),'All'] ; % create a field with the name but replace the word "cell" and everything after with the word "All"
            if sum(strcmp(AllField,ForIgorFields))<1 ; % if no fields with that name exist 
                ForIgorNew.(AllField) = ForIgorNew.(ForIgorFields{a}) ; % make that field the same as the cell field
                ForIgorFields = fieldnames(ForIgorNew) ; % put that fieldname in
            else % if it does exist
                if size(ForIgorNew.(AllField),2) == size(ForIgorNew.(ForIgorFields{a}),2) ; % if they are concatinable
                    ForIgorNew.(AllField) = [ForIgorNew.(AllField); ForIgorNew.(ForIgorFields{a})] ; % concatinate the next field on
                else
                    ForIgorNew.(AllField) = 'THESE ARE DIFFERENT LENGTHS!' ;
                end
            end
        end
    end
end

% get average of population data and remove "All" fields
if strcmpi('average',StatString) ; % if average is requested
    ForIgorFields = fieldnames(ForIgorNew) ; % field names in structure
    for a=1:length(ForIgorFields) ; % for every field
        Breaki = strfind(ForIgorFields{a},'All') ; % look for the string 'All'
        if ~isempty(Breaki) && Breaki+3>length(ForIgorFields{a}); % if such a string exists and it is not followed by another string (i.e. it is not "AllSem", etc)
            MeanField = [ForIgorFields{a},'Mean'] ; % create a field with the name and add the word Mean
            ForIgorNew.(MeanField) = nanmean(ForIgorNew.(ForIgorFields{a}),1) ; % take the mean of the data and put into the new field 
            ForIgorNew = rmfield(ForIgorNew,ForIgorFields{a}) ; % remove the All field
        end
    end
end
 
% get sem of population data and remove "All" fields
if strcmpi('sem',StatString) ; % if average is requested
    ForIgorFields = fieldnames(ForIgorNew) ; % field names in structure
    for a=1:length(ForIgorFields) ; % for every field
        Breaki = strfind(ForIgorFields{a},'All') ; % look for the string 'All'
        if ~isempty(Breaki) && Breaki+3>length(ForIgorFields{a}); % if such a string exists and it is not followed by another string (i.e. it is not "AllSem", etc)
            SemField = [ForIgorFields{a},'Sem'] ; % create a field with the name and add the word Mean
            ForIgorNew.(SemField) = nanstd(ForIgorNew.(ForIgorFields{a}),[],1)./sqrt(size(ForIgorNew.(ForIgorFields{a}),1)) ; % take the mean of the data and put into the new field 
            ForIgorNew = rmfield(ForIgorNew,ForIgorFields{a}) ; % remove the All field
        end
    end
end

% get fraction of change in pairwise spike correlation and covariation +ei vs. -ei (e.g. spikeCorr-ei/spikeCorr+ei)
if strcmpi('FracChange',StatString) ;
    ForIgorFields = fieldnames(ForIgorNew) ;
    for a=1:length(ForIgorFields) ; % for every field
        potential = strfind(ForIgorFields{a},'snCorrCoefsMcPp') ; % look for this string
        if ~isempty(potential) ; % if it exists
            Breaki = strfind(ForIgorFields{a},'AllMean') ; % look for the string 'AllMean'
                if ~isempty(Breaki) && Breaki+7>length(ForIgorFields{a}); % if such a string exists and it is not followed by another string (i.e. it is not "AllMeanHiMom", etc
                    
                    denominator_string = ['snCorrCoefsPcPp',ForIgorFields{a}(length('snCorrCoefsPcPp')+1:end)]  ; % then find the denominator    
                    for b=1:length(ForIgorFields) ; % for every field
                        if strcmp(ForIgorFields{b},denominator_string) == 1 ; % if the string matches
                            identifier = ['snCorrCoefs',ForIgorFields{a}(length('snCorrCoefsPcPp')+1:end),'PC'] ; % make new field name
                            ForIgorNew.(identifier) = (ForIgor.(ForIgorFields{a})-ForIgor.(ForIgorFields{b}))/ForIgor.(ForIgorFields{a}) *100 ; % calculate new field
                        end
                    end
                end
        end
    end
    
    ForIgorFields = fieldnames(ForIgorNew) ;
    for a=1:length(ForIgorFields) ; % for every field
        potential = strfind(ForIgorFields{a},'snCovsMcPp') ; % look for this string
        if ~isempty(potential) ; % if it exists
            Breaki = strfind(ForIgorFields{a},'AllMean') ; % look for the string 'AllMean'
                if ~isempty(Breaki) && Breaki+7>length(ForIgorFields{a}); % if such a string exists and it is not followed by another string (i.e. it is not "AllMeanHiMom", etc
                    
                    denominator_string = ['snCovsPcPp',ForIgorFields{a}(length('snCovsPcPp')+1:end)]  ; % then find the denominator    
                    for b=1:length(ForIgorFields) ; % for every field
                        if strcmp(ForIgorFields{b},denominator_string) == 1 ; % if the string matches
                            identifier = ['snCovs',ForIgorFields{a}(length('snCovsPcPp')+1:end),'PC'] ; % make new field name
                            ForIgorNew.(identifier) = (ForIgor.(ForIgorFields{a})-ForIgor.(ForIgorFields{b}))/ForIgor.(ForIgorFields{a}) *100 ; % calculate new field
                        end
                    end
                end
        end
    end
end
                        
        

        
    





    