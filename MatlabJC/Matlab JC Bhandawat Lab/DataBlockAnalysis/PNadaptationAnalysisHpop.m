function ForIgorNew = PNadaptationAnalysisHpop(ForIgor) 

% this function will analyze data across multiple cell outputs from
% PNadaptationAnalysisH
% JC 1/11/12

ForIgorNew = ForIgor ;

% concatinate all population data into "All" fields

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


% divide up the spike rate data by pairwise within cells

ForIgorFields = fieldnames(ForIgorNew) ; % field names in structure
for a=1:length(ForIgorFields) ; % for every field
    Breaki = strfind(ForIgorFields{a},'All') ; % look for the string 'All'
    if ~isempty(Breaki) && Breaki+3>length(ForIgorFields{a}); % if such a string exists and it is not followed by another string (i.e. it is not "AllSem", etc)
        if strncmp(ForIgorFields{a},'SRminBgMeanBg',13) | strncmp(ForIgorFields{a},'SRminBgSemBg',13) ; % if the first 13 characters of the string are ....
            for b=1:length(ForIgor.LogConcentrationAllMean) ; % for each pulse strength used
                identifier = [ForIgorFields{a},'DbP',num2str(abs(ForIgor.LogConcentrationAllMean(b)))] ; % create a field with the name DbP (divy by pulse)
                ForIgorNew.(identifier) = ForIgorNew.(ForIgorFields{a})(:,b) ; % grab column that corresponds
            end
        end
    end
end

                    
                
 
% fit with sigmoid (fiting with sigmoid is prob not good idea until
% saturation points at both low and high end are found

    
    




