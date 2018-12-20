function FieldNamesCell = FieldNameCellMaker(FrontEndString,Trials)

%This function will make a cell array of strings that specify field names
%for use with "DataStruct2Mat.m".  It simply concatinates the "Trials" vector with
%the "FrontEndString".

%JC 1/3/12

for a = 1:length(Trials) ; % for each Trial
    if Trials(a)<10 ;
        TrialString = ['0',num2str(Trials(a))] ;
    else
        TrialString = num2str(Trials(a)) ;
    end
    
    FieldNamesCell{a} = [FrontEndString,TrialString] ;
end
