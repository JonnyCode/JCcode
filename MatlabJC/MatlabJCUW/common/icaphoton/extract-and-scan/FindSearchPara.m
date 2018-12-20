function ReturnSearchPara = FindSearchPara(Condition, SearchPara)
% ReturnSearchPara.m
%
% ReturnSearchPara = FindSearchPara(Condition, SearchPara)
%
% Given an EpochCondition and a parameter in SearchCrit, this function
% will find the value of the parameter in the condition.
%  Created:  09/06/01  GDF
%

FindFlag = 0;
CheckFind = 1;
if isfield(Condition,'SearchCrit');
    NumSearchCrit = length(Condition.SearchCrit);
    for cnt = 1: NumSearchCrit
        if strcmp(Condition.SearchCrit(cnt), SearchPara);
            %disp('found')
            FindFlag = FindFlag + 1;
            CheckFind = 0;
            break
        else
            FindFlag = FindFlag + 1;
        end
    end

    ReturnSearchPara = Condition.SearchPara(FindFlag);

    if CheckFind == 1
        fprintf(1,'The search parameter was not found in SearchCrit \n');
        ReturnSearchPara = [];
    end
else
    disp('SearchCrit is not a field in this Epoch or FamilyCondition')
    ReturnedSearchPara = [];
end


