function [ReturnSearchPara,ReturnedField,FindFlag] = FindFamPara(Condition, SearchPara)
% ReturnSearchPara = FindFamCue(Condition, SearchPara)
%
% Given a FamilyCondition (ie. CellInfo.FamilyCondition(2)) and a parameter in 
% FamilyCueGuide or FamilyStepGuide, this function will return the value of the 
% parameter in the condition, as well as which field and row it was found in.  
%  Created:  09/06/01  MKMK
%		

FindFlag = 0;
CheckFind = 1;
NumSearchCrit = length(Condition.FamilyCueGuide);
for cnt = 1: NumSearchCrit
	if strcmp(Condition.FamilyCueGuide(cnt), SearchPara)
		% In case there is more than one cue, use FindFlag to determine which 
		% row corresponds to the cue we are trying to find.
		FindFlag = FindFlag + 1;
		CheckFind = 0;
		ReturnSearchPara = unique(Condition.FamilyCues(FindFlag,:));
		ReturnedField = 'FamilyCues';
		break
	else
		FindFlag = FindFlag + 1;
	end
end
if CheckFind == 1
	% FamilyStep must be only one parameter
	if strcmp(Condition.FamilyStepGuide, SearchPara)
		CheckFind = 0;
		ReturnSearchPara = unique(Condition.FamilyStep);
		ReturnedField = 'FamilyStep';
	end
end
if CheckFind == 1
	fprintf(1,'The search parameter was not found in FamilyCueGuide or FamilyStepGuide \n');
	ReturnSearchPara = [];	
end
if FindFlag == 0
	FindFlag = 1;
end
